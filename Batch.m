classdef Batch < Reactor
    % This subclass defines a Batch Reactor
    % Features:
    %   - Fix the parameters particular to a Batch Reactor
    %   - Compute the output stream of individual Batch given its V and
    %   reaction time
    %   - Compute the reaction time required in a Batch Reactor to get a
    %   desired product
    % =========================================================================
    % Isabela Fons Moreno-Palancas
    % Created: March 14, 2020. Last update: April 1, 2020
    % Corrected: October 1, 2026 (T-101)
    % =========================================================================
    
    properties
        timeBatch
        timeBatch_Units = 's' ;
    end
    
    properties (Hidden = true) % This property is not displayed on the property list
        heatArray = [] ; % Stores the value of Q during the reaction time >> Useful to compute OPEX
    end
    
    methods
        function R = set.timeBatch(R,timeBatch)
            R.timeBatch = timeBatch ;
        end
        
        function R = set.timeBatch_Units(R,timeBatch_Units)
            R.timeBatch_Units = timeBatch_Units ;
        end
        
        function [Product,R] = compute_output(R,Feed,RS)
            %% @compute_output computes the evolution of the composition P and T with time
            % Features:
            %   - Momentum balance is neglected for liquid fase and computed for gas phase
            % considering the ideal gas law.Fix the parameters of thereactive system by hand or using Hysys
            % =========================================================================
            % Isabela Fons Moreno-Palancas
            % Last update: March 18, 2020
            % =========================================================================%
            
            %% Definition of the initial conditions to solve the ODE system
            moles_inlet = Feed.molarFlow/(1+R.bypassRatio) ;
            InitialConditions = [moles_inlet Feed.T Feed.P] ;
            
            % Call to the ODE function
            [t,y] = ode45(@ode_outputBatch , [0 R.timeBatch] , InitialConditions ) ;
            
            % Computation of the OUTPUT VARIABLES
            % Mass balance in the mixer
            moles_beforeMix = y(end,1:RS.nComponents) ;
            moles_out = moles_beforeMix + moles_inlet*R.bypassRatio ;
            % Momentum balance in the mixer
            P_out = y(end,RS.nComponents+2) ;
            % Energy balane in the mixer
            T_beforeMix = y(end,(RS.nComponents+1)) ;
            
            componentCp_beforeMix = RS.compute_HeatCapacity(T_beforeMix,P_out) ;
            componentCp_bypass    = RS.compute_HeatCapacity(Feed.T,Feed.P) ;
            
            T_out = (componentCp_beforeMix*moles_beforeMix'*T_beforeMix + componentCp_bypass*(moles_inlet'*R.bypassRatio)*Feed.T)/(componentCp_beforeMix*moles_beforeMix'+componentCp_bypass*(moles_inlet*R.bypassRatio)') ;
            
            % Definition of the product stream
            Product = Stream ;
            Product.molarFlow = moles_out ;
            Product.molarFlow_Units = Feed.molarFlow_Units ;
            Product.T = T_out ;
            Product.P = P_out ;
            Product.phase = Feed.phase ;
            Product.viscosity = Feed.viscosity ;
            if strcmp(Product.phase,'L')
                Product.density = Feed.density ;
            end
            
            %%  Plots
            if strcmp(R.activatePlots,'on') == 1
                
                BATCH_molarProfile = figure ;
                BATCH_temperatureProfile = figure ;
                BATCH_pressureProfile = figure ;
                
                figure(BATCH_molarProfile)
                plot(t,y(:,1:RS.nComponents))
                xlabel(string(['Time (',R.timeBatch_Units,')'])),ylabel(string(['Moles (', Feed.molarFlow_Units,')'])),title('Evolution of the composition')
                figure(BATCH_temperatureProfile)
                plot(t,y(:,RS.nComponents+1))
                xlabel(string(['Time (',R.timeBatch_Units,')'])),ylabel('T (K)'),title('Temperature Profile')
                figure(BATCH_pressureProfile)
                plot(t,y(:,RS.nComponents+2)/1000)
                xlabel(string(['Time (',R.timeBatch_Units,')'])),ylabel('P (kPa)'),title('Pressure Profile')
            end
            %% Compute the derivatives of the ODE system
            function dydt = ode_outputBatch(t,y)
                
                moles = y(1:RS.nComponents) ;
                T = y(RS.nComponents+1) ;
                P = y(RS.nComponents+2) ;
                
                Rg  = 8.314 ; %Pa·m^3/mol/K
                
                % Rate of reaction
                concentration = moles/R.V ; %mol/m^3
                RS = RS.computeRate(concentration,T) ;    
                constant_WtoV = (1-R.porosityCatalyst)*R.densityCatalyst ; % constant_WtoV is a conversion factor to change from mol/(time·kg cat) to mol/(time·m^3 reactor)
                r_i = constant_WtoV*RS.r_i ; %[1 x nReactions]
                
                % Heat exchange
                componentCp = RS.compute_HeatCapacity(T,P) ;
                DH = RS.DHref + (componentCp*RS.stochiometricMatrix')*(T-RS.Tref) ; %J/mol [1xnReactions]
                if strcmp(R.heatMode,'Isothermal') == 1
                    Q = R.V*(r_i*DH') ;
                elseif strcmp(R.heatMode,'Adiabatic') == 1
                    Q = 0 ;
                elseif strcmp(R.heatMode,'Other') == 1
                    if isempty(R.outletUtilityTemperature)
                        meanTemperatureDifference = R.inletUtilityTemperature-T ;
                    else
                        deltaTIn = R.inletUtilityTemperature-T ;
                        deltaTOut = R.outletUtilityTemperature-T ;
                        % At a temperature cross (or equal end differences),
                        % use the arithmetic mean because the LMTD is singular.
                        if deltaTIn*deltaTOut <= 0 || deltaTIn == deltaTOut
                            meanTemperatureDifference = (deltaTIn+deltaTOut)/2 ;
                        else
                            meanTemperatureDifference = ...
                                (deltaTIn-deltaTOut)/log(deltaTIn/deltaTOut) ;
                        end
                    end
                    Q = R.U*R.heatTransferArea*meanTemperatureDifference ;
                end
                
                R.heatArray = [ R.heatArray ; Q ] ; 
                
                % Mass, energy and momentum balances
                dndt = R.V * r_i * RS.stochiometricMatrix ;
                dTdt = (Q - R.V*(r_i*DH'))/(componentCp*moles) ;
                if Feed.phase == 'L'
                    dPdt = 0 ;
                elseif Feed.phase == 'G'
                    dPdt = (Rg/R.V)*(T*sum(dndt) + sum(moles)*dTdt) ;
                end
                
                dydt(1:RS.nComponents) = dndt' ;
                dydt((RS.nComponents)+1) = dTdt ;
                dydt((RS.nComponents)+2) = dPdt ;
                dydt = dydt' ;
            end
        end
        
        function Y = compute_timeBatch(R,Feed,RS,desired_conversion,keyComponentIndex)
            %% @compute_timeBatch computes reaction time as a function of conversion
            % Features:
            %   - This function can only be applied to systems with one reaction
            %   - Momentum balance is neglected
            % =========================================================================
            % Isabela Fons Moreno-Palancas
            % Last update: March 27, 2020
            % =========================================================================%
            %%
            % Definition of the INITIAL CONDITIONS to solve the ODE system
            moles_inlet = Feed.molarFlow/(1+R.bypassRatio) ;
            InitialConditions = [0 Feed.T Feed.P] ;
            
            % Definition of the KEY COMPONENT
            if rank(RS.stochiometricMatrix) == 1
                alpha = RS.stochiometricMatrix(1,:);
            end
            moles_key_inlet = moles_inlet(keyComponentIndex) ;
            alpha_key = alpha(keyComponentIndex) ;
            
            %Call to the ODE function
            [X,y] = ode45(@ode_timeBatch , [0 desired_conversion] , InitialConditions) ;
            
            % Computation of the OUTPUT VARIABLES
            R.timeBatch = y(end,1);
            P_out = y(end,3);
            moles_beforeMix = moles_inlet + alpha*(moles_key_inlet*X(end)/alpha_key) ;
            moles_out = moles_beforeMix + moles_inlet*R.bypassRatio ;
            T_beforeMix = y(end,2);
            % Computation of the molar heat capacity for each component.
            componentCp_beforeMix = RS.compute_HeatCapacity(T_beforeMix,P_out) ;
            componentCp_bypass    = RS.compute_HeatCapacity(Feed.T,Feed.P) ;

            T_out = (componentCp_beforeMix*moles_beforeMix'*T_beforeMix + componentCp_bypass*(moles_inlet'*R.bypassRatio)*Feed.T)/(componentCp_beforeMix*moles_beforeMix'+componentCp_bypass*(moles_inlet*R.bypassRatio)') ;
            
            % Definition of the product stream
            Product = Stream ;
            Product.molarFlow = moles_out ;
            Product.molarFlow_Units = Feed.molarFlow_Units ;
            Product.T = T_out ;
            Product.P = P_out ;
            Product.phase = Feed.phase ;
            Product.viscosity = Feed.viscosity ;
            if strcmp(Product.phase,'L')
                Product.density = Feed.density ;
            end
            
            Y = {Product,R} ;
            
            %% Compute the derivatives of the ODE system
            function dydX = ode_timeBatch(X,y)
                
                t = y(1) ;
                T = y(2) ;
                P = y(3) ;
                
                % Rate of reaction
                extent = moles_key_inlet * X / -alpha_key ;
                moles = moles_inlet + alpha*extent ;
                concentration = moles/R.V ;
                RS = RS.computeRate(concentration,T) ;
                r_i = RS.r_i ; %[1 x nReactions]
                if length(r_i) == 2
                    r_global = r_i(1) - r_i(2) ; % r_global = r_dir - r_inv
                else
                    r_global = r_i ;
                end
                % Mass balance
                constant_WtoV = (1-R.porosityCatalyst)*R.densityCatalyst ;
                dtdX = moles_key_inlet / (constant_WtoV * r_global * R.V * -alpha_key );
                
                % Energy balance
                if strcmp(R.heatMode,'Isothermal') == 1
                    dTdX = 0 ;
                else
                    if strcmp(R.heatMode,'Adiabatic') == 1
                        Q = 0 ;
                    elseif strcmp(R.heatMode,'Other') == 1
                        Q = R.U * R.heatTransferArea * (R.inletUtilityTemperature - T) ;
                    end
                    
                    componentCp = RS.compute_HeatCapacity(T,P) ;
                    DCp = alpha * componentCp' ;
                    DH = RS.DHref + DCp*(T-RS.Tref) ;
                    dTdt = (Q - DH * (moles_key_inlet/-alpha_key) * (1/dtdX))/(moles_inlet*componentCp' + DCp*X*moles_key_inlet/-alpha_key) ;
                    dTdX = dTdt*dtdX ;
                end
                
                %Momentum balance
                dPdX = 0 ;
                
                dydX = [dtdX ; dTdX ; dPdX ] ;
            end
        end
    end
end

