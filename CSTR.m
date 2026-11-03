classdef CSTR < Reactor
    % This subclass defines a CSTR
    % Features:
    %   - Fix the parameters particular to a CSTR
    %   - Compute the output stream of individual CSTRs
    % =========================================================================
    % Isabela Fons Moreno-Palancas
    % Created: March 14, 2020. Last update: April 20, 2020
    % Corrected: November 2026 (energy balance sign, heat exchange mode,
    % gas volumetric flow, solver tolerances and residual scaling)
    % =========================================================================
    properties (Hidden = true) % This property is not displayed on the property list
        heatFlux % Stores the value of Q >> Useful to compute OPEX
    end
    
    methods
        
        function [Product,R] = compute_output(R,Feed,RS)
            %% @compute_output computes the composition, T and P of a CSTR
            % Features:
            %   - Pressure of the feed and product streams is assumed to be equal
            % =========================================================================
            % Isabela Fons Moreno-Palancas
            % Last update: March 27, 2020
            % =========================================================================%
            
            %%
            Guess = [Feed.molarFlow , Feed.T, Feed.P] ;
            % Dimensionless residuals so that the material, energy,
            % temperature and pressure equations have comparable weight
            % in fsolve. The physical balances are unchanged.
            convergenceTolerance = 1e-9 ;
            moles_inlet = Feed.molarFlow/(1+R.bypassRatio) ; % mol/s
            flowScale = max(max(abs(moles_inlet)),1e-8) ; % mol/s
            materialScale = max(abs(moles_inlet),flowScale) ; % mol/s
            inletCp = RS.compute_HeatCapacity(Feed.T,Feed.P) ; % J/(mol*K)
            energyScale = max(sum(abs(moles_inlet).*abs(inletCp))* ...
                max(abs(Feed.T),1),1) ; % W
            temperatureScale = max(abs(Feed.T),1) ; % K
            pressureScale = max(abs(Feed.P),1) ; % Pa

            typicalX = abs(Guess) ;
            typicalX(typicalX == 0) = flowScale ;
            options = optimoptions('fsolve','Display','none', ...
                'FunctionTolerance',1e-12,'StepTolerance',1e-12,'OptimalityTolerance',1e-12, ...
                'TypicalX',typicalX(:)) ; % column: fsolve scales internally by columns
            [y,~,exitflag] = fsolve(@fsolveCSTR,Guess,options) ;

            % Exact specifications should not retain roundoff introduced by
            % the scaling. Recheck all balances at the snapped state.
            if strcmp(R.heatMode,'Isothermal')
                y(RS.nComponents+1) = Feed.T ;
            end
            y(RS.nComponents+2) = Feed.P ;
            residual = fsolveCSTR(y) ;
            maximumScaledResidual = norm(residual,inf) ;
            if ~isfinite(maximumScaledResidual) || ...
                    maximumScaledResidual > convergenceTolerance
                warning('CSTR:notConverged', ...
                    ['fsolve did not converge (exitflag %d, maximum scaled ' ...
                    'residual %.3e).'],exitflag,maximumScaledResidual) ;
            end

            if strcmp(R.heatMode,'Adiabatic')
                R.heatFlux = 0 ;
            elseif strcmp(R.heatMode,'Other')
                R.heatFlux = computeHeatFlux(y(RS.nComponents+1)) ;
            else
                [r_i,DH,componentCp] = reactionProperties(y) ;
                moles_inlet = Feed.molarFlow/(1+R.bypassRatio) ;
                R.heatFlux = moles_inlet*componentCp' * ...
                    (y(RS.nComponents+1)-Feed.T) + R.V*(r_i*DH') ;
            end
            
            % Mass balance in the mixer
            moles_beforeMix = y(1:RS.nComponents) ;
            moles_bypass = Feed.molarFlow*R.bypassRatio/(1+R.bypassRatio) ;
            moles_out = moles_beforeMix + moles_bypass ;
            % Momentum balance in the mixer
            P_out = y(RS.nComponents+2) ;
            % Energy balance in the mixer
            T_beforeMix = y(RS.nComponents+1) ;
            componentCp_beforeMix = RS.compute_HeatCapacity(T_beforeMix,P_out) ;
            componentCp_bypass    = RS.compute_HeatCapacity(Feed.T,Feed.P) ;
            T_out = (componentCp_beforeMix*moles_beforeMix'*T_beforeMix + componentCp_bypass*moles_bypass'*Feed.T)/(componentCp_beforeMix*moles_beforeMix'+componentCp_bypass*moles_bypass');
            
            % Definition of the product stream
            Product = Stream ;
            Product.molarFlow = moles_out ;
            Product.molarFlow_Units = Feed.molarFlow_Units ;
            Product.T = T_out ;
            Product.P = P_out ;
            Product.phase = Feed.phase ;
            Product.viscosity = Feed.viscosity ;
            Product.volumetricFlow_Units = Feed.volumetricFlow_Units ;
            Product.concentration_Units = Feed.concentration_Units ;
            Product.volumetricFlow = [] ;
            Product.concentration = [] ;
            if strcmp(Product.phase,'L')
                Product.volumetricFlow = Feed.volumetricFlow ; 
                Product.density = Feed.density ;
            end
            %% Solve the NonLinear System
            function y = fsolveCSTR(x)
                
                moles = x(1:RS.nComponents) ;
                T = x(RS.nComponents+1) ;
                P = x(RS.nComponents+2) ;
                
                [r_i,DH,componentCp] = reactionProperties(x) ;
                r_j = r_i*RS.stochiometricMatrix ; %[1 x nComponents]
                               
                % Mass balance
                materialResidual = moles_inlet - moles + r_j*R.V ;
                y(1:RS.nComponents) = materialResidual./materialScale ;

                %Energy balance
                if strcmp(R.heatMode,'Isothermal') == 1
                    y(RS.nComponents+1) = (T-Feed.T)/temperatureScale ;
                else
                    if strcmp(R.heatMode,'Adiabatic') == 1
                        Q = 0;
                    elseif strcmp(R.heatMode,'Other') == 1
                        Q = computeHeatFlux(T) ;
                    end

                    % Table 1 of the ReactorApp article: sensible heat plus
                    % V*sum(r_i*DH_i), minus heat entering the reactor.
                    energyResidual = moles_inlet*componentCp' * ...
                        (T-Feed.T) + R.V*(r_i*DH') - Q ;
                    y(RS.nComponents+1) = energyResidual/energyScale ;
                end

                % Momentum balance
                y(RS.nComponents+2) = (P-Feed.P)/pressureScale ;
                
                y = y';
                
            end

            function [r_i,DH,componentCp] = reactionProperties(x)
                moles = x(1:RS.nComponents) ;
                T = x(RS.nComponents+1) ;
                P = x(RS.nComponents+2) ;

                if Feed.phase == 'L'
                    Qv = Feed.volumetricFlow ;
                elseif Feed.phase == 'G'
                    Qv = sum(moles)*8.314*T/P ; % m^3/s
                end
                concentration = moles./Qv ; % mol/m^3
                RS = RS.computeRate(concentration,T) ;
                constant_WtoV = (1-R.porosityCatalyst)*R.densityCatalyst ;
                r_i = constant_WtoV*RS.r_i ; % mol/(m^3*s)
                componentCp = RS.compute_HeatCapacity(T,P) ; % J/(mol*K)
                DH = RS.DHref + componentCp*RS.stochiometricMatrix' * ...
                    (T-RS.Tref) ; % J/mol
            end

            function Q = computeHeatFlux(T)
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
                Q = R.U*R.heatTransferArea*meanTemperatureDifference ; % W
            end
        end
        
    end
end

