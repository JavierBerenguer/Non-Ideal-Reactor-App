classdef PFR < Reactor
    % This subclass defines a PFR
    % Features:
    %   - Fix the parameters particular to a PFR
    %   - Compute the output stream of individual PFR
    % =========================================================================
    % Isabela Fons Moreno-Palancas
    % Last update: April 1, 2020
    % Corrected: October 1, 2026 (T-101)
    % =========================================================================
    
    properties
        L
        diameterTubes = 0.1 ; % m
        nTubes = 1 ;
        pressureDropEqn = 'Pipe' ;
        particleDiameter
    end
    
    properties (Hidden = true) % This property is not displayed on the property list
        heatArray = [] ; % Stores the value of dQdL along the reactor >> Useful to compute OPEX
    end
    
    methods
          
        function R = set.L(R,L)
            R.L = L ;
        end
        
        function R = set.diameterTubes(R,diameterTubes)
            R.diameterTubes = diameterTubes ;
        end
        
        function R = set.nTubes(R,nTubes)
            R.nTubes = nTubes ;
        end
        
        function R = set.pressureDropEqn(R,pressureDropEqn)
            pressureDropEqn_options = {'Pipe','Ergun'};
            switch pressureDropEqn
                case pressureDropEqn_options
                    R.pressureDropEqn = pressureDropEqn ;
                otherwise
                    error('Invalid pressure drop equation name') ;
            end
        end
        
        function [Product,R] = compute_output(R,Feed,RS)
            %% @compute_output computes the composition, T and P of a PFR
            %
            % =========================================================================
            % Isabela Fons Moreno-Palancas
            % Last update: April 16, 2020
            % =========================================================================%
            %% Length takes precedence; otherwise derive it from reactor volume.
            if ~isempty(R.L)
                R.V = R.nTubes*pi*(R.diameterTubes/2)^2*R.L ;
            else
                R.L = R.V/(R.nTubes*pi*(R.diameterTubes/2)^2) ;
            end
            
            % Operations to compute the initial conditions
            moles_inlet = Feed.molarFlow/(1+R.bypassRatio) ;
            moles_inlet_tube = moles_inlet/R.nTubes ;
            
            InitialConditions = [moles_inlet_tube Feed.T Feed.P R.inletUtilityTemperature] ;
            absoluteTolerance = max(1e-10*abs(InitialConditions),1e-12) ;
            options = odeset('RelTol',1e-8,'AbsTol',absoluteTolerance) ;
            [L,y] = ode45(@odePFR,[0 R.L],InitialConditions,options) ;

            % computeCost integrates heatArray on a uniform 0..L mesh.
            heatMesh = linspace(0,R.L,201)' ;
            heatStates = interp1(L,y,heatMesh,'pchip') ;
            R.heatArray = zeros(size(heatMesh)) ;
            for iHeat = 1:numel(heatMesh)
                R.heatArray(iHeat) = heatPerLength(heatStates(iHeat,:)) ;
            end
            
            % Mass balance in the mixer
            moles_beforeMix = y(end,1:RS.nComponents)*R.nTubes ;
            moles_out = moles_beforeMix + moles_inlet*R.bypassRatio ;
            % Momentum balance in the mixer
            P_out = y(end,RS.nComponents+2);
            % Energy balance in the mixer
            T_beforeMix = y(end,(RS.nComponents+1)) ; % Temperature of the stream leaving the reactor before entering the mixer
            componentCp_beforeMix = RS.compute_HeatCapacity(T_beforeMix,P_out) ;
            componentCp_bypass    = RS.compute_HeatCapacity(Feed.T,Feed.P) ;
            T_out = (componentCp_beforeMix*moles_beforeMix'*T_beforeMix + componentCp_bypass*(moles_inlet'*R.bypassRatio)*Feed.T)/(componentCp_beforeMix*moles_beforeMix'+ componentCp_bypass*(moles_inlet*R.bypassRatio)') ;
            
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
            
            %%  Plots
            if strcmp(R.activatePlots,'on') == 1
                
                PFR_molarProfile = figure ;
                PFR_temperatureProfile = figure ;
                PFR_pressureProfile = figure ;
                
                figure(PFR_molarProfile)
                plot(L,y(:,1:RS.nComponents))
                xlabel('L (m)'),ylabel(string(['Molar Flow (', Feed.molarFlow_Units,')'])),title('Molar Flows Profile')
                if ~isempty(RS.componentNames)
                    legend(RS.componentNames)
                end
                figure(PFR_temperatureProfile)
                plot(L,y(:,RS.nComponents+1))
                xlabel('L (m)'),ylabel('T (K)'),title('Temperature Profile')
                figure(PFR_pressureProfile)
                plot(L,y(:,RS.nComponents+2)/1000)
                xlabel('L (m)'),ylabel('P (kPa)'),title('Pressure Profile')
                
            end
            %% Compute the derivatives of the ODE system
            
            function dydL = odePFR(~,y)
                
                moles = y(1:RS.nComponents);
                T = y(RS.nComponents+1);
                P = y(RS.nComponents+2);
                utilityT = y(RS.nComponents+3);
                
                Rg = 8.314; %J/mol/K
                
                %Rate of reaction
                if strcmp(Feed.phase, 'L')
                    Qv = Feed.volumetricFlow ;
                elseif strcmp(Feed.phase, 'G')
                    Qv = sum(moles) * Rg * T/P ; %m^3/s
                end
                concentration = moles/Qv ; %mol/m^3
                RS = RS.computeRate(concentration,T) ;
                constant_WtoV = (1-R.porosityCatalyst)*R.densityCatalyst ; % constant_WtoV is a conversion factor to change from mol/(time·kg cat) to mol/(time·m^3 reactor)
                r_i = constant_WtoV*RS.r_i ; %[1 x nReactions]
                
                %Mass balance
                crossSectionalArea = (R.V/R.L)/R.nTubes ; %m^2
                dndL = crossSectionalArea*r_i*RS.stochiometricMatrix ;
                
                %Energy balance
                componentCp = RS.compute_HeatCapacity(T,P) ;
                DH = RS.DHref + (componentCp*RS.stochiometricMatrix')*(T-RS.Tref) ; %J/mol [1xnReactions]
                
                if strcmp(R.heatMode,'Isothermal') == 1
                    dQdL = crossSectionalArea*(r_i*DH') ;
                else
                    % Heat exchange < Energy Balance
                    if strcmp(R.heatMode,'Adiabatic') == 1
                        dQdL = 0 ;
                    elseif strcmp(R.heatMode,'Other') == 1
                        % dQdL = 2* R.nTubes * R.U * pi * R.diameterTubes *(utilityT-T); %No hace falta multiplicar por 2 porque se está usando el diámetro de los tubos y no el radio
                        dQdL = R.nTubes * R.U * pi * R.diameterTubes *(utilityT-T);
                    end
                end
                
                dTdL = (dQdL - crossSectionalArea*(r_i*DH'))/(componentCp*moles) ;
                
                % Momentum balance
                if strcmp(R.pressureMode,'Constant') == 1
                    dPdL = 0 ;
                else
                    if strcmp(Feed.phase,'L')
                        density = Feed.density ;
                    elseif strcmp(Feed.phase,'G')
                        density = (RS.componentMw*moles/Qv)/1000 ; %kg/m^3
                    end
                    v = Qv/crossSectionalArea ; %m/s
                    
                    if strcmp(R.pressureDropEqn,'Pipe')
                        % Pipe reactor
                        Re = density*v*R.diameterTubes/Feed.viscosity ;
                        f = 0.046*Re^(-0.2) ;
                        if Re > 4000
                            a = 0.96 ;
                        else
                            a = 0.5 ;
                        end
                        dPdL = (2*f*Qv/R.diameterTubes+ Rg/a/P*(dTdL*sum(moles)+T*sum(dndL)))/(Rg*T/(a*P^2)*sum(moles)-(crossSectionalArea^2/(RS.componentMw*moles/1000))) ; %Se divide entre 1000 porque Mw viene en g/mol
                    elseif strcmp(R.pressureDropEqn,'Ergun')
                        % Bed reactor. Fogler Page 178.
                        G = density*v ; % kg/m^2/s
                        dPdL = -G/density/R.particleDiameter * (1-R.porosityCatalyst)/(R.porosityCatalyst^3) * ( 150*(1-R.porosityCatalyst)*Feed.viscosity/R.particleDiameter + 1.75*G ) ;
                    end
                    
                end
                
                if isempty(R.outletUtilityTemperature)
                    dutilityTdL = 0 ;
                else
                    dutilityTdL = ...
                        (R.outletUtilityTemperature-R.inletUtilityTemperature)/R.L ;
                end
                
                dydL(1:RS.nComponents) = dndL' ;
                dydL((RS.nComponents)+1) = dTdL ;
                dydL((RS.nComponents)+2) = dPdL ;
                dydL((RS.nComponents)+3) = dutilityTdL ;
                dydL = dydL' ;
                
            end

            function dQdL = heatPerLength(state)
                moles = state(1:RS.nComponents) ;
                T = state(RS.nComponents+1) ;
                P = state(RS.nComponents+2) ;
                utilityT = state(RS.nComponents+3) ;
                crossSectionalArea = (R.V/R.L)/R.nTubes ; % m^2

                if strcmp(Feed.phase,'L')
                    Qv = Feed.volumetricFlow ;
                else
                    Qv = sum(moles)*8.314*T/P ; % m^3/s
                end
                RS = RS.computeRate(moles/Qv,T) ;
                constant_WtoV = (1-R.porosityCatalyst)*R.densityCatalyst ;
                r_i = constant_WtoV*RS.r_i ;
                componentCp = RS.compute_HeatCapacity(T,P) ;
                DH = RS.DHref + (componentCp*RS.stochiometricMatrix') * ...
                    (T-RS.Tref) ; % J/mol

                if strcmp(R.heatMode,'Isothermal')
                    dQdL = crossSectionalArea*(r_i*DH') ;
                elseif strcmp(R.heatMode,'Adiabatic')
                    dQdL = 0 ;
                else
                    dQdL = R.nTubes*R.U*pi*R.diameterTubes*(utilityT-T) ;
                end
            end
        end
    end
end

