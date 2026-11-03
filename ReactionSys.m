classdef ReactionSys
% This class defines the reactive system. 
% Features:
%   - Fix the parameters of thereactive system by hand or using Hysys
%   - Includes the function that computes the rate of reaction for a given reactive system
% =========================================================================
% Isabela Fons Moreno-Palancas
% Last update: July 7, 2026
% Corrected: October 1, 2026 (T-101)
% =========================================================================
    properties
        componentNames
        componentFormula
        componentMw        
        componentHeatOfFormation
        componentCp
%       The property "componentCp" is likewise decomposed into a struct with the following properties
%             componentCp.option
%             componentCp.UserValues
%             componentCp.Average
%             componentCp.Function
%             componentCp.FunctionWithP
        
        stochiometricMatrix
        
        k0
        k0_units
        Ea 
        k0_denominator = 0 ;
        Ea_denominator = 0 ;
        partialOrders_denominator = 0 ; 
        
        Tref = 273.15 + 25 ; 
        DHref  
        r_i
        userDefinedKinetics 
        userDefinedKineticsExpressions = strings(0, 1)
        userDefinedKineticsSource = ''
        userDefinedKineticsFunctionName = ''
        userKineticsTimeUnit = 's'
        userKineticsConcentrationUnit = 'mol/m^3'
    end
    
    properties (Dependent)
        nReactions
        nComponents
        partialOrders 
    end
    
    methods

        function RS = setHysysProperties(RS,varargin)
            % DESCRIPTION: Constructor function of the class ReactionSys to obtain all the data available from Aspen Hysys
            % M. FILES REQUIRED:
            %   - call_DataBase.m
            % Hysys FILES REQUIRED (if it's not any of the listed below, specify it)
            %   - ComponentDataBase.xml or Styrene_reactor_v5.xml
            
            Component = call_DataBase(varargin) ;
            RS.componentFormula = {} ;
            RS.componentNames   = {} ;
            CpStruct.Function = {} ;
            CpStruct.FunctionWithP = {} ;
            for i = 1:size(Component,2)
                RS.componentFormula(i)     = mat2cell(Component(i).Formula,1) ;
                RS.componentNames(i)       = mat2cell(Component(i).Name,1) ;
                RS.componentMw(i)          = Component(i).MolecularWeightValue ;
                CpStruct.Average(i)        = Component(i).MolarHeatCapacityAverage ;
                CpStruct.Function{i}       = Component(i).MolarHeatCapacityFunction ;
                CpStruct.FunctionWithP{i}  = Component(i).MolarHeatCapacityFunction_withPressure ;
                RS.componentHeatOfFormation(i)   = Component(i).HeatOfFormation ;
            end
            RS.componentCp = CpStruct ;
        end

        function RS = set.stochiometricMatrix(RS,stochiometricMatrix)
            RS.stochiometricMatrix = stochiometricMatrix ;
        end
        
        function nComponents = get.nComponents(RS)
        % DESCRIPTION: In this function the user provides the stochiometric matrix and the number of components is automatically computed
            nComponents = size(RS.stochiometricMatrix,2) ;
        end
        
        function nReactions = get.nReactions(RS)
        % User provides the stochiometric matrix and the number of reactions is automatically computed
            nReactions = size(RS.stochiometricMatrix,1) ;
        end
        
        function partialOrders = get.partialOrders(RS) 
        % User provides the stochiometric matrix and partial orders assuming the reactions are elementary are automatically computed
            partialOrders = (RS.stochiometricMatrix<0).*abs(RS.stochiometricMatrix);
        end 

        function RS = set.DHref(RS,DHref)
            RS.DHref = DHref ;
        end
        
        function DHref = get.DHref(RS)
            % User provides the stochiometric matrix and enthalpy change of each reaction is automatically computed
            if isempty(RS.DHref) && isempty(RS.componentHeatOfFormation)
                DHref = ones(1,RS.nReactions) ;
            elseif length(RS.componentHeatOfFormation) == RS.nComponents
                DHref = RS.componentHeatOfFormation * RS.stochiometricMatrix' ;
            elseif isempty(RS.DHref)== 0
                DHref = RS.DHref ;
            end
        end

        function RS = set.componentMw(RS,componentMw)
            RS.componentMw = componentMw ;
        end
        
        function RS = set.componentCp(RS,componentCp)
        % In case it is not desired to conect with Hysys, the "set" function for that properties is enabled
            if isa(componentCp,'double')
                RS.componentCp.UserValues = componentCp ;
                RS.componentCp.option = 'User defined' ;
            elseif isa(componentCp,'struct')
                RS.componentCp = componentCp ;
            end
        end
        
        function componentCp = get.componentCp(RS)
            if isempty(RS.componentCp)
                componentCp.UserValues = ones(1,RS.nComponents) ;
                componentCp.option = 'User defined' ;
            else
                componentCp = RS.componentCp ;
            end
        end
        
        function componentCp = compute_HeatCapacity(RS,T,P)
        % Calculate the molar heat capacity for each component.
            componentCp = zeros(1,RS.nComponents) ;
            for i = 1:RS.nComponents
                if strcmp(RS.componentCp.option,'Average')
                    componentCp(i) = RS.componentCp.Average(i) ;
                elseif strcmp(RS.componentCp.option,'Cp = f(T)')
                    equation = RS.componentCp.Function{i} ;
                    componentCp(i) = equation(T) ;
                elseif strcmp(RS.componentCp.option,'Cp = f(T,P)')
                    equation = RS.componentCp.FunctionWithP{i} ;
                    componentCp(i) = equation(T,P) ;
                elseif strcmp(RS.componentCp.option,'User defined')
                    componentCp(i) = RS.componentCp.UserValues(i) ;
                end
            end
        end

        function RS = set.k0(RS,k0)
            RS.k0 = k0 ;
        end
        
        function RS = set.Ea(RS,Ea)
            RS.Ea = Ea ;
        end
        
        function Ea = get.Ea(RS)
            if isempty(RS.Ea)
                Ea = zeros(1,RS.nReactions) ;
            else
                Ea = RS.Ea ;
            end
        end
        
        function RS = set.k0_denominator(RS,k0_denominator)
            RS.k0_denominator = k0_denominator ;
        end
        
        function RS = set.Ea_denominator(RS,Ea_denominator)
            RS.Ea_denominator = Ea_denominator ;
        end

        function RS = set.partialOrders_denominator(RS,partialOrders_denominator)
            RS.partialOrders_denominator = partialOrders_denominator ;
        end
        
        function RS = set.Tref(RS,Tref)
            RS.Tref = Tref ;
        end
        
        function RS = set.k0_units(RS,k0_units)
            if ischar(k0_units)
                RS.k0_units = k0_units;
            end
        end
               
        function RS = computeRate(RS,concentration,T)

            if size(concentration) == [RS.nComponents,1]
                concentration = concentration' ;
            end

            if isempty(RS.userDefinedKinetics)
                % DESCRIPTION: function to compute the rate of reaction using the Langmuir-Hinshelwood equation (same as Aspen Hysys).
                %   Inputs : concentrations and temperature
                %   Output : r_i [ 1 x nReactions ]

                numerator = ones(1,RS.nReactions) ;
                for i = 1:RS.nReactions
                    k = RS.k0(i).*exp(-RS.Ea(i)/8.314/T) ;
                    numerator(i) = k * prod(concentration.^RS.partialOrders(i,:)) ;
                end

                denominator = 1 ;
                for ii = 1:length(RS.k0_denominator)
                    k_denominator = RS.k0_denominator(ii).*exp(-RS.Ea_denominator(ii)/8.314/T) ;
                    f = prod(concentration.^RS.partialOrders_denominator(ii,:)) ;
                    denominator = denominator + k_denominator*f ;
                end

                RS.r_i = numerator/denominator ;

            else
                % Specify rate equation as a function handle variable.
                % INPUT ARGUMENTS:
                %   - concentration as a vector of [1 x nComponents]
                %   - T
                % OUTPUT ARGUMENTS:
                %   - Rate of reaction as a vector of [1 x nReactions]
                RS.r_i = RS.userDefinedKinetics(concentration,T) ;

            end
        end

    end

    methods (Static)

        function RS = fromSimpleKinetics(type, params)
            % fromSimpleKinetics - Factory method to create ReactionSys from
            % the simplified kinetics used in the non-ideal reactor UI.
            %
            % All kinetics are isothermal (Ea = 0, T is ignored).
            % The ReactionSys is configured so that computeRate(concentration, T)
            % returns the correct r_i for each kinetics type.
            %
            % Usage:
            %   RS = ReactionSys.fromSimpleKinetics('1st_order', struct('k', 0.1))
            %   RS = ReactionSys.fromSimpleKinetics('2nd_order', struct('k', 0.01, 'CA0', 1000))
            %   RS = ReactionSys.fromSimpleKinetics('michaelis_menten', struct('a', 1, 'b', 0.5, 'CA0', 1000))
            %   RS = ReactionSys.fromSimpleKinetics('reversible', struct('kf', 0.1, 'kr', 0.05, 'CA0', 1000))
            %   RS = ReactionSys.fromSimpleKinetics('parallel', struct('k1', 0.1, 'n1', 2, 'k2', 0.05, 'n2', 1, 'CA0', 1000))
            %   RS = ReactionSys.fromSimpleKinetics('custom', struct('rateFcn', @(CA) 0.5*CA/(1+0.5*CA), 'CA0', 1000))
            %
            % Inputs:
            %   type   - string: '1st_order', '2nd_order', 'michaelis_menten',
            %            'reversible', 'parallel', 'custom'
            %   params - struct with kinetics parameters (see examples above)
            %
            % Output:
            %   RS - ReactionSys object configured for the selected kinetics.
            %        concentration vector convention:
            %          1st/2nd order: [CA]
            %          michaelis_menten: [CA]
            %          reversible: [CA, CB]
            %          parallel: [CA, CB, CC]
            %          custom: [CA]
            %
            % After creation, use:
            %   RS = RS.computeRate(concentration, T) ;
            %   rA = RS.r_i(1) ;  % rate of reaction 1 (positive = production rate per stoich)
            %
            % Internal units (SI):
            %   concentration: mol/m^3 | k(1st): 1/s | k(2nd): m^3/(mol*s)

            RS = ReactionSys ;

            switch type

                case '1st_order'
                    % A -> Products, -rA = k * CA
                    % stoichiometry: [-1]  (1 reaction, 1 component)
                    % partialOrders auto-computed from stoich: [1]
                    k = params.k ;
                    RS.stochiometricMatrix = [-1] ;
                    RS.k0 = k ;
                    RS.Ea = 0 ;

                case '2nd_order'
                    % A -> Products, -rA = k * CA^2
                    % Use userDefinedKinetics because auto partialOrders
                    % from stoich would give order 1, not 2.
                    k = params.k ;
                    RS.stochiometricMatrix = [-1] ;
                    RS.userDefinedKinetics = @(C, T) k * C(1)^2 ;

                case 'michaelis_menten'
                    % A -> Products, -rA = a * CA / (1 + b * CA)
                    a = params.a ;
                    b = params.b ;
                    RS.stochiometricMatrix = [-1] ;
                    RS.userDefinedKinetics = @(C, T) a * C(1) / (1 + b * C(1)) ;

                case 'reversible'
                    % A <-> B, -rA = kf*CA - kr*CB
                    % 2 components: [A, B], 1 reaction: [-1, +1]
                    kf = params.kf ;
                    kr = params.kr ;
                    RS.stochiometricMatrix = [-1, 1] ;
                    RS.userDefinedKinetics = @(C, T) kf * C(1) - kr * C(2) ;

                case 'parallel'
                    % A -> B  with r1 = k1 * CA^n1
                    % A -> C  with r2 = k2 * CA^n2
                    % 3 components: [A, B, C], 2 reactions
                    % stoich: [-1 1 0 ; -1 0 1]
                    k1 = params.k1 ; n1 = params.n1 ;
                    k2 = params.k2 ; n2 = params.n2 ;
                    RS.stochiometricMatrix = [-1 1 0 ; -1 0 1] ;
                    RS.userDefinedKinetics = @(C, T) ...
                        [k1 * max(C(1),0)^n1, k2 * max(C(1),0)^n2] ;

                case 'custom'
                    % -rA = rateFcn(CA), user-provided function handle
                    rateFcn = params.rateFcn ;
                    RS.stochiometricMatrix = [-1] ;
                    RS.userDefinedKinetics = @(C, T) rateFcn(C(1)) ;

                otherwise
                    error('ReactionSys:fromSimpleKinetics', ...
                        'Unknown kinetics type: "%s". Valid types: 1st_order, 2nd_order, michaelis_menten, reversible, parallel, custom.', type) ;
            end
        end

    end
end


