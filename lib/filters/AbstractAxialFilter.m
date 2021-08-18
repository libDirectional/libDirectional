classdef (Abstract) AbstractAxialFilter < AbstractFilter
    % Abstract base class for filters on the hypersphere with antipodal symmetry
    
    properties
        compositionOperator % denoted by "\oplus" in our papers
        compositionOperatorDerivative
        % Dimension (dim=2 -> circle/complex numbers, dim=4 -> quaternions)
    end
    
    methods ( Access = protected )
        function setCompositionOperator(this)
            switch this.dim
                case 2
                    this.compositionOperator = @(a,b) complexMultiplication(a,b);
                    this.compositionOperatorDerivative = ...
                        @(a, b) [b(1) -b(2)  a(1) -a(2);
                                 b(2)  b(1)  a(2)  a(1)];
                case 4
                    this.compositionOperator = @(a,b) quaternionMultiplication(a,b);
                    this.compositionOperatorDerivative = ...
                        @(a, b) [b(1) -b(2) -b(3) -b(4)  a(1) -a(2) -a(3) -a(4);
                                 b(2)  b(1)  b(4) -b(3)  a(2)  a(1) -a(4)  a(3);
                                 b(3) -b(4)  b(1)  b(2)  a(3)  a(4)  a(1) -a(2);
                                 b(4)  b(3) -b(2)  b(1)  a(4) -a(3)  a(2)  a(1)];
                otherwise
                    error('invalid dimension');
            end
        end
    end
    
    methods
        function est = getPointEstimate(this)
            est = this.getEstimate().mode;
        end
    end
       
end
