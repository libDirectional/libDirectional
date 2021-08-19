classdef CustomHypercylindricalDistribution < CustomLinBoundedDistribution & AbstractHypercylindricalDistribution
    % Hypercylindrical distribution with custom pdf.
    
    methods
        function this = CustomHypercylindricalDistribution(f_, boundD_, linD_)
            arguments
                f_ (1,1) function_handle
                boundD_ (1,1) {mustBeInteger,mustBePositive}
                linD_ (1,1) {mustBeInteger,mustBePositive}
            end
            % Constructor, it is the user's responsibility to ensure that f is a valid
            % hypertoroidal density and takes arguments of the same form as
            % .pdf, i.e., it needs to be vectorized.
            % 
            % Parameters:
            %   f_ (function handle)
            %       pdf of the distribution
            %   dim_ (scalar)
            %       dimension of the hypertorus
            this@CustomLinBoundedDistribution(f_, boundD_, linD_);
        end
    end
    
    methods (Static)
        function chhd = fromDistribution(dist)
            arguments
                dist (1,1) AbstractHypercylindricalDistribution
            end
            chhd = CustomHypercylindricalDistribution(@(xa)dist.pdf(xa), dist.linD, dist.boundD);
        end
    end
end
