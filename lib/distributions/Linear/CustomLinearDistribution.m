classdef CustomLinearDistribution < AbstractLinearDistribution & CustomDistribution
    % Hypertoroidal distribution with custom pdf.
    
    methods
        function this = CustomLinearDistribution(f_, dim_)
            % Constructor, it is the user's responsibility to ensure that f is a valid
            % linear density and takes arguments of the same form as
            % .pdf, i.e., it needs to be vectorized.
            % 
            % Parameters:
            %   f_ (function handle)
            %       pdf of the distribution
            %   dim_ (scalar)
            %       dimension of the hypertorus
            this@CustomDistribution(f_, dim_);
        end  
    end
    
    methods (Static)
        function chd = fromDistribution(dist)
            % Creates a CustomLinearDistribution from some other distribution
            %
            % Parameters:
            %   dist (AbstractLinearDistribution)
            %       distribution to convert
            % Returns:
            %   chd (CustomLinearDistribution)
            %       CustomLinearDistribution with identical pdf
            arguments
                dist AbstractLinearDistribution
            end            
            chd = CustomLinearDistribution(@(xa)dist.pdf(xa), dist.dim);
        end
    end
end
