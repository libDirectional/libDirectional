classdef CustomHypertoroidalDistribution < AbstractHypertoroidalDistribution & CustomDistribution
    % Hypertoroidal distribution with custom pdf.
    
    methods
        function this = CustomHypertoroidalDistribution(f_,dim_)
            % Constructor, it is the user's responsibility to ensure that f is a valid
            % hypertoroidal density and takes arguments of the same form as
            % .pdf, i.e., it needs to be vectorized.
            % 
            % Parameters:
            %   f_ (function handle)
            %       pdf of the distribution
            %   dim_ (scalar)
            %       dimension of the hypertorus
            this@CustomDistribution(f_,dim_);
        end

        function ccd = toCustomCircular(this)
            % Convert to a custom circular distribution (only in 1D case)
            %
            % Returns:
            %   ccd (CustomCircularDistribution)
            %       CustomCircularDistribution with same parameters
            assert(this.dim == 1);
            ccd = CustomCircularDistribution(this.f);
        end
        
        function ctd = toCustomToroidal(this)
            % Convert to a custom toroidal distribution (only in 2D case)
            %
            % Returns:
            %   ctd (CustomToroidalDistribution)
            %       CustomToroidalDistribution with same parameters
            assert(this.dim == 2);
            ctd = CustomToroidalDistribution(this.f);
        end        
    end
    
    methods (Static)
        function chd = fromDistribution(dist)
            % Creates a CustomHypertoroidalDistribution from some other distribution
            %
            % Parameters:
            %   dist (AbstractHypertoroidalDistribution)
            %       distribution to convert
            % Returns:
            %   chd (CustomHypertoroidalDistribution)
            %       CustomHypertoroidalDistribution with identical pdf
            assert(isa(dist, 'AbstractHypertoroidalDistribution'));
            
            chd = CustomHypertoroidalDistribution(@(xa)dist.pdf(xa),dist.dim);
        end
    end
end
