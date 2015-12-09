classdef CustomHypertoroidalDistribution < AbstractHypertoroidalDistribution
    % Hypertoroidal distribution with custom pdf.
    
    properties
        f
    end
    
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
            assert(isa(f_, 'function_handle'));
            assert(isscalar(dim_) && dim_>=1);
            this.dim = dim_;
            this.f = f_;
        end
        
        function p = pdf(this, xa)
            % Evaluate pdf at each column of xa.
            %
            % Parameters:
            %   xa (dim x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            assert(size(xa,1)==this.dim);
            p = this.f(xa);
        end
        
        function hd = shift(this, shiftAngles)
            % Shift distribution by the given angles
            %
            % Parameters:
            %   shiftAngles (dim x 1 column vector) 
            %       angles to shift by
            % Returns:
            %   hd (CustomHypertoroidalDistribution)
            %       shifted distribution            
            assert(all(size(shiftAngles) == [this.dim, 1]));
            
            hd = this;
            hd.f = @(xa) this.f(xa-repmat(shiftAngles,[1,size(xa,2)]));
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
