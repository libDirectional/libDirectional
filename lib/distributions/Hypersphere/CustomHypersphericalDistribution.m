classdef CustomHypersphericalDistribution < AbstractHypersphericalDistribution & CustomDistribution
    % Hyperspherical distribution with custom pdf.
    
    methods
        function this = CustomHypersphericalDistribution(f_,dim_)
            % Constructor, it is the user's responsibility to ensure that f is a valid
            % hyperspherical density and takes arguments of the same form as
            % .pdf, i.e., it needs to be vectorized.
            % 
            % Parameters:
            %   f_ (function handle)
            %       pdf of the distribution
            %   dim_ (scalar)
            %       dimension of the real space in which the hypersphere is
            %       embedded
            this@CustomDistribution(f_,dim_);
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
            assert(isequal(size(p),[1,size(xa,2)])); % Validate output format of pdf is as expected
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
            %   chd (CustomHypersphericalDistribution)
            %       CustomHypersphericalDistribution with identical pdf
            assert(isa(dist, 'AbstractHypersphericalDistribution'));
            
            chd = CustomHypersphericalDistribution(@(xa)dist.pdf(xa),dist.dim);
        end
    end
end
