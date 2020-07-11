classdef HypertoroidalUniformDistribution < AbstractHypertoroidalDistribution & AbstractUniformDistribution
    % Uniform distribution on the hypertorus.
    
    methods
        function this = HypertoroidalUniformDistribution(dim)
            this.dim = dim;
        end
        
        function m = trigonometricMoment(this, n)
            % Calculate n-th trigonometric moment analytically
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (dim x 1 column vector)
            %       n-th trigonometric moment (complex entries)
            m = double(n==0); % 1 for n=0, 0 otherwise
            m = repmat(m, this.dim, 1);
        end
        
        function result = entropy(this)
            % Calculates the entropy analytically 
            %
            % Returns:
            %   result (scalar)
            %       entropy of the distribution
            result = this.dim * log(2*pi);
        end
        
        function m = circularMean(~)
            % Calculate the circular mean
            %
            % Returns:
            %   m (scalar)
            %       circular mean in [0, 2pi)
            warning('MEAN:UNDEFINED','Circular uniform distribution does not have a unique mean');
            m = NaN;
        end
        
        function s = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       n samples on the hypertorus
            s = 2*pi*rand(this.dim,n);
        end
        
        function hu = shift(this, shiftAngles)
            % Shift distribution by the given angles
            %
            % Parameters:
            %   shiftAngles (dim x 1 column vector) 
            %       angles to shift by
            % Returns:
            %   hd (HypertoroidalUniformDistribution)
            %       shifted distribution
            assert(all(size(shiftAngles) == [this.dim, 1]));
            
            hu = this;
        end
    end
       
    methods (Sealed)
        function r = integral(this, l, r)
            % Calculates the integral of the pdf from l to r
            %
            % Parameters:
            %   l (dim x 1 column vector)
            %       left bound of integral in each dimension, default 0
            %   r (dim x 1 column vector)
            %       right bound of integral in each dimension, default 2*pi      
            if nargin < 2;  l = zeros(this.dim,1); end
            if nargin < 3;  r = 2*pi*ones(this.dim,1); end
            assert(all(size(l) == [this.dim, 1]));
            assert(all(size(r) == [this.dim, 1]));     
            
            volume = prod(r - l);
            r = 1/(2*pi)^this.dim * volume;
        end
    end 
end