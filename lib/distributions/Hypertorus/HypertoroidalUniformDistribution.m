classdef HypertoroidalUniformDistribution < AbstractHypertoroidalDistribution
    % Uniform distribution on the hypertorus.
    
    methods
        function this = HypertoroidalUniformDistribution(dim)
            this.dim = dim;
        end
        
        function vals = pdf(this,xa)
            vals = 1/(2*pi)^this.dim*ones(1,size(xa,2));
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
            result = log(2*pi);
        end
        
        function m = circularMean(this)
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
    end 
end