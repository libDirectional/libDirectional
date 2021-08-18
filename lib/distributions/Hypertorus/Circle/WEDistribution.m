classdef WEDistribution < AbstractCircularDistribution
    % Wrapped exponential distribution
    %
    % see Sreenivasa Rao Jammalamadaka and Tomasz J. Kozubowski, "New 
    % Families of Wrapped Distributions for Modeling Skew Circular Data",
    % COMMUNICATIONS IN STATISTICS - Theory and Methods, Vol. 33, No. 9, 
    % pp. 2059-2074, 2004
    
    properties
        lambda
    end
    
    methods
        function this = WEDistribution(lambda_)
            % Constructor
            assert(isscalar(lambda_));
            assert(lambda_>0);
            this.lambda = lambda_;
        end
        
        function p = pdf(this, xa)
            % Evaluate pdf at each column of xa
            %
            % Parameters:
            %   xa (1 x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            assert(size(xa,1)==1);
            xa=mod(xa,2*pi);
            p = this.lambda*exp(-this.lambda*xa)/(1-exp(-2*pi*this.lambda));
        end
                
        function m = trigonometricMoment(this, n)
            % Calculate n-th trigonometric moment analytically
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (scalar)
            %       n-th trigonometric moment (complex number)
            m = 1/(1-1i*n/this.lambda);
        end
                
        function s = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (1 x n)
            %       n samples on the circle
            s = mod(exprnd(this.lambda,1,n),2*pi);
        end
        
        function result = entropy(this)
            % Calculates the entropy analytically 
            %
            % Returns:
            %   result (scalar)
            %       entropy of the distribution
            beta = exp(2*pi*this.lambda);
            result = 1 + log((beta-1)/this.lambda) - beta/(beta-1)*log(beta);
        end
    end
end

