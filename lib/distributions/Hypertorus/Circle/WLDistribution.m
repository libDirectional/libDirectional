classdef WLDistribution < AbstractCircularDistribution
    % Wrapped Laplace distribution
    %
    % see Sreenivasa Rao Jammalamadaka and Tomasz J. Kozubowski, 
    % "New Families of Wrapped Distributions for Modeling Skew Circular Data",
    % COMMUNICATIONS IN STATISTICS - Theory and Methods, Vol. 33, No. 9, 
    % pp. 2059-2074, 2004
    
    properties
        lambda
        kappa
    end
    
    methods
        function this = WLDistribution(lambda_, kappa_)
            % Constructor
            assert(isscalar(lambda_));
            assert(isscalar(kappa_));
            assert(lambda_>0);
            assert(kappa_>0);
            this.lambda = lambda_;
            this.kappa = kappa_;
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
            p = this.lambda*this.kappa/(1+this.kappa^2) * ( exp(-this.lambda*this.kappa*xa)/(1-exp(-2*pi*this.lambda*this.kappa)) + exp(this.lambda/this.kappa*xa)/(exp(2*pi*this.lambda/this.kappa) - 1));  
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
            m = 1/(1-1i*n/this.lambda/this.kappa) / (1 + 1i * n /(this.lambda/this.kappa));
        end
    end
end

