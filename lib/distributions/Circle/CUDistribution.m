classdef CUDistribution < AbstractCircularDistribution
    % Circular uniform distribution
    %
    % see Sreenivasa Rao Jammalamadaka and A. SenGupta, "Topics in Circular 
    % Statistics", 2001, Sec. 2.2.1, page 33.
    
    methods
        function this = CUDistribution()
            % Constructor. No real initilization required. 
            % Strictly speaking, objects of this class serve no other
            % purpose than to allow users to always work with distribution
            % objects (e.g. for conversion to Fourier representation)
        end
    end
    
    methods (Static) % Make them accessible via the class as well because they are independent from parameters
        function p = pdf(xa)
            % Evaluate pdf at each column of xa
            %
            % Parameters:
            %   xa (1 x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            assert(size(xa,1)==1);
            p=ones(1,size(xa,2))/(2*pi);
        end
        
        function val = cdf(xa,startingPoint)
            % Evaluate cumulative distribution function 
            %
            % Parameters:
            %   xa (1 x n)
            %       points where the cdf should be evaluated
            %   startingPoint (scalar)
            %       point where the cdf is zero (starting point can be
            %       [0,2pi) on the circle, default 0
            % Returns:
            %   val (1 x n)
            %       cdf evaluated at columns of xa
            if nargin <=2
                startingPoint = 0;
            end
            val = (xa - startingPoint)/(2*pi);
            val (val<0) = val (val<0) + 1;
        end
                
        function m = trigonometricMoment(n)
            % Calculate n-th trigonometric moment analytically
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (scalar)
            %       n-th trigonometric moment (complex number)
            m=double(n==0); % 1 for n=0, 0 otherwise
        end
        
        function result = integral(l,r)
            % Calculates the integral of the pdf from l to r analytically
            %
            % Parameters:
            %   l (scalar)
            %       left bound of integral, default 0
            %   r (scalar)
            %       right bound of integral, default 2*pi
            % Returns:
            %   result (scalar)
            %       value of the integral
            if nargin < 1
                l = 0;
            end
            if nargin < 2
                r = 2*pi;
            end
            result = (r-l)/2/pi;
        end
        
        function m = circularMean()
            % Calculate the circular mean
            %
            % Returns:
            %   m (scalar)
            %       circular mean in [0, 2pi)
            warning('MEAN:UNDEFINED','Circular uniform distribution does not have a unique mean');
            m=NaN;
        end
        
        function result = entropy()
            % Calculates the entropy analytically 
            %
            % Returns:
            %   result (scalar)
            %       entropy of the distribution
            result = log(2*pi);
        end
        
        function s = sample(n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (1 x n)
            %       n samples on the circle
            s=2*pi*rand(1,n);
        end
    end
end

