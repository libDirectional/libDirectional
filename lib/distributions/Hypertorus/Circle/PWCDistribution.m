classdef PWCDistribution < AbstractCircularDistribution
    % Piecewise constant (i.e. discrete) distribution, similar to a
    % histogram.
    %
    % Gerhard Kurz, Florian Pfaff, Uwe D. Hanebeck,
    % Discrete Recursive Bayesian Filtering on Intervals and the Unit Circle
    % Proceedings of the 2016 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI 2016), 
    % Baden-Baden, Germany, September 2016.
    
    properties
        w % weights as row vector
    end
    
    methods
        function this = PWCDistribution(w_)
            % Constructor
            %
            % Parameters:
            %   w (1 x n)
            %       weight for each interval            
            assert(size(w_,1)==1);
            this.w = w_/(mean(w_)*2*pi); % normalize automatically
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
            xa = mod(xa,2*pi);
            idx = floor(xa/2/pi*length(this.w));
            p = this.w(1+idx);
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
            interv = zeros(1,length(this.w));
            for j=1:length(this.w)
                l = PWCDistribution.leftBorder(j, length(this.w));
                r = PWCDistribution.rightBorder(j, length(this.w));
                c = PWCDistribution.intervalCenter(j,length(this.w));
                interv(j) = this.pdf(c) * (exp(1i*n*r) - exp(1i*n*l)); %integral from l to r
            end
            m = -1i/n*sum(interv);
        end
        
        function e = entropy(this)
            % Calculates the entropy analytically
            %
            % Returns:
            %   e (scalar)
            %       entropy of the distribution            
            e = -2*pi/length(this.w) * sum(this.w .* log(this.w));
        end
    end
        
    methods (Static)
        function w = calculateParametersNumerically(pdf, n)
            % Calculates the weights from a pdf
            %
            % Parameters:
            %   pdf (function handle)
            %       pdf of a circular density
            %   n (scalar)
            %       number of discretization steps
            % Returns:
            %   w (1 x n)
            %       weights of the corresponding PWC
            assert(isa(pdf, 'function_handle'));
            assert(n>=1);
            w = zeros(1,n);
            for j=1:n
                l = PWCDistribution.leftBorder(j,n);
                r = PWCDistribution.rightBorder(j,n);
                w(j) = integral(pdf, l, r);
            end
        end
        
        function l = leftBorder(m, n)
            % Calculates the left border of the m-th interval for a total of n intervals
            %
            % Parameters:
            %   n (scalar)
            %       total number of intervals
            %   m (scalar)
            %       number of intervals whose border is to be computed
            % Returns:
            %   m (scalar)
            %       left border of the m-th interval
            assert(1 <= m && m<=n);
            l = 2*pi/n*(m-1);
        end
        
        function r = rightBorder(m, n)
            % Calculates the right border of the m-th interval for a total of n intervals
            %
            % Parameters:
            %   n (scalar)
            %       total number of intervals
            %   m (scalar)
            %       number of intervals whose border is to be computed
            % Returns:
            %   m (scalar)
            %       right border of the m-th interval            
            assert(1 <= m && m<=n);
            r = 2*pi/n*m;
        end
        
        function c = intervalCenter(m, n)
            % Calculates the center of the m-th interval for a total of n intervals
            %
            % Parameters:
            %   n (scalar)
            %       total number of intervals
            %   m (scalar)
            %       number of intervals whose center is to be computed
            % Returns:
            %   m (scalar)
            %       center of the m-th interval            
            assert(1 <= m && m<=n);
            c = 2*pi/n*(m-0.5);
        end
    end    
end

