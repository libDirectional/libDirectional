classdef CircularUniformDistribution < AbstractCircularDistribution & HypertoroidalUniformDistribution
    % Circular uniform distribution
    %
    % see Sreenivasa Rao Jammalamadaka and A. SenGupta, "Topics in Circular 
    % Statistics", 2001, Sec. 2.2.1, page 33.
    
    methods
        function this = CircularUniformDistribution()
            % Constructor
            this@HypertoroidalUniformDistribution(1);
        end
        
        function cu = shift(~,~)
            cu = CircularUniformDistribution;
        end
    end
    
    methods (Static) % Make them accessible via the class as well because they are independent of parameters        
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
    end
end

