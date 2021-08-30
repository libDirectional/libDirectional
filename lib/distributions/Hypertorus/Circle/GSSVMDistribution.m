classdef GSSVMDistribution < AbstractCircularDistribution
    properties
        mu (1,1) double {mustBeReal,mustBeFinite,mustBeNonNan}
        kappa (1,1) double {mustBeNonnegative}
        n (1,1) double {mustBeInteger}% skewnewssParameter
    end
    
    methods
        function this = GSSVMDistribution(mu, kappa, n)
            % Constructor
            this.mu = mod(mu,2*pi);
            this.kappa = kappa;
            this.n = n;
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
            switch this.n
                case 1
                    error('This n is not yet supported');
                case 2
                    error('This n is not yet supported');
                case 3
                    error('This n is not yet supported');
                case 4
                    error('This n is not yet supported');
                otherwise
                    error('This n is not yet supported');
            end
        end
        
        function dist = shift(this, angle)
            % Shift distribution by the given angle
            %
            % Parameters:
            %   shiftAngles (scalar)
            %       angle to shift by
            % Returns:
            %   hd (VMDistribution)
            %       shifted distribution
            assert(isscalar (angle));
            dist = this;
            dist.mu = this.mu+angle;
        end
        
    end
    
end

