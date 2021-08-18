classdef GSSCDistribution < AbstractCircularDistribution
    properties
        mu (1,1) double {mustBeReal,mustBeFinite,mustBeNonNan}
        rho (1,1) double {mustBeNonnegative,mustBeGreaterThan(rho,-0.5),mustBeLessThan(rho,0.5)}
        n (1,1) double {mustBeInteger}% skewnewssParameter
        lambda (1,1) double = 1
    end
    
    methods
        function this = GSSCDistribution(mu_, rho_, n_, lambda_)
            arguments
                mu_ (1,1) double {mustBeReal,mustBeFinite,mustBeNonNan}
                rho_ (1,1) double {mustBeNonnegative,mustBeGreaterThan(rho_,-0.5),mustBeLessThan(rho_,0.5)}
                n_ (1,1) double {mustBeInteger}% skewnewssParameter
                lambda_ (1,1) double = 1
            end
            % Constructor
            this.mu = mod(mu_,2*pi);
            this.rho = rho_;
            this.n = n_;
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
            cosTerm = 1+2*this.rho*cos(xa-this.mu);
            sinTerm = 1+this.lambda*sin(xa-this.mu);
            switch this.n
                case 1
                    p = 1/(2*pi)* cosTerm*sinTerm;
                case 2
                    p = 1/(2*pi)* 2/(2+this.lambda^2)*cosTerm.*sinTerm.^2;
                case 3
                    p = 1/(2*pi)* 2/(2+3*this.lambda^2)*cosTerm.*sinTerm.^3;
                case 4
                    p = 1/(2*pi)* 8/(8+24*this.lambda^2+3*this.lambda^4)*cosTerm.*sinTerm.^4;
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
            dist.mu=dist.mu+angle;
        end
        
    end
    
end

