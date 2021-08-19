classdef ToroidalVMRivestDistribution < AbstractToroidalDistribution
    % Bivariate von Mises, version by Rivest with two correlation
    % parameters
    %
    % Corresponds to A = [alpha, 0; 0, beta]
    %
    % Rivest, L.-P. 
    % A Distribution for Dependent Unit Vectors
    % Communications in Statistics - Theory and Methods, 1988, 17, 461-483
    
    properties
        mu      % 2 x 1 location parameter 
        kappa   % 2 x 1 concentration paramter
        alpha   % scalar correaltion parameter
        beta    % scalar correaltion parameter
        C       % normalization constant
    end
    
    methods
        function this = ToroidalVMRivestDistribution(mu_, kappa_, alpha_, beta_)
            % Constructor
            %
            % Parameters:
            %   mu_ (2 x 1)
            %       location parameter
            %   kappa_ (2 x 1)
            %       concentration parameter (>=0)
            %   lambda_ (scalar)
            %       correlation parameter
            assert(size(mu_,1)==2);
            assert(size(mu_,2)==1);
            assert(size(kappa_,1)==2);
            assert(size(kappa_,2)==1);
            assert(size(alpha_,1)==1);
            assert(size(alpha_,2)==1);
            assert(size(beta_,1)==1);
            assert(size(beta_,2)==1);
            assert(kappa_(1) >= 0);
            assert(kappa_(2) >= 0);
            
            this.mu = mod(mu_, 2*pi);
            this.kappa = kappa_;
            this.alpha = alpha_;
            this.beta = beta_;
            
            this.C = 1/this.normConst();
        end
        
        function Cinv = normConst(this)
            % Normalization constant
            % 
            % Returns:
            %   Cinv (scalar)
            %       inverse of the normalization constant 
            %
            s = @(j,l) besseli(j, this.kappa(1)) * besseli(l, this.kappa(2)) * besseli((j+l)/2, (this.alpha + this.beta)/2) * besseli((j-l)/2, (this.alpha - this.beta)/2);
            %todo: choose number of summands in infinite series
            %automatically
            n = 10;
            total = 0;
            for j=-n:n
                for l=-n:n
                    if mod(j+l ,2)==0 %todo could be optimized
                        total = total + s(j,l);
                    end
                end
            end
            Cinv = 4*pi^2 * total;
        end
        
        function p = pdf(this, xa)
            % Evaluate pdf at each column of xa
            %
            % Parameters:
            %   xa (2 x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            assert(size(xa, 1) == 2);
            
            p = this.C * exp ( ...
                  this.kappa(1) * cos(xa(1,:) - this.mu(1)) ...
                + this.kappa(2) * cos(xa(2,:) - this.mu(2)) ...
                + this.alpha * cos(xa(1,:) - this.mu(1)) .* cos(xa(2,:) - this.mu(2)) ... 
                + this.beta * sin(xa(1,:) - this.mu(1)) .* sin(xa(2,:) - this.mu(2)));
        end
        
        function m = trigonometricMoment(this,n)
            % Calculate n-th trigonometric moment, i.e., 
            % E([e^(inx_1); e^(inx_2)])
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (2x1)
            %       n-th trigonometric moment (complex vector)      
            if n == 1
                % derived by Gerhard Kurz
                m = 10;
                total1 = 0;
                total2 = 0;
                s1 = @(j,l) (besseli(j+1, this.kappa(1)) + besseli(j-1, this.kappa(1))) * besseli(l, this.kappa(2)) * besseli((j+l)/2, (this.alpha + this.beta)/2) * besseli((j-l)/2, (this.alpha - this.beta)/2);
                s2 = @(j,l) besseli(j, this.kappa(1)) * (besseli(l+1, this.kappa(2)) + besseli(l-1, this.kappa(2))) * besseli((j+l)/2, (this.alpha + this.beta)/2) * besseli((j-l)/2, (this.alpha - this.beta)/2);
                for j=-m:m
                    for l=-m:m
                        if mod(j+l ,2)==0 %todo could be optimized
                            total1 = total1 + s1(j,l);
                            total2 = total2 + s2(j,l);
                        end
                    end
                end
                
                m(1,1) = this.C * 2 * pi^2 * total1 *exp(1i * n * this.mu(1));
                m(2,1) = this.C * 2 * pi^2 * total2 *exp(1i * n * this.mu(2));
            else
                m = this.trigonometricMomentNumerical(n);
            end
        end
        
        function rhoc = circularCorrelationJammalamadaka(this)
            % Calculates Jammalamadaka's correlation coefficient numerically
            % see Jammalamadaka, S. R. & Sarma, Y. 
            % A Correlation Coefficient for Angular Variables 
            % Statistical Theory and Data Analysis II, 1988, 349-364
            %
            % Returns:
            %   rhoc (scalar)
            %       Jammalamadaka's correlation coefficient        
            
            m = 10;
            total0 = 0;            
            total1 = 0;
            total2 = 0;
            a = this.alpha;
            b = this.beta;
            k = this.kappa;
            s0 = @(j,l) besseli(j, k(1)) * besseli(l, k(2)) ...
                * ((besseli((j+l+2)/2, (a+b)/2) + besseli((j+l-2)/2, (a+b)/2))* besseli((j-l)/2, (a-b)/2) ...
                - besseli((j+l)/2, (a+b)/2) *(besseli((j-l+2)/2, (a-b)/2) + besseli((j-l-2)/2, (a-b)/2)));
            s1 = @(j,l) (besseli(j+2, k(1))/2 + besseli(j, k(1)) + besseli(j-2, k(1))/2 ) * besseli(l, k(2)) * besseli((j+l)/2, (a + b)/2) * besseli((j-l)/2, (a - b)/2);
            s2 = @(j,l) besseli(j, k(1)) * (besseli(l+2, k(2))/2 + besseli(l, k(2)) + besseli(l-2, k(2))/2) * besseli((j+l)/2, (a + b)/2) * besseli((j-l)/2, (a - b)/2);
            for j=-m:m
                for l=-m:m
                    if mod(j+l ,2)==0 %todo could be optimized
                        total0 = total0 + s0(j,l);
                        total1 = total1 + s1(j,l);
                        total2 = total2 + s2(j,l);
                    end
                end
            end
            
            EsinAsinB = this.C * pi^2 * total0;
            
            EsinAsquared = 1-this.C* 2*pi^2 * total1;
            EsinBsquared = 1-this.C* 2*pi^2 * total2;
            
            rhoc = EsinAsinB/sqrt(EsinAsquared * EsinBsquared);
        end
        
        function tvm = toToroidalVMMatrixDistribution(this)
            tvm = ToroidalVMMatrixDistribution(this.mu, this.kappa, diag([this.alpha, this.beta]));
        end
        
        function tvm = shift(this, shiftAngles)
            % Shift distribution by the given angles
            %
            % Parameters:
            %   shiftAngles (dim x 1 column vector) 
            %       angles to shift by
            % Returns:
            %   hd (ToroidalVMRivestDistribution)
            %       shifted distribution
            assert(all(size(shiftAngles) == [this.dim, 1]));
            
            tvm = this;
            tvm.mu = mod(this.mu+shiftAngles,2*pi);
        end        
    end
end

