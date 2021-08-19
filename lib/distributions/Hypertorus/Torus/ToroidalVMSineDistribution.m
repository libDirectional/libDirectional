classdef ToroidalVMSineDistribution < AbstractToroidalDistribution
    % Bivariate von Mises, sine version
    %
    % see http://en.wikipedia.org/wiki/Bivariate_von_Mises_distribution
    % and Singh H, Hnizdo V, Demchuk E (2002) 
    % Probabilistic model for two dependent circular variables
    % Biometrika 89: 719�723
    
    properties
        mu      % 2 x 1 location parameter 
        kappa   % 2 x 1 concentration paramter
        lambda  % scalar correlation parameter
        C       % normalization constant
    end
    
    methods
        function this = ToroidalVMSineDistribution(mu_, kappa_, lambda_)
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
            assert(size(lambda_,1)==1);
            assert(size(lambda_,2)==1);
            assert(kappa_(1) >= 0);
            assert(kappa_(2) >= 0);
            
            this.mu = mod(mu_, 2*pi);
            this.kappa = kappa_;
            this.lambda = lambda_;
            
            this.C = 1/this.normConst();
        end
        
        function Cinv = normConst(this)
            % Normalization constant
            % 
            % Returns:
            %   Cinv (scalar)
            %       inverse of the normalization constant 
            %
            % Formula from Singh 2002, Theorem 1, eq. (2.3)
            s = @(m) nchoosek(2*m,m) * (this.lambda^2/4/this.kappa(1)/this.kappa(2))^m * besseli(m,this.kappa(1))*besseli(m,this.kappa(2));
            %todo: choose number of summands in infinite series
            %automatically
            Cinv = 4*pi^2 * sum(arrayfun(s, 0:10));
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
                + this.lambda * sin(xa(1,:) - this. mu(1)).*sin(xa(2,:) - this. mu(2)));
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
            if n==1
                % Sing 2002, Theorem 2b
                s1 = @(m) nchoosek(2*m,m) * (this.lambda^2/4/this.kappa(1)/this.kappa(2))^m * besseli(m+1,this.kappa(1))*besseli(m,this.kappa(2));
                s2 = @(m) nchoosek(2*m,m) * (this.lambda^2/4/this.kappa(1)/this.kappa(2))^m * besseli(m,this.kappa(1))*besseli(m+1,this.kappa(2));
                s = @(m) nchoosek(2*m,m) * (this.lambda^2/4/this.kappa(1)/this.kappa(2))^m * besseli(m,this.kappa(1))*besseli(m,this.kappa(2));
                s1sum = sum(arrayfun(s1, 0:10));
                s2sum = sum(arrayfun(s2, 0:10));
                ssum = sum(arrayfun(s, 0:10));
                m(1,1) = s1sum/ssum*exp(1i * n * this.mu(1));
                m(2,1) = s2sum/ssum*exp(1i * n * this.mu(2));
            else
                m = this.trigonometricMomentNumerical(n);
            end
        end
        
        function rho = circularCorrelationJammalamadaka(this)
            % Circular correlation as defined by Jammalamadaka, closed-form solution derived by
            % Gerhard Kurz
            sinAsinB = @(m) nchoosek(2*m,m) * m * (this.lambda^2/4/this.kappa(1)/this.kappa(2))^m * besseli(m,this.kappa(1))*besseli(m,this.kappa(2));
            %cosAsquared = @(m) nchoosek(2*m,m) * (this.lambda^2/4/this.kappa(1)/this.kappa(2))^m * besseli(m+2,this.kappa(1))*besseli(m,this.kappa(2));
            cosAsquared = @(m) nchoosek(2*m,m) * (this.lambda^2/4/this.kappa(2))^m  * (besseli(m+2,this.kappa(1))/this.kappa(1)^m + besseli(m+1,this.kappa(1))/this.kappa(1)^(m+1)) *besseli(m,this.kappa(2));
            %cosBsquared = @(m) nchoosek(2*m,m) * (this.lambda^2/4/this.kappa(1)/this.kappa(2))^m * besseli(m,this.kappa(1))*besseli(m+2,this.kappa(2));
            cosBsquared = @(m) nchoosek(2*m,m) * (this.lambda^2/4/this.kappa(1))^m  * (besseli(m+2,this.kappa(2))/this.kappa(2)^m + besseli(m+1,this.kappa(2))/this.kappa(2)^(m+1)) *besseli(m,this.kappa(1));
            EsinAsinB = 8*pi^2*this.C/this.lambda*sum(arrayfun(sinAsinB, 0:10));
            EsinAsquared = 1 - this.C * 4 * pi^2 * sum(arrayfun(cosAsquared,0:10));
            EsinBsquared = 1 - this.C * 4 * pi^2 * sum(arrayfun(cosBsquared,0:10));
            rho = EsinAsinB/sqrt(EsinAsquared * EsinBsquared);
        end

        function tvm = toToroidalVMMatrixDistribution(this)
            tvm = ToroidalVMMatrixDistribution(this.mu, this.kappa, [0, 0; 0, this.lambda]);
        end
        
        function tvm = toToroidalVMRivestDistribution(this)
            tvm = ToroidalVMRivestDistribution(this.mu, this.kappa, 0, this.lambda);
        end
                                   
        function twn = toToroidalWN(this)
            % Convert to ToroidalWNDistribution by method presented in Singh 2002, eq (2.1)
            % 
            % Returns:
            %   twn (ToroidalWNDistribution)
            %       TWN with similar pdf as proposed by Singh 2002
            if this.lambda^2>= this.kappa(1)*this.kappa(2)
                warning('conversion may not be meaningful');
            end
            si1squared = this.kappa(2)/(this.kappa(1)*this.kappa(2) - this.lambda^2);
            si2squared = this.kappa(1)/(this.kappa(1)*this.kappa(2) - this.lambda^2);
            rho = this.lambda/sqrt(this.kappa(1)*this.kappa(2));
            si12 = sqrt(si1squared)*sqrt(si2squared)*rho;
            C_ = [si1squared, si12; si12, si2squared];
            twn = ToroidalWNDistribution(this.mu,C_);
        end
        
        function dist = marginalizeTo1D(this, dimension)
            % Get marginal distribution in first or second dimension, i.e., 
            % f(x_1) or f(x_2), respectively.
            %
            % Parameters:
            %   dimension (scalar)
            %       the marginal in which dimension to calculate (1 or 2),
            %       the other dimension is marginalized out
            % Returns:
            %   dist (CustomCircularDistribution)
            %       marginal distribution (marginals are NOT VM distributed in general)
            assert(dimension==1 || dimension==2);
            %see Singh 2002, eq (2.2)
            other = 3-dimension;
            f = @(x) 2*pi * this.C * besseli(0,sqrt(this.kappa(other)^2 + this.lambda^2 *sin(x - this.mu(dimension)).^2 )) .* exp(this.kappa(dimension) .* cos(x - this.mu(dimension)));
            dist = CustomCircularDistribution(f);
        end
        
        function tvm = shift(this, shiftAngles)
            % Shift distribution by the given angles
            %
            % Parameters:
            %   shiftAngles (dim x 1 column vector) 
            %       angles to shift by
            % Returns:
            %   hd (ToroidalVMSineDistribution)
            %       shifted distribution
            assert(all(size(shiftAngles) == [this.dim, 1]));
            
            tvm = this;
            tvm.mu = mod(this.mu+shiftAngles,2*pi);
        end        
    end
    
end

