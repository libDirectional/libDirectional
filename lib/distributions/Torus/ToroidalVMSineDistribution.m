classdef ToroidalVMSineDistribution < AbstractToroidalDistribution
    % Bivariate von Mises, sine version
    %
    % see http://en.wikipedia.org/wiki/Bivariate_von_Mises_distribution
    % and Singh H, Hnizdo V, Demchuk E (2002) 
    % Probabilistic model for two dependent circular variables
    % Biometrika 89: 719–723
    
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
        
    end
    
end

