classdef ToroidalVMCosineDistribution < AbstractToroidalDistribution
    % Bivariate von Mises, cosine version, (model with positive
    % interaction)
    %
    % Corresponds to A = [-kappa3, 0; 0, -kappa3]
    %
    % Mardia, K. V.; Taylor, C. C. & Subramaniam, G. K. 
    % Protein Bioinformatics and Mixtures of Bivariate von Mises Distributions for Angular Data Biometrics, 
    % Blackwell Publishing Inc, 2007, 63, 505-512
    %
    % Mardia, K. V. & Frellsen, 
    % J. Hamelryck, T.; Mardia, K. & Ferkinghoff-Borg, J. (Eds.) 
    % Statistics of Bivariate von Mises Distributions Bayesian Methods in Structural Bioinformatics, 
    % Springer Berlin Heidelberg, 2012, 159-178
    
    properties
        mu      % 2 x 1 location parameter 
        kappa   % 2 x 1 concentration paramter
        kappa3  % scalar correlation parameter, we do not put kappa3 into the kappa vector because it plays a fairly different role
        C       % normalization constant
    end
    
    methods
        function this = ToroidalVMCosineDistribution(mu_, kappa_, kappa3_)
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
            assert(size(kappa3_,1)==1);
            assert(size(kappa3_,2)==1);
            assert(kappa_(1) >= 0);
            assert(kappa_(2) >= 0);
            
            this.mu = mod(mu_, 2*pi);
            this.kappa = kappa_;
            this.kappa3 = kappa3_;
            
            this.C = 1/this.normConst();
        end
        
        function Cinv = normConst(this)
            % Normalization constant
            % 
            % Returns:
            %   Cinv (scalar)
            %       inverse of the normalization constant 
            %
            s = @(p) besseli(p, this.kappa(1))*besseli(p, this.kappa(2))*besseli(p, -this.kappa3); %todo why -kappa3?
            %todo: choose number of summands in infinite series
            %automatically
            Cinv = 4*pi^2 * (s(0) + 2*sum(arrayfun(s, 1:10)));
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
                - this.kappa3 * cos(xa(1,:) - this. mu(1) - xa(2,:) + this. mu(2)));
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
                % derived by Gerhard Kurz
                s1 = @(m) (besseli(m+1, this.kappa(1)) + besseli(m-1, this.kappa(1))) * besseli(m, this.kappa(2)) * besseli(m, -this.kappa3) ;
                s2 = @(m) besseli(m, this.kappa(1)) * (besseli(m+1, this.kappa(2)) + besseli(m-1, this.kappa(2))) * besseli(m, -this.kappa3) ;
                s = @(p) besseli(p, this.kappa(1))*besseli(p, this.kappa(2))*besseli(p, -this.kappa3);
                s1sum = s1(0)/2 + sum(arrayfun(s1, 1:10));
                s2sum = s2(0)/2 + sum(arrayfun(s2, 1:10));
                ssum = s(0)+ 2*sum(arrayfun(s, 1:10));
                m(1,1) = s1sum/ssum*exp(1i * n * this.mu(1));
                m(2,1) = s2sum/ssum*exp(1i * n * this.mu(2));
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
            
            %derived by Gerhard Kurz
            scc = @(m) (besseli(m+1, this.kappa(1)) + besseli(m-1, this.kappa(1))) * (besseli(m+1, this.kappa(2)) + besseli(m-1, this.kappa(2))) * besseli(m, -this.kappa3) ;
            EcosAcosB = this.C * 4*pi^2 *( besseli(1,this.kappa(1))* besseli(1, this.kappa(2))*besseli(0, -this.kappa3) + 1/2*sum(arrayfun(scc, 1:10)));
            sccss = @(m) besseli(m, this.kappa(1))*besseli(m, this.kappa(2)) * (besseli(m+1, -this.kappa3)+besseli(m-1, -this.kappa3));
            EcosAcosBsinAsinB = this.C * 4*pi^2 * (besseli(0, this.kappa(1)) * besseli(0, this.kappa(2)) * besseli(1, -this.kappa3) + sum(arrayfun(sccss, 1:10)));
            EsinAsinB = EcosAcosBsinAsinB - EcosAcosB;
            
            s1 = @(m) (besseli(m+2, this.kappa(1))/2 + besseli(m, this.kappa(1)) + besseli(m-2, this.kappa(1))/2) * besseli(m, this.kappa(2)) * besseli(m, -this.kappa3);
            s2 = @(m) (besseli(m+2, this.kappa(2))/2 + besseli(m, this.kappa(2)) + besseli(m-2, this.kappa(2))/2) * besseli(m, this.kappa(1)) * besseli(m, -this.kappa3);
            EsinAsquared = 1 - this.C * 4*pi^2* (1/2 * (besseli(0, this.kappa(1)) + besseli(2,this.kappa(1))) * besseli(0, this.kappa(2)) * besseli(0,this.kappa3) + sum(arrayfun(s1, 1:10)));
            EsinBsquared = 1 - this.C * 4*pi^2* (1/2 * (besseli(0, this.kappa(2)) + besseli(2,this.kappa(2))) * besseli(0, this.kappa(1)) * besseli(0,this.kappa3) + sum(arrayfun(s2, 1:10)));
            
            rhoc = EsinAsinB/sqrt(EsinAsquared * EsinBsquared);
        end
        
        function tvm = toToroidalVMMatrixDistribution(this)
            tvm = ToroidalVMMatrixDistribution(this.mu, this.kappa, -this.kappa3*eye(2,2));
        end
        
        function tvm = toToroidalVMRivestDistribution(this)
            tvm = ToroidalVMRivestDistribution(this.mu, this.kappa, -this.kappa3, -this.kappa3);
        end
                                   
        function twn = toToroidalWN(this)
            % Convert to ToroidalWN by method presented in Mardia 2012, eq.
            % (6.15)
            % 
            % Returns:
            %   twn (ToroidalWNDistribtution)
            %       TWN with similar pdf as proposed by Mardia 2012
                        
            if this.kappa3^2 >= (this.kappa(1) - this.kappa3) *(this.kappa(2) - this.kappa3) || this.kappa(1)-this.kappa3 <= 0 || this.kappa(2)-this.kappa3 <= 0
                warning('conversion may not be meaningful');
            end
            
            Cinv = [this.kappa(1) - this.kappa3, this.kappa3;
                this.kappa3, this.kappa(2)-this.kappa3];
            twn = ToroidalWNDistribution(this.mu,inv(Cinv));
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
            %see Mardia 2007, eq. (5)
            other = 3-dimension;
            f = @(x) 2*pi * this.C * besseli(0,sqrt(this.kappa(other)^2 + this.kappa3^2 - 2*this.kappa(other)*this.kappa3 *cos(x - this.mu(dimension)) )) .* exp(this.kappa(dimension) .* cos(x - this.mu(dimension)));
            dist = CustomCircularDistribution(f);
        end
        
        function tvm = shift(this, shiftAngles)
            % Shift distribution by the given angles
            %
            % Parameters:
            %   shiftAngles (dim x 1 column vector) 
            %       angles to shift by
            % Returns:
            %   hd (ToroidalVMCosineDistribution)
            %       shifted distribution
            assert(all(size(shiftAngles) == [this.dim, 1]));
            
            tvm = this;
            tvm.mu = mod(this.mu+shiftAngles,2*pi);
        end        
    end
    
end

