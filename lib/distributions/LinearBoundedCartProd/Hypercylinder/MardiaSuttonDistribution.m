classdef MardiaSuttonDistribution < AbstractHypercylindricalDistribution
    % Gauss von Mises Distribution
    % 
    % Mardia, K. V. & Sutton, T. W. 
    % A Model for Cylindrical Variables with Applications
    % Journal of the Royal Statistical Society. Series B (Methodological), 
    % Wiley for the Royal Statistical Society, 1978, 40, pp. 229-233
    
    properties
        mu      %linear mu
        mu0     %circular mu
        kappa   %circular concentration
        rho1
        rho2
        sigma   %linear uncertainty
    end
    
    methods
        function this = MardiaSuttonDistribution(mu_, mu0_, kappa_, rho1_, rho2_, sigma_)
            % Constructor
            % check parameters
            assert(isscalar(mu_), 'mu must be scalar'); 
            assert(isscalar(mu0_), 'mu0 must be scalar'); 
            assert(isscalar(rho1_), 'rho1 must be scalar'); 
            assert(isscalar(rho2_), 'rho2 must be scalar'); 
            assert(isscalar(kappa_) && kappa_>0, 'kappa has to be a positive scalar');
            
            % assign parameters
            this.mu = mu_;
            this.mu0 = mod(mu0_, 2*pi);
            this.kappa = kappa_;
            this.rho1 = rho1_;
            this.rho2 = rho2_;
            this.sigma = sigma_;
            
            this.linD = 1;
            this.boundD = 1;
        end
        
        function p = pdf(this, xa)
            % Evaluate pdf at each column of xa
            %
            % Parameters:
            %   xa (d x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            
            assert(size(xa,1) ==  2);
            [muc, sigmac] = this.getMuSigma(xa(1,:));
            p = 1/(2*pi*besseli(0,this.kappa)) * exp(this.kappa * cos(xa(1,:) - this.mu0)) .* normpdf(xa(2,:), muc, sigmac);
        end
        
        function [muc, sigmac] = getMuSigma(this, xaPart)
            muc = this.mu + this.sigma * sqrt(this.kappa) * (this.rho1 * (cos(xaPart) -  cos(this.mu0)) + this.rho2 * (sin(xaPart) - sin(this.mu0)));
            rho = sqrt(this.rho1^2+this.rho2^2);
            sigmac = sqrt(this.sigma^2 * (1 - rho^2));
        end
              
        function m = mode(this)
            % Determines the mode of the distribution, i.e., the point
            % where the pdf is largest.
            %
            % Returns:
            %   m (linD + boundD x 1 column vector)
            %       the mode
            m = [this.mu0; this.mu];
        end
                                                            
        function C = hybridCovarianceNumerical(this)
            % Calculates covariance of [x1, x2, .., x_linD, cos(x_(linD+1), sin(x_(linD+1))]
            % numerically
            %
            % Returns:
            %   C (linD+2  x linD+2  matrix)
            %       covariance matrix of [x1, x2, .., x_linD, cos(x_(linD+1), sin(x_(linD+1))]
            s = this.sample(1000);
            S = [cos(s(1,:));sin(s(1,:));s(2:end,:)];
            C = cov(S');
        end
        
        function s = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (d x n matrix)
            %       n samples on R^(d-1) x [0,2pi)
            %
            assert(isscalar(n));
            vm = VMDistribution(this.mu0, this.kappa);            
            sVM = vm.sample(n);
            [muc, sigmac] = this.getMuSigma(sVM);
            sGauss = normrnd(muc, repmat(sigmac,1,n));
            s = [sVM; sGauss];
        end
        
        function C = linearCovariance(this)
            % Computes covariance of linear dimensions
            %
            % Returns:
            %   C (linD x linD)
            %       covariance matrix            
            C = this.sigma^2;
        end
        
        function vm = marginalizeLinear(this)
            vm = VMDistribution(this.mu0, this.kappa);
        end
    end   
end
