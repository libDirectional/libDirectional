classdef GaussVMDistribution < AbstractHypercylindricalDistribution
    % Gauss von Mises Distribution
    % 
    % Horwood, J. T. & Poore, A. B. 
    % Gauss von Mises Distribution for Improved Uncertainty Realism in Space Situational Awareness 
    % SIAM/ASA Journal on Uncertainty Quantification, 2014, 2, 276-304
    
    properties
        mu
        P
        alpha
        beta
        Gamma
        kappa
        A
    end
    
    methods
        function this = GaussVMDistribution(mu_, P_, alpha_, beta_, Gamma_, kappa_)
            % Constructor
            % check parameters
            assert(size(mu_,2) == 1, 'mu must be a column vector'); 
            assert(size(P_,1) == size(mu_,1) && size(P_,2) == size(mu_,1), 'P and mu must have matching size');
            assert(all(all(P_ == P_')), 'P must be symmetric');
            assert(all(eig(P_) > 0), 'P must be positive definite');
            assert(isscalar(alpha_));
            assert(all(size(beta_) == size(mu_)), 'size of beta must match size of mu');
            assert(size(Gamma_,1) == size(mu_,1) && size(Gamma_,2) == size(mu_,1), 'Gamma and mu must have matching size');
            assert(all(all(Gamma_ == Gamma_')), 'Gamma must be symmetric');
            assert(isscalar(kappa_) && kappa_>0, 'kappa has to be a positive scalar');
            
            % assign parameters
            this.mu = mu_;
            this.P= P_;
            this.alpha = mod(alpha_, 2*pi); % not really necessary but might simplify testing
            this.beta = beta_;
            this.Gamma = Gamma_;
            this.kappa = kappa_;
            this.A = chol(P_);
            
            this.linD = size(mu_,1);
            this.boundD = 1;
            this.dim = this.linD + this.boundD;
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
            
            assert(size(xa,1) ==  this.linD + 1);
            p = mvnpdf(xa(2:end, :)', this.mu', this.P)' .* exp(this.kappa * cos( xa(1,:) - this.getTheta(xa(2:end,:)))) /(2*pi*besseli(0,this.kappa));
        end
        
        function Theta = getTheta(this, xa)
            % xa is only linear part
            z = this.A\(xa-repmat(this.mu,1,size(xa,2)));
            Theta = repmat(this.alpha,1,size(xa,2)) + this.beta' * z + 0.5 * sum((chol(this.Gamma)*z).^2, 1); %todo this only works for pos. semidefinite Gamma, general solution?
        end
               
        function m = mode(this)
            % Determines the mode of the distribution, i.e., the point
            % where the pdf is largest.
            %
            % Returns:
            %   m (linD + boundD x 1 column vector)
            %       the mode
            %
            % See Horwood, 4.1 
            m = [this.alpha; this.mu];
        end
                                             
        function mu = hybridMoment(this)
            % Calculates mean of [x1, x2, .., x_linD, cos(x_(linD+1), sin(x_(linD+1))]
            %
            % Returns:
            %   mu (linD+2 x 1)
            %       expectation value of [x1, x2, .., x_linD, cos(x_(linD+1), sin(x_(linD+1))]
            
            % see Horwood, 4.3
            M = eye(this.linD) - 1i * this.Gamma;
            eiphi = 1/sqrt(det(M)) * besselratio(0,this.kappa) * exp(1i * this.alpha  - 0.5 * this.beta' / M * this.beta );
            mu = [real(eiphi); imag(eiphi); this.mu];
            
            
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
            sGauss = mvnrnd(this.mu, this.P, n)';
            vm = VMDistribution(0, this.kappa);
            sVM = mod(vm.sample(n) + this.getTheta(sGauss), 2*pi);
            s = [sVM;sGauss];
        end
        
        function [d,w] = sampleDeterministicHorwood(this)
            %todo 
            % Horwood 5.1
            dim = this.linD;
            
            % solution for canonical form (mu = 0, P=I,
            % alpha=beta=Gamma=0)
            B = @(p, kappa) 1-besseli(p,kappa)/besseli(0,kappa);
            xi = sqrt(3);
            eta = acos(B(2,this.kappa)/2/B(1,this.kappa) - 1);
            wxi0 = 1/6;
            weta0 = B(1,this.kappa)^2/(4*B(1,this.kappa) - B(2,this.kappa));
            w00 = 1 - 2*weta0 - 2*dim*wxi0;
            N00 = zeros(dim+1,1);
            Neta0 = [zeros(dim, 2); -eta eta];
            Nxi0 = zeros(dim+1, 2*dim);
            for i=1:dim
                Nxi0(i,2*i-1) = -xi;
                Nxi0(i,2*i) = xi;
            end
            d = [N00 Neta0 Nxi0];
            w = [w00 weta0 weta0 repmat(wxi0,1,2*dim)];
            
            % transform back
            d(1,:) = mod(d(1,:) + this.getTheta(d(2:end,:)),2*pi);
            d(2:end,:) = this.A * d(2:end,:) + repmat(this.mu,1,2*dim+3);
        end

        function gauss = toGaussian(this)
            % Convert to Gaussian
            %
            % Returns:
            %   gauss (GaussianDistribution)
            %       GaussianDistribution with same parameters         
            
            % Uses approximation only valid for large kappa and small
            % Gamma, see Horwood, 4.6
            
            mtmp = [this.alpha; this.mu];
%             Atmp = [this.A, zeros(this.linD,1) ; this.beta', 1/sqrt(this.kappa)];
            Atmp = [1/sqrt(this.kappa),this.beta';zeros(this.linD,1),this.A];
            Ptmp = Atmp*Atmp';
            gauss = GaussianDistribution(mtmp, Ptmp);
        end
        
        function C = linearCovariance(this)
            % Computes covariance of linear dimensions
            %
            % Returns:
            %   C (linD x linD)
            %       covariance matrix            
            C = this.P;
        end
               
        function gauss = marginalizeCircular(this)
            gauss = GaussianDistribution(this.mu, this.P);
        end
    end   
end
