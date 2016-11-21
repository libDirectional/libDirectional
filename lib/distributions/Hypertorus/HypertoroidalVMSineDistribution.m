classdef HypertoroidalVMSineDistribution < AbstractHypertoroidalDistribution
    % Mardia, K. V.; Hughes, G.; Taylor, C. C. & Singh, H. 
    % A Multivariate von Mises Distribution with Applications to Bioinformatics
    % Canadian Journal of Statistics, Wiley Online Library, 2008, 36, 99-109
    
    properties
        mu
        kappa
        Lambda
        T
    end
    
    methods
        function this = HypertoroidalVMSineDistribution(mu_, kappa_, Lambda_)
            % Constructor
            assert(size(mu_,1) >= 1);
            assert(size(mu_,2) == 1);
            assert(all(size(kappa_) == size(mu_)));
            assert(all(kappa_ >=0));
            assert(size(Lambda_,1) == size(mu_,1));
            assert(size(Lambda_,2) == size(mu_,1));
            assert(all(all(Lambda_ == Lambda_')));
            assert(max(abs(diag(Lambda_))) == 0);
            
            this.dim = size(mu_,1);
            this.mu = mod (mu_, 2*pi);
            this.kappa = kappa_;
            this.Lambda = Lambda_;
            switch this.dim
                case 1
                    this.T = 2*pi*besseli(0,this.kappa); 
                case 2
                    s = @(m) nchoosek(2*m,m) * (this.Lambda(1,2)^2/4/this.kappa(1)/this.kappa(2))^m * besseli(m,this.kappa(1))*besseli(m,this.kappa(2));
                    this.T = 4*pi^2 * sum(arrayfun(s, 0:10));
                    % todo: choose number of summands in infinite series
                    % automatically
                otherwise
                    warning('unable to compute normalization constant');
                    this.T = 1; %todo compute normconst
            end
        end
        
        function p = pdf(this, xa)
            % Evaluate pdf at each column of xa.
            %
            % Parameters:
            %   xa (dim x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location            
            assert(size(xa,1) == this.dim);
           
            c = cos(xa - repmat(this.mu,1,size(xa,2)));
            s = sin(xa - repmat(this.mu,1,size(xa,2)));
            p = exp(this.kappa' * c + diag(s' * this.Lambda * s /2)') / this.T; % todo can we avoid computing the whole matrix if only diagonal is needed?
        end
        
        function hd = shift(this, shiftAngles)
            % Shift distribution by the given angles
            %
            % Parameters:
            %   shiftAngles (dim x 1 column vector) 
            %       angles to shift by
            % Returns:
            %   hd (HypertoroidalVMSineDistribution)
            %       shifted distribution            
            assert(all(size(shiftAngles) == [this.dim, 1]));
            
            hd = this;
            hd.mu = mod(this.mu+shiftAngles, 2*pi);
        end        
        
        function g = toGaussian(this)
            SigmaInv = diag(this.kappa) - this.Lambda;
            g = GaussianDistribution(this.mu, inv(SigmaInv));
        end
        
        function vm = toVM(this)
            % Convert to a VM distribution (only in 1D case)
            %
            % Returns:
            %   vm (VMDistribution)
            %       VMDistribution with same parameters
            assert(this.dim == 1);
            vm = VMDistribution(this.mu, this.kappa);
        end
        
        function tvm = toToroidalVMSine(this)
            % Convert to a toroidal VM sine distribution (only in 2D case)
            %
            % Returns:
            %   tvm (ToroidalVMSineDistribution)
            %       ToroidalVMSineDistribution with same parameters
            assert(this.dim == 2);
            tvm = ToroidalVMSineDistribution(this.mu, this.kappa, this.Lambda(1,2));
        end
        
        %todo conditionals
        
        %todo paramterestimation
    end
    
end

