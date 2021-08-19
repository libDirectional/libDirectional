classdef HypertoroidalWNDistribution < AbstractHypertoroidalDistribution
    % Hypertoroidal generalization of the wrapped normal distribution.
    %
    properties
        mu (:,1) double
        C (:,:) double % Cannot enforce nonempty because we do not want to define a default value
    end
    
    methods
        function this = HypertoroidalWNDistribution(mu_, C_)
            % Constructor
            %
            % Parameters:
            %   mu_ (dim x 1)
            %       mean vector of Gauss before wrapping
            %   C_ (dim x dim)
            %       covariance matrix of Gauss before wrapping
            arguments
                mu_ (:,1) double
                C_ (:,:) double {mustBeNonempty}
            end
            
            % Check parameters
            assert(size(C_,1)==size(C_,2), 'C must be dim x dim');
            assert(issymmetric(C_), 'C must be symmetric');
            assert(all(eig(C_)>0), 'C must be positive definite');
            assert(size(mu_,1)==size(C_,2), 'mu must be dim x 1'); 
            
            % Assign parameters
            this.mu=mod(mu_,2*pi);
            this.C=C_;
            this.dim=size(mu_,1);
        end
        
        function p = pdf(this, xa, m)
            % Evaluate pdf at each column of xa.
            %
            % Parameters:
            %   xa (dim x n)
            %       n locations where to evaluate the pdf
            %   m
            %       -m*2pi to +m*2pi offsets will be considered
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            arguments
                this (1,1) HypertoroidalWNDistribution
                xa double
                m (1,1) {mustBeInteger,mustBePositive} = 3
            end
            assert(size(xa,1)==this.dim);
            
            % Do at most 1000 at once to prevent out of memory
            if size(xa,2)>1000
                p=NaN(1,size(xa,2));
                for i=1:size(xa,2)/1000
                    p(1,i*1000-999:i*1000)=this.pdf(xa(:,i*1000-999:i*1000));
                end
                p(1,i*1000+1:end)=this.pdf(xa(:,i*1000+1:end));
                return
            end
            offsets = arrayfun(@(i){-m:m},1:this.dim);
            offsetMatrices = cell(length(offsets),1);
            [offsetMatrices{:}] = ndgrid(offsets{:});
            offsetsAsVectors = 2*pi*cell2mat(cellfun(@(mat){mat(:)'},offsetMatrices));
            allVals = mvnpdf(kron(xa',ones(size(offsetsAsVectors,2),1))+repmat(offsetsAsVectors',size(xa,2),1),repmat(this.mu',size(offsetsAsVectors,2)*size(xa,2),1),this.C);
            p = sum(reshape(allVals,(1+2*m)^this.dim,[]));
        end
        
        function s = sample(this, n)
            % Obtain n samples from the distribution.
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       n samples on the torus
            arguments
                this (1,1) HypertoroidalWNDistribution
                n (1,1) double {mustBeInteger,mustBePositive}
            end
            s = mod(mvnrnd(this.mu, this.C, n),2*pi)'; % sample multivariate normal distribution, then wrap
        end
        
        function m = trigonometricMoment(this, n)
            % Calculate n-th trigonometric moment, i.e., 
            % E([e^(inx_1); e^(inx_2); ..., e^(inx_dim)])
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (dim x 1)
            %       n-th trigonometric moment (complex vector)
            arguments
                this (1,1) HypertoroidalWNDistribution
                n (1,1) double {mustBeInteger}
            end
            m = exp(arrayfun(@(i)1i*n*this.mu(i)-n^2*this.C(i,i)/2,1:this.dim).');
        end
        
        function gauss = toGaussian(this)
            % Convert to dim-dimensional Gaussian.
            %
            % Returns:
            %   gauss (GaussianDistribution)
            %       GaussianDistribution with same parameters
            %
            % This is a simple conversion that just keeps the parameters.
            % For large uncertanties, better conversions are possible.
            arguments
                this (1,1) HypertoroidalWNDistribution
            end
            gauss = GaussianDistribution(this.mu, this.C);
        end
        
        function wn = toWN(this)
            % Convert to a WN distribution (only in 1D case)
            %
            % Returns:
            %   wn (WNDistribution)
            %       WNDistribution with same parameters
            arguments
                this (1,1) HypertoroidalWNDistribution
            end
            assert(this.dim == 1);
            wn = WNDistribution(this.mu, sqrt(this.C));
        end
        
        function wn = toCircular(this)
            arguments
                this HypertoroidalWNDistribution
            end
            wn = WNDistribution(this.mu, this.C);
        end
        
        function twn = toToroidalWN(this)
            % Convert to a toroidal WN distribution (only in 2D case)
            %
            % Returns:
            %   twn (ToroidalWNDistribution)
            %       ToroidalWNDistribution with same parameters
            arguments
                this (1,1) HypertoroidalWNDistribution
            end
            assert(this.dim == 2);
            twn = ToroidalWNDistribution(this.mu, this.C);
        end
        
        function m = mode(this)
            % Determines the mode of the distribution, i.e., the point
            % where the pdf is largest.
            %
            % Returns:
            %   m (vector)
            %       the mode
            arguments
                this (1,1) HypertoroidalWNDistribution
            end
            m = this.mu;
        end
        
        function hd = shift(this, shiftAngles)
            % Shift distribution by the given angles
            %
            % Parameters:
            %   shiftAngles (dim x 1 column vector) 
            %       angles to shift by
            % Returns:
            %   hd (HypertoroidalWNDistribution)
            %       shifted distribution
            assert(all(size(shiftAngles) == [this.dim, 1]));
            
            hd = this;
            hd.mu = mod(this.mu+shiftAngles,2*pi);
        end
        
        function htvm = toHypertoroidalVM(this)
            % unpublished method
            
            % d is structured as 
            %    P
            %    Q
            %
            % W is structured as
            %    P P'    P Q'
            %    Q P'    Q Q'
            
            % preparation
            c = zeros(this.dim,1);
            for i=1:this.dim
                c(i) = exp(-1/2 * this.C(i,i)); % E(cos(theta_i))
            end
            
            % compute d
            d = zeros(this.dim*(this.dim+1)/2,1);
            % kappa part (P)
            for i=1:this.dim
                d(i) = c(i);
            end
            % Lambda part (Q)
            ind = this.dim+1;
            for i=1:this.dim
                for j=i+1:this.dim
                    d(ind) = c(i)*c(j) * sinh(this.C(i,j)); %E(sin(theta_i) sin(theta_j))
                    %why is there a 2 in eq. (2.1)?
                    ind = ind+1;
                end
            end
            
            % compute W
            W =  zeros(this.dim*(this.dim+1)/2);
            % kappa part (P*P')
            for i=1:this.dim
                W(i,i) = (1-c(i)^4)/2; % E(sin(theta_i)^2)
            end
            % Lambda part
            %todo P*Q'
            for l=1:this.dim
                ind = this.dim+1;
                for i=1:this.dim
                    for j=i+1:this.dim
                        W(l, ind) = 0; % E(- sin(theta(l) * something ) 
                    end
                end
            end
            
            %todo mirror
            
            %Q*Q'
            %diagonal entries
            ind = this.dim+1;
            for i=1:this.dim
                for j=i+1:this.dim
                    W(ind, ind) = (1-c(i)^4*c(j)^4 * cosh(4 * this.C(i,j)))/2; %E(cos^2(theta_i) sin^2(theta_j) + sin^2(theta_i) cos^2(theta_j))
                    ind = ind+1;
                end
            end
            %todo off diagonal parts
            
            
            param = W \ d;
            
            % construct htvm
            kappa = param(1:this.dim);
            Lambda = zeros(this.dim, this.dim); %todo fill into Lambda matrix
            htvm = HypertoroidalVMSineDistribution(this.mu, kappa, Lambda);
        end
    end
    
end
