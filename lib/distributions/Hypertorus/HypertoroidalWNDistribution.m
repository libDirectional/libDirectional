classdef HypertoroidalWNDistribution < AbstractHypertoroidalDistribution
    % Hypertoroidal generalization of the wrapped normal distribution.
    %
    properties
        mu
        C
    end
    
    methods
        function this = HypertoroidalWNDistribution (mu_, C_)
            % Constructor
            %
            % Parameters:
            %   mu_ (dim x 1)
            %       mean vector of Gauss before wrapping
            %   C_ (dim x dim)
            %       covariance matrix of Gauss before wrapping
            
            % Check parameters
            assert(size(C_,1)==size(C_,2), 'C must be dim x dim');
            assert(isequal(C_,C_'), 'C must be symmetric');
            assert(all(eig(C_)>0), 'C must be positive definite');
            assert(isequal(size(mu_),[size(C_,2),1]), 'mu must be dim x 1'); 
            
            % Assign parameters
            this.mu=mod(mu_,2*pi);
            this.C=C_;
            this.dim=size(mu_,1);
        end
        
        function p = pdf(this,xa,m)
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
            assert(size(xa,1)==this.dim);
            if nargin==2,m=3;end
            
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
            gauss = GaussianDistribution(this.mu, this.C);
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
            assert(all(size(shiftAngles) == size(this.mu)));
            
            hd = this;
            hd.mu = mod(this.mu+shiftAngles,2*pi);
        end
    end
    
end
