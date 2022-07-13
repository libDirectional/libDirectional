classdef HypercylindricalWNDistribution < AbstractHypercylindricalDistribution
    % Partially wrapped normal distribution 
    %
    % Reference:
    % Gerhard Kurz
    % Directional Estimation for Robotic Beating Heart Surgery
    % Karlsruhe Institute of Technology, 
    % Intelligent Sensor-Actuator-Systems Laboratory, 2015      
    
    properties
        mu (:,1) double
        C (:,:) double
    end
    
    methods
        function this = HypercylindricalWNDistribution(mu_, C_, boundD_)
            arguments
                mu_ (:,1) double
                C_ (:,:) double
                boundD_ (1,1) {mustBeInteger,mustBeNonnegative}
            end
            assert(all(size(C_) == [size(mu_,1),size(mu_,1)]), 'C must match size of mu');
            assert(all(all(C_ == C_')), 'C must be symmetric');
            assert(all(eig(C_) > 0), 'C must be positive definite');
            assert(boundD_<= size(mu_,1));
            
            this.boundD = boundD_;
            this.linD = size(mu_,1) - boundD_;
            this.dim = this.linD + this.boundD;
            
            % assign parameters
            this.mu = mu_;
            this.mu(1:boundD_) = mod(this.mu(1:boundD_),2*pi);
            this.C = C_;
        end
        
        function p = pdf(this, xa, m)
            % Evaluate pdf at each column of xa
            %
            % Parameters:
            %   xa (d x n)
            %       n locations where to evaluate the pdf
            %   m (scalar)
            %       number of summands in each direction
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            arguments
                this (1,1) HypercylindricalWNDistribution
                xa (:,:) double
                m (1,1) {mustBeInteger, mustBeNonnegative} = 3
            end
            assert(size(xa,1)==this.dim);
            % Do at most 1000 at once to prevent out of memory. We recurse
            % here into this very function!
            if size(xa,2)>1000
                p = NaN(1,size(xa,2));
                for i=1:size(xa,2)/1000
                    p(1,i*1000-999:i*1000)=this.pdf(xa(:,i*1000-999:i*1000));
                end
                p(1,i*1000+1:end)=this.pdf(xa(:,i*1000+1:end));
                return
            end
            
            offsets = arrayfun(@(i){-m:m},1:this.boundD);
            offsetMatrices = cell(length(offsets),1);
            [offsetMatrices{:}] = ndgrid(offsets{:});
            offsetsAsVectorsOnlyPeriodic = 2*pi*cell2mat(cellfun(@(mat){mat(:)'},offsetMatrices));
            offsetsAsVectorsWithLin = [offsetsAsVectorsOnlyPeriodic; zeros(this.linD, size(offsetsAsVectorsOnlyPeriodic,2))];
            evalPdfAt = repmat(offsetsAsVectorsWithLin,[1,size(xa,2)])+kron(xa,ones(1,size(offsetsAsVectorsWithLin,2)));
            allps = mvnpdf(evalPdfAt',this.mu',this.C);
            p = sum(reshape(allps,[],size(xa,2)),1);
        end
        
        function m = mode(this)
            % Determines the mode of the distribution, i.e., the point
            % where the pdf is largest.
            %
            % Returns:
            %   m (linD + boundD x 1 column vector)
            %       the mode
            arguments
                this (1,1) HypercylindricalWNDistribution
            end
            m = this.mu;
        end

        function dist = setMode(this, newMode)
            arguments
                this (1,1) HypercylindricalWNDistribution
                newMode (:,1) double
            end
            dist = this;
            dist.mu = newMode;
        end
        
        function mu = hybridMoment(this)
            % Calculates mean of [x1, x2, .., x_linD, cos(x_(linD+1), sin(x_(linD+1)), ..., cos(x_(linD+boundD), sin(x_(linD+boundD))]
            %
            % Returns:
            %   mu (linD+2 x 1)
            %       expectation value of [x1, x2, .., x_linD, cos(x_(linD+1), sin(x_(linD+1)), ..., cos(x_(linD+boundD), sin(x_(linD+boundD))]
            arguments
                this (1,1) HypercylindricalWNDistribution
            end
            mu = NaN(2*this.boundD+this.linD,1);
            mu(2*this.boundD+1:end, :) = this.mu(this.boundD+1:end,:);
            for i=1:this.boundD
                mu(2*i-1,:) = cos(this.mu(i))*exp(-this.C(i,i)/2);
                mu(2*i,:) =  sin(this.mu(i))*exp(-this.C(i,i)/2);
            end
        end
        
        function mean = hybridMean(this)
            arguments
                this (1,1) HypercylindricalWNDistribution
            end
            mean = this.mu;
        end
        
        function pwn = convolve(this, pwn2)
            % Calculate convolution of two PWN distributions (exact)
            %
            % Parameters:
            %   pwn2 (HypercylindricalWNDistribution)
            %       distribution to convolve with
            % Returns:
            %   pwn (HypercylindricalWNDistribution)
            %       convolution of this and pwn2
            arguments
                this (1,1) HypercylindricalWNDistribution
                pwn2 (1,1) HypercylindricalWNDistribution
            end
            
            assert(this.linD == pwn2.linD);
            assert(this.boundD == pwn2.boundD);
            
            % Copy this to preserve class when inheriting method
            pwn = this;
            pwn.mu = this.mu+pwn2.mu;
            pwn.mu(1:this.boundD) = mod(pwn.mu(1:this.boundD), 2*pi);
            pwn.C = this.C + pwn2.C;
        end
        
        function s = sample(this, n)
            arguments
                this (1,1) HypercylindricalWNDistribution
                n (1,1) {mustBeInteger, mustBePositive}
            end
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (3 x n matrix)
            %       n samples on [0,2pi) x R^2
            %
            % sample multivariate normal distribution, then wrap
            s = mvnrnd(this.mu, this.C, n)';
            s(1:this.boundD,:) = mod(s(1:this.boundD,:),2*pi);
        end
        
        function gauss = toGaussian(this)
            % Convert to a Gaussian with linD+boundD dimensions
            %
            % Returns:
            %   gauss (GaussianDistribution)
            %       GaussianDistribution with same parameters
            %
            % this is a simple conversion that just keeps the parameters
            % for large uncertanties, better conversions are possible
            arguments
                this (1,1) HypercylindricalWNDistribution
            end
            gauss = GaussianDistribution(this.mu, this.C);
        end
        
        function C = linearCovariance(this)
            % Computes covariance of linear dimensions
            %
            % Returns:
            %   C (linD x linD)
            %       covariance matrix
            arguments
                this (1,1) HypercylindricalWNDistribution
            end
            C = this.C(this.boundD+1:end,this.boundD+1:end);
        end        
        
        function gauss = marginalizeCircular(this)
            arguments
                this (1,1) HypercylindricalWNDistribution
            end
            gauss = GaussianDistribution(this.mu(this.boundD+1:end), this.C(this.boundD+1:end,this.boundD+1:end));
        end
        
        function wn = marginalizeLinear(this)
            arguments
                this (1,1) HypercylindricalWNDistribution
            end
            wn = HypertoroidalWNDistribution(this.mu(1:this.boundD), this.C(1:this.boundD,1:this.boundD));
        end        
    end   
end

