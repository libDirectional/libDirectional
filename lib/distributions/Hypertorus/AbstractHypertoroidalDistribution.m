classdef (Abstract) AbstractHypertoroidalDistribution < AbstractPeriodicDistribution
    % Base class for distributions on the hypertorus [0,2pi)^dim (Cartesian
    % product of dim circles)
    % Convention for Dimension: dim=1 is a circle, dim=2 a torus.
    
    methods
        function p = plot(this, varargin)
            % Create an appropriate plot of the pdf
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to plot/surf command
            % Returns:
            %   p (scalar)
            %       plot handle
            resolutionIndex = find(cellfun(@(c)isequal(c,'resolution'),varargin));
            if ~isempty(resolutionIndex)
                resolution = varargin{resolutionIndex+1};
                varargin([resolutionIndex,resolutionIndex+1]) = [];
            else
                resolution = 128;
            end
            switch this.dim
                case 1
                    theta = linspace(0,2*pi,resolution);
                    ftheta = this.pdf(theta);
                    p = plot(theta, ftheta,varargin{:});
                case 2
                    step = 2*pi/resolution;
                    [alpha,beta] = meshgrid(0:step:2*pi,0:step:2*pi);
                    f = this.pdf([alpha(:)'; beta(:)']);
                    f = reshape(f,size(alpha,1), size(alpha,2));
                    p = surf(alpha, beta, f, varargin{:});
                case 3
                    stepCirc = 0.5;
                    [X,Y,Z] = sphere(4);
                    clf 
                    hold on
                    color = jet;
                    [gridx,gridy,gridz]=ndgrid(0:stepCirc:2*pi);
                    fgrid=reshape(this.pdf([gridx(:)';gridy(:)';gridz(:)']),size(gridx));
                    fmax=max(fgrid(:));
                    sizes=0.5*stepCirc*fgrid/fmax;
                    arrayfun(@(x,y,z,currSize) surf(currSize*X+x,currSize*Y+y,currSize*Z+z, 'facecolor', color(1+floor(currSize*126/stepCirc),:)),...
                        gridx(sizes>0.01),gridy(sizes>0.01),gridz(sizes>0.01),sizes(sizes>0.01)); % Only use grid points for which size>0.01
                    hold off
                    xlabel('x_1')
                    setupAxisCircular('x','y','z')
                    ylabel('x_2')
                    zlabel('x_3')
                    view(40,20)
                    grid
                otherwise
                    error('Plotting for this dimension is currently not supported');
            end
        end
        
        function m = circularMean(this)
            % Calculate the (componentwise) circular mean
            %
            % Warning: atan2(0,0) returns 0, not NaN!
            %
            % Returns:
            %   m (dim x 1)
            %       circular mean in [0, 2pi)^2
            a = this.trigonometricMoment(1);
            m = mod(atan2(imag(a),real(a)), 2*pi);
        end
        
        function m = meanDirection(this)
            % More general name for circular mean.
            m = circularMean(this);
        end
        
        function m = mode(this)
            m = modeNumerical(this);
        end
        
        function m = modeNumerical(this)
            m = fminunc(@(x)-this.pdf(x),zeros(this.dim,1));
        end
        
        function m = trigonometricMoment(this, n)
            % Calculate n-th trigonometric moment, i.e., 
            % E([e^(inx_1); e^(inx_2);...^e(inx_d)])
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (d x 1)
            %       n-th trigonometric moment (complex vector)             
            m = trigonometricMomentNumerical(this, n);
        end
        
        function result = integral(this, l, r)
            % Calculates the integral of the pdf from l to r
            %
            % Parameters:
            %   l (dim x 1 column vector)
            %       left bound of integral in each dimension, default 0
            %   r (dim x 1 column vector)
            %       right bound of integral in each dimension, default 2*pi
            arguments
                this (1,1) AbstractHypertoroidalDistribution
                l (:,1) double = zeros(this.dim,1);
                r (:,1) double = 2*pi*ones(this.dim,1);
            end
            assert(size(l,1)==this.dim);
            assert(size(r,1)==this.dim);
            
            result = this.integralNumerical(l, r);
        end
        
        function result = integralNumerical(this, l, r)
            % Numerically calculates the integral of the pdf from l to r
            %
            % Parameters:
            %   l (dim x 1 column vector)
            %       left bound of integral in each dimension, default 0
            %   r (dim x 1 column vector)
            %       right bound of integral in each dimension, default 2*pi 
            arguments
                this (1,1) AbstractHypertoroidalDistribution
                l (:,1) double = zeros(this.dim,1);
                r (:,1) double = 2*pi*ones(this.dim,1);
            end
            assert(size(l,1)==this.dim);
            assert(size(r,1)==this.dim);     
            
            switch this.dim
                case 1
                    result = integral(@(x) this.pdf(x), l, r);
                case 2
                    f = @(x,y) reshape(this.pdf([x(:)';y(:)']), size(x));
                    result = integral2(f, l(1), r(1), l(2), r(2));
                case 3
                    f = @(x,y,z) reshape(this.pdf([x(:)';y(:)';z(:)']), size(x));
                    result = integral3(f, l(1), r(1), l(2), r(2), l(3), r(3));
                otherwise
                    error('Numerical moment calculation for this dimension is currently not supported');
            end
        end
        
        function m = trigonometricMomentNumerical(this, n)
            % Calculate n-th trigonometric moment numerically, i.e., 
            % E([e^(inx_1); e^(inx_2);...^e(inx_d)])
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (dx1)
            %       n-th trigonometric moment (complex vector)  
            switch this.dim
                case 1
                    m = integral(@(x) this.pdf(x).*exp(1i*n*x), 0, 2*pi);
                case 2
                    f1 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*exp(1i*n*x(:)'), size(x));
                    f2 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*exp(1i*n*y(:)'), size(x));
                    m(1,1) = integral2(f1, 0, 2*pi, 0, 2*pi);
                    m(2,1) = integral2(f2, 0, 2*pi, 0, 2*pi);
                case 3
                    f1 = @(x,y,z) reshape(this.pdf([x(:)';y(:)';z(:)']).*exp(1i*n*x(:)'), size(x));
                    f2 = @(x,y,z) reshape(this.pdf([x(:)';y(:)';z(:)']).*exp(1i*n*y(:)'), size(x));
                    f3 = @(x,y,z) reshape(this.pdf([x(:)';y(:)';z(:)']).*exp(1i*n*z(:)'), size(x));
                    m(1,1) = integral3(f1, 0, 2*pi, 0, 2*pi, 0, 2*pi);
                    m(2,1) = integral3(f2, 0, 2*pi, 0, 2*pi, 0, 2*pi);
                    m(3,1) = integral3(f3, 0, 2*pi, 0, 2*pi, 0, 2*pi);
                otherwise
                    error('Numerical moment calculation for this dimension is currently not supported');
            end
        end
        
        function dist = toCircular(this)
            % Converts to circular distribution. Can be improved by
            % overwriting it.
            arguments
                this AbstractHypertoroidalDistribution
            end
            assert(this.dim==1, 'Can only convert distributions of dimension 1.')
            dist = CustomCircularDistribution(@(x)this.pdf(x));
        end
        
        function e = entropy(this)
            % Calculates the entropy analytically if possible, 
            % fall back to numerical calculation by default
            %
            % Returns:
            %   e (scalar)
            %       entropy of the distribution
            e = this.entropyNumerical();
        end
        
        function e = entropyNumerical(this)
            % Calculates the entropy numerically
            %
            % Returns:
            %   e (scalar)
            %       entropy of the distribution            
            switch this.dim
                case 1
                    e = -integral(@(x) this.pdf(x).*log(this.pdf(x)), 0, 2*pi);
                case 2
                    f = @(x,y) reshape(this.pdf([x(:)';y(:)']).*log(this.pdf([x(:)';y(:)'])) , size(x));
                    e = -integral2(f, 0, 2*pi, 0, 2*pi);
                case 3
                    f = @(x,y,z) reshape(this.pdf([x(:)';y(:)';z(:)']).*log(this.pdf([x(:)';y(:)';z(:)'])), size(x));
                    e = -integral3(f, 0, 2*pi, 0, 2*pi, 0, 2*pi);
                otherwise
                    error('Numerical moment calculation for this dimension is currently not supported.');
            end
        end
        
        function mu = mean2dimD(this)
            % Calculates E([cos(x1), sin(x1), cos(x2), sin(x2), ...cos(xd), sin(xd)])
            %
            % Returns:
            %   mu (2*dim x 1)
            m = this.trigonometricMoment(1);
            mu = reshape([real(m)';imag(m)'],[],1);
        end
        
        function l = logLikelihood(this, samples, weights)
            % Calculates the log-likelihood of the given samples
            %
            % Parameters:
            %   samples (dim x n)
            %       n samples on the torus
            %   weights (1 x n)
            %       weight for each sample (by default, all weights are 1)
            % Returns:
            %   l (scalar)
            %       log-likelihood of obtaining the given samples
            assert(size(samples,1)==this.dim);
            assert(size(samples,2)>=1);
            
            if nargin>2
                assert(size(weights,1)==1);
                assert(size(weights,2)==size(samples,2));
                l = sum(weights.*log(this.pdf(samples)));
            else
                l = sum(log(this.pdf(samples)));
            end
        end
        
        function s = sampleMetropolisHastings(this, n, proposal, startPoint, burnIn, skipping)
            arguments
                this (1,1) AbstractDistribution
                n (1,1) {mustBePositive,mustBeInteger}
                % Proposal is HypercylindricalWN with zero mean and
                % identity matrix. Because there is no correlation, using randn
                % can be used.
                proposal (1,1) function_handle = @(x) mod(x + randn(this.dim,1),2*pi)
                startPoint = this.meanDirection()
                burnIn (1,1) double = 10
                skipping (1,1) double = 5
            end
            s = sampleMetropolisHastings@AbstractDistribution(this, n, proposal, startPoint, burnIn, skipping);
        end
        
        function hd = shift(this, shiftAngles)
            % Shifts the distribution by shiftAngles. If not overloaded,
            % this will return a CustomHypertoroidalDistribution.
            %
            % Parameters:
            %   shiftAngles (dim x 1 column vector)
            %       angles to shift by
            % Return:
            %   hd (CustomHypertoroidalDistribution)
            %       shifted distribution
            arguments
                this (1,1) AbstractHypertoroidalDistribution
                shiftAngles (:,1) double
            end
            assert(size(shiftAngles,1)==this.dim);
            hd = CustomHypertoroidalDistribution(@(xa) this.pdf(xa), this.dim);
            hd.shiftBy = shiftAngles;
        end
        
        function d = squaredDistanceNumerical(this, other)
            % Numerically calculates the squared distance to another distribution
            %
            % Parameters:
            %   other (AbstractCircularDistribution)
            %       distribution to compare with
            % Returns:
            %   d (scalar)
            %       integral over [0,2pi)^dim of squared difference of pdfs
            assert(isa(other, 'AbstractHypertoroidalDistribution'));
            assert(this.dim==other.dim, 'Cannot compare distributions with different number of dimensions');
            
            switch this.dim
                case 1
                    d=integral(@(x) (this.pdf(x)-other.pdf(x)).^2, 0, 2*pi);
                case 2
                    d=integral2(...
                        @(x,y)reshape((this.pdf([x(:)';y(:)'])-other.pdf([x(:)';y(:)'])).^2,size(x)),...
                        0,2*pi,0,2*pi);
                case 3
                    d=integral3(...
                        @(x,y,z)reshape((this.pdf([x(:)';y(:)';z(:)'])-other.pdf([x(:)';y(:)';z(:)'])).^2,size(x)),...
                        0,2*pi,0,2*pi,0,2*pi);
                otherwise
                    error('Numerical calculation of squared distance is currently not supported for this dimension.')
            end
        end
        
        function kld = kldNumerical(this, other)
            % Numerically calculates  the Kullback-Leibler divergence to another distribution
            % Be aware that the kld is not symmetric.
            %
            % Parameters:
            %   other (AbstractHypertoroidalDistribution)
            %       distribution to compare with
            % Returns:
            %   kld (scalar)
            %       kld of this distribution to other distribution
            assert(isa(other, 'AbstractHypertoroidalDistribution'));
            assert(this.dim==other.dim, 'Cannot compare distributions with different number of dimensions');
            
            switch this.dim
                case 1
                    kld=integral(@(x) this.pdf(x).*log(this.pdf(x)./other.pdf(x)), 0, 2*pi);
                case 2
                    kld=integral2(...
                        @(x,y)reshape(this.pdf([x(:)';y(:)']).*log(this.pdf([x(:)';y(:)'])./other.pdf([x(:)';y(:)'])),size(x)),...
                        0,2*pi,0,2*pi);
                case 3
                    kld=integral3(...
                        @(x,y,z)reshape(this.pdf([x(:)';y(:)';z(:)']).*log(this.pdf([x(:)';y(:)';z(:)'])./other.pdf([x(:)';y(:)';z(:)'])),size(x)),...
                        0,2*pi,0,2*pi,0,2*pi);
                otherwise
                    error('Numerical calculation of kld is currently not supported for this dimension.')
            end
        end
        
        function dist = hellingerDistanceNumerical(this, other)
            % Numerically calculates the Hellinger distance to another
            % distribution.
            %
            % Parameters:
            %   other (AbstractHypertoroidalDistribution)
            %       distribution to compare with
            % Returns:
            %   dist (scalar)
            %       hellinger distance of this distribution to other distribution
            assert(isa(other, 'AbstractHypertoroidalDistribution'));
            assert(this.dim==other.dim,'Cannot compare distributions with different number of dimensions');
            
            switch this.dim
                case 1
                    dist=0.5*integral(@(x)(sqrt(this.pdf(x))-sqrt(other.pdf(x))).^2,0,2*pi);
                case 2
                    dist=0.5*integral2(...
                        @(x,y)reshape((sqrt(this.pdf([x(:)';y(:)']))-sqrt(other.pdf([x(:)';y(:)']))).^2,size(x)),...
                        0,2*pi,0,2*pi);
                case 3
                    dist=0.5*integral3(...
                        @(x,y,z)reshape((sqrt(this.pdf([x(:)';y(:)';z(:)']))-sqrt(other.pdf([x(:)';y(:)';z(:)']))).^2,size(x)),...
                        0,2*pi,0,2*pi,0,2*pi);
                otherwise
                    error('Numerical calculation of Hellinger distance is currently not supported for this dimension.')
            end
        end
        
        function dist = totalVariationDistanceNumerical(this, other)
            % Numerically calculates the total varation distance to another distribution
            %
            % Parameters:
            %   other (AbstractHypertoroidalDistribution)
            %       distribution to compare with
            % Returns:
            %   dist (scalar)
            %       total variation distance of this distribution to other distribution
            assert(isa(other, 'AbstractHypertoroidalDistribution'));
            assert(this.dim==other.dim, 'Cannot compare distributions with different number of dimensions');
            
            switch this.dim
                case 1
                    dist=0.5*integral(@(x)abs(this.pdf(x)-other.pdf(x)),0,2*pi);
                case 2
                    dist=0.5*integral2(...
                        @(x,y)reshape(abs(this.pdf([x(:)';y(:)'])-other.pdf([x(:)';y(:)'])),size(x)),...
                        0,2*pi,0,2*pi);
                case 3
                    dist=0.5*integral3(...
                        @(x,y,z)reshape(abs(this.pdf([x(:)';y(:)';z(:)'])-other.pdf([x(:)';y(:)';z(:)'])),size(x)),...
                        0,2*pi,0,2*pi,0,2*pi);
                otherwise
                    error('Numerical calculation of total variation distance is currently not supported for this dimension.')
            end
        end
        
        function s = getManifoldSize(this)
            s = (2*pi)^(this.dim);
        end
        
    end
end
