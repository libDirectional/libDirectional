% The Bingham Distribution.
% This class represents a d-dimensional Bingham distribution.
%
%   Notation:
%   In this class, d represents the dimension of the distribution.
%   Currently, there is no support for uniform Bingham distributions.
%
% see
% C. Bingham, "An antipodally symmetric distribution on the sphere",
% The Annals of Statistics, vol. 2, no. 6, pp. 1201-1225, Nov. 1974.

classdef BinghamDistribution < AbstractHypersphericalDistribution
    
    properties
        Z (:,1) double {mustBeNonpositive}       % Concentrations as a vector
        M (:,:) double  % Rotation matrix
        F (1,1) double     % Normalization constant
        dF (1,:) double     % Partial derivates of F
    end
    
    properties (Constant)
        S2 = AbstractHypersphericalDistribution.computeUnitSphereSurface(2) % Circle length
    end
    
    methods
        function B = BinghamDistribution(Z_, M_)
            % Constructs a Bingham distribution object.
            % Parameters:
            %   Z_ (d x 1 column vector)
            %       concentration parameters (have to be increasing, and last
            %       entry has to be zero)
            %   M_ (d x d matrix)
            %       orthogonal matrix that describes rotation
            % Returns:
            %   B (BinghamDistribution)
            %       an object representing the constructed distribution
            arguments
                Z_ (:,1) double
                M_ (:,:) double
            end
            B.dim = size(M_,1);
            
            % Check Dimensions
            assert(size(M_,2) == B.dim, 'M is not square');
            assert(size(Z_,1) == B.dim, 'Z has wrong number of rows');
            assert(size(Z_,2) == 1, 'Z needs to be column vector');
            
            % Enforce last entry of Z to be zero
            assert(Z_(B.dim, 1) == 0, 'last entry of Z needs to be zero');
            
            % Enforce z1<=z2<=...<=z(d-1)<=0=z(d)
            assert(all(Z_(1:end-1) <= Z_(2:end)), 'values in Z have to be ascending');
            
            %enforce that M is orthogonal
            epsilon = 0.001;
            assert (max(max(M_*M_' - eye(B.dim,B.dim))) < epsilon, 'M is not orthogonal');
            
            B.Z = Z_;
            B.M = M_;
            B.F = BinghamDistribution.computeF(B.Z);
            B.dF = BinghamDistribution.computeDF(B.Z);
        end
        
        function p = pdf(this, xa)
            % Evaluates pdf at each column of xa
            % Parameters:
            %   xa (d x n matrix)
            %       each column represents one of the n points in R^d that the
            %       pdf is evaluated at; can be just a (d x 1) vector as well
            % Returns:
            %   p (1 x n row vector)
            %       values of the pdf at each column of xa
            arguments
                this (1,1) BinghamDistribution
                xa (:,:) double
            end
            assert(size(xa,1) == this.dim);
            
            C = this.M * diag(this.Z) * this.M';
            p = 1/this.F * exp(sum(xa.*(C*xa)));
        end
        
        function meanDirection(~)
            error('Due to their symmetry, the mean direction is undefined for Bingham distributions.');
        end
        
        function B = multiply(this, B2)
            % Computes the product of two Bingham pdfs
            % This method makes use of the fact that the Bingham distribution
            % is closed under Bayesian inference. Thus, the product of two
            % Bingham pdfs is itself the pdf of a Bingham distribution. This
            % method computes the parameters of the resulting distribution.
            %
            % Parameters:
            %   B2 (BinghamDistribution)
            %       second Bingham Distribution
            % Returns:
            %   B (BinghamDistribution)
            %       Bingham distribution representing this*B2 (after
            %       renormalization)
            arguments
                this (1,1) BinghamDistribution
                B2 (1,1) BinghamDistribution
            end
            assert(isa(B2, 'BinghamDistribution'));
            if this.dim == B2.dim
                C = this.M * diag(this.Z) * this.M' + B2.M * diag(B2.Z) * B2.M'; % new exponent
                
                C = 0.5*(C+C'); % Ensure symmetry of C, asymmetry may arise as a consequence of a numerical instability earlier.
                
                [V,D] = eig(C); % Eigenvalue decomposition
                [D, order] = sort(diag(D),'ascend');  % sort eigenvalues
                V = V(:,order);
                Z_ = D;
                Z_ = Z_-Z_(end); % last entry should be zero
                M_ = V;
                B = BinghamDistribution(Z_,M_);
            else
                error('dimensions do not match');
            end
        end
        
        function B = compose(this, B2)
            % Compose two Bingham distributions
            % Using Moment Matching based approximation, we compose two Bingham
            % distributions. The mode of the new distribution should be the
            % quaternion multiplication of the original modes; the uncertainty
            % should be larger than before
            %
            % Parameters:
            %   B2 (BinghamDistribution)
            %       second Bingham Distribution
            % Returns:
            %   B (BinghamDistribution)
            %       Bingham distribution representing the convolution
            arguments
                this (1,1) BinghamDistribution
                B2 (1,1) BinghamDistribution
            end
            B1=this;
            B1S = B1.moment();
            B2S = B2.moment();
            if this.dim==2 && B2.dim==2
                % for complex numbers
                % derived from complex multiplication
                % Gerhard Kurz, Igor Gilitschenski, Simon Julier, Uwe D. Hanebeck,
                % Recursive Bingham Filter for Directional Estimation Involving 180 Degree Symmetry
                % Journal of Advances in Information Fusion, 9(2):90 - 105, December 2014.
                a11 = B1S(1,1);
                a12 = B1S(1,2);
                a22 = B1S(2,2);
                b11 = B2S(1,1);
                b12 = B2S(1,2);
                b22 = B2S(2,2);
                
                S(1,1) = a11*b11 - 2*a12*b12 + a22*b22;
                S(1,2) = a11*b12 - a22*b12 - a12*b22 + a12*b11;
                S(2,1) = S(1,2);
                S(2,2) = a11*b22 + 2*a12*b12 + a22*b11;
                
                B = BinghamDistribution.fitToMoment(S);
            elseif this.dim==4 && B2.dim==4
                % adapted from Glover's C code in libBingham, see also
                % Glover, J. & Kaelbling, L. P.
                % Tracking 3-D Rotations with the Quaternion Bingham Filter
                % MIT, 2013
                
                a11 = B1S(1,1);
                a12 = B1S(1,2);
                a13 = B1S(1,3);
                a14 = B1S(1,4);
                a22 = B1S(2,2);
                a23 = B1S(2,3);
                a24 = B1S(2,4);
                a33 = B1S(3,3);
                a34 = B1S(3,4);
                a44 = B1S(4,4);
                
                b11 = B2S(1,1);
                b12 = B2S(1,2);
                b13 = B2S(1,3);
                b14 = B2S(1,4);
                b22 = B2S(2,2);
                b23 = B2S(2,3);
                b24 = B2S(2,4);
                b33 = B2S(3,3);
                b34 = B2S(3,4);
                b44 = B2S(4,4);
                
                %can be derived from quaternion multiplication
                S(1,1) = a11*b11 - 2*a12*b12 - 2*a13*b13 - 2*a14*b14 + a22*b22 + 2*a23*b23 + 2*a24*b24 + a33*b33 + 2*a34*b34 + a44*b44;
                S(1,2) = a11*b12 + a12*b11 + a13*b14 - a14*b13 - a12*b22 - a22*b12 - a13*b23 - a23*b13 - a14*b24 - a24*b14 - a23*b24 + a24*b23 - a33*b34 + a34*b33 - a34*b44 + a44*b34;
                S(2,1) = S(1,2);
                S(1,3) = a11*b13 + a13*b11 - a12*b14 + a14*b12 - a12*b23 - a23*b12 - a13*b33 + a22*b24 - a24*b22 - a33*b13 - a14*b34 - a34*b14 + a23*b34 - a34*b23 + a24*b44 - a44*b24;
                S(3,1) = S(1,3);
                S(1,4) = a11*b14 + a12*b13 - a13*b12 + a14*b11 - a12*b24 - a24*b12 - a22*b23 + a23*b22 - a13*b34 - a34*b13 - a23*b33 + a33*b23 - a14*b44 - a24*b34 + a34*b24 - a44*b14;
                S(4,1) = S(1,4);
                S(2,2) = 2*a12*b12 + a11*b22 + a22*b11 + 2*a13*b24 - 2*a14*b23 + 2*a23*b14 - 2*a24*b13 - 2*a34*b34 + a33*b44 + a44*b33;
                S(2,3) = a12*b13 + a13*b12 + a11*b23 + a23*b11 - a12*b24 + a14*b22 - a22*b14 + a24*b12 + a13*b34 - a14*b33 + a33*b14 - a34*b13 + a24*b34 + a34*b24 - a23*b44 - a44*b23;
                S(3,2) = S(2,3);
                S(2,4) = a12*b14 + a14*b12 + a11*b24 + a12*b23 - a13*b22 + a22*b13 - a23*b12 + a24*b11 - a14*b34 + a34*b14 + a13*b44 + a23*b34 - a24*b33 - a33*b24 + a34*b23 - a44*b13;
                S(4,2) = S(2,4);
                S(3,3) = 2*a13*b13 + 2*a14*b23 - 2*a23*b14 + a11*b33 + a33*b11 - 2*a12*b34 + 2*a34*b12 - 2*a24*b24 + a22*b44 + a44*b22;
                S(3,4) = a13*b14 + a14*b13 - a13*b23 + a23*b13 + a14*b24 - a24*b14 + a11*b34 + a12*b33 - a33*b12 + a34*b11 + a23*b24 + a24*b23 - a12*b44 - a22*b34 - a34*b22 + a44*b12;
                S(4,3) = S(3,4);
                S(4,4) = 2*a14*b14 - 2*a13*b24 + 2*a24*b13 + 2*a12*b34 - 2*a23*b23 - 2*a34*b12 + a11*b44 + a22*b33 + a33*b22 + a44*b11;
                
                B = BinghamDistribution.fitToMoment(S);
            else
                error('unsupported dimension');
            end
        end
        
        function s = sample(this, n)
            % Stocahastic sampling
            % Fall back to Kent's method by default
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       one sample per column
            arguments
                this (1,1) BinghamDistribution
                n (1,1) {mustBeInteger, mustBePositive}
            end
            s = sampleKent(this, n);
        end
        
        function s = sampleKent(this, n)
            % Generate samples from Bingham distribution using rejection
            % sampling based on a angular central Gaussian
            %
            % Kent, J. T.; Ganeiber, A. M. & Mardia, K. V.
            % A New Method to Simulate the Bingham and Related Distributions in Directional Data Analysis with Applications
            % arXiv preprint arXiv:1310.8110, 2013
            arguments
                this (1,1) BinghamDistribution
                n (1,1) {mustBeInteger, mustBePositive}
            end
            s = zeros(this.dim, n);
            i = 1;
            A = - this.M * diag(this.Z) * this.M'; % Kent uses a minus sign here!
            q = this.dim;
            
            % compute b
            bfun = @(b) sum(1./(b-2*this.Z)) - 1; % use a minus sign before  2*this.z because Kent's matrix is negative
            b = fsolve(bfun, 1, optimset('display', 'none'));
            Omega = eye(this.dim) + 2*A/b;
            %efficiency = 1/(exp(-(q-b)/2) * (q/b)^(q/2))
            %Mb = 1/this.F * exp(-(q-b)/2) * (q/b)^(q/2) * det(Omega)^(-1/2);
            Mbstar = exp(-(q-b)/2) * (q/b)^(q/2);
            nReject = 0;
            fbingstar = @(x) exp(-x' * A * x);
            facgstar = @(x) (x' * Omega * x)^(-q/2);
            while(i<=n)
                % draw x from angular central Gaussian
                y = mvnrnd(zeros(this.dim,1), inv(Omega))';
                x = y/norm(y);
                % check rejection
                W = rand(1);
                if W < fbingstar(x) /(Mbstar * facgstar(x))
                    s(:,i) = x;
                    i = i + 1;
                else
                    nReject = nReject + 1;
                end
            end
            %nReject
        end
        
        function X = sampleGlover(this, n)
            % Generate samples from Bingham distribution
            % based on Glover's implementation in libBingham
            % uses Metropolis-Hastings
            % see http://en.wikipedia.org/wiki/Metropolis-Hastings_algorithm
            %
            % The implementation has a bug because it just repeats the
            % previos sample if a sample is rejected.
            %
            % Parameters:
            %   n (scalar)
            %       number of samples to generate
            % Returns:
            %   X (dimx n matrix)
            %       generated samples (one sample per column)
            arguments
                this (1,1) BinghamDistribution
                n (1,1) {mustBeInteger, mustBePositive}
            end
            burnin = 5;
            samplerate = 10;
            
            x = this.mode();
            z = sqrt(-1./(this.Z - 1));
            
            target = this.pdf(x);  % target
            proposal = acgpdf_pcs(x', z, this.M);  % proposal
            
            X2 = (randn(n*samplerate+burnin,this.dim).*repmat(z',[n*samplerate+burnin,1]))*this.M'; % sample Gaussian
            X2 = X2 ./ repmat(sqrt(sum(X2.^2,2)), [1 this.dim]); % normalize
            
            Target2 = this.pdf(X2');
            Proposal2 = acgpdf_pcs(X2, z, this.M);
            
            nAccepts = 0;
            X = zeros(size(X2));
            for i=1:n*samplerate+burnin
                a = Target2(i) / target * proposal / Proposal2(i);
                if a > rand()
                    x = X2(i,:);
                    proposal = Proposal2(i);
                    target =  Target2(i);
                    nAccepts = nAccepts + 1;
                end
                X(i,:) = x;
            end
            
            %accept_rate = num_accepts / (n*sample_rate + burn_in)
            
            X = X(burnin+1:samplerate:end,:)';
        end
        
        function m = mode(this)
            % Calculate the mode of a Bingham distribution
            % Returns:
            %   m (column vector)
            %       mode of the distribution (note that -m is the mode as well)
            arguments
                this (1,1) BinghamDistribution
            end
            m = this.M(:,end); %last column of M
        end
        
        function S = moment(this)
            arguments
                this (1,1) BinghamDistribution
            end
            % Calculate scatter/covariance matrix of Bingham distribution
            % Returns:
            %   S (d x d matrix)
            %     scatter/covariance matrix in R^d
            D = diag(this.dF/this.F);
            % the sum of the diagonal of D is always 1, however this may
            % not be the case because dF and F are calculated using
            % approximations
            D = D / sum(diag(D));
            S = this.M * D * this.M';
            S = (S+S')/2; % enforce symmetry
        end
        
        function [samples, weights] = sampleDeterministic(this, lambda)
            % Computes deterministic samples of a Bingham distribution.
            % The computed samples represent the current Bingham distribution.
            % They are choosen in a deterministic way, which is a circular
            % adaption of the UKF.
            %
            % see
            % Igor Gilitschenski, Gerhard Kurz, Simon J. Julier, Uwe D. Hanebeck,
            % Unscented Orientation Estimation Based on the Bingham Distribution
            % IEEE Transactions on Automatic Control, January 2016.
            %
            % Important: The paper discusses samples from both modes of the
            % Bingham. This code only samples from one mode. If desired, the samples
            % from the other mode can be obtained by mirroring:
            % [s,w] = bd.sampleDeterministic
            % s2 = [s, -s]
            % w2 = [w, w]/2
            %
            % Parameters:
            %   lambda (scalar or string)
            %       weighting parameter in [0,1] or the string 'uniform' (only
            %       in 2D)
            % Returns:
            %   samples (d x ...  matrix)
            %     generated samples (one sample per column)
            %   weights (1 x ... vector)
            %     weight > 0 for each sample (weights sum to one)
            arguments
                this (1,1) BinghamDistribution
                lambda = []
            end
            if isempty(lambda) % default value for lambda
                if this.dim == 2
                    lambda = 'uniform';
                else
                    lambda = 0.5;
                end
            end
            if strcmp(lambda,'uniform') && this.dim == 2
                % uniform weights can only be guaranteed for d=2
                B = BinghamDistribution(this.Z, eye(2,2));
                S = B.moment();
                alpha = asin(sqrt(S(1,1)*3/2));
                samples = [ 0, 1;
                    sin(alpha),  cos(alpha);
                    -sin(alpha), cos(alpha)];
                samples = this.M*samples';
                weights = [1/3, 1/3, 1/3];
            else
                assert (lambda>=0 && lambda <=1);
                B = BinghamDistribution(this.Z, eye(this.dim,this.dim));
                S = B.moment();
                samples = zeros(2*this.dim-1,this.dim);
                weights = zeros(1, 2*this.dim-1);
                p = zeros(1, this.dim-1);
                alpha = zeros(1, this.dim-1);
                samples(1,end) = 1; %sample at mode
                for i=1:this.dim-1
                    p(i) = S(i,i) + (1-lambda)*(S(end,end)/(this.dim-1));
                    alpha(i) = asin(sqrt(S(i,i)/p(i)));
                    samples(2*i,end) = cos(alpha(i));
                    samples(2*i+1,end) = cos(alpha(i));
                    samples(2*i,i) = sin(alpha(i));
                    samples(2*i+1,i) = -sin(alpha(i));
                    weights(1,2*i) = p(i)/2;
                    weights(1,2*i+1) = p(i)/2;
                end
                weights(1) = 1-sum(weights(2:end)); % = lambda*S(4,4)
                samples = this.M*samples';
            end
        end
        
        function [samples, weights] = sampleOptimalQuantization(this, N)
            % Computes optimal quantization of the
            % Bingham distribution using 2*N samples.
            %
            % Parameters:
            %   N(scalar)
            %       number of samples on half circle
            % Returns:
            %   samples (1 x 2N)
            %       2N samples on the circle, parameterized as [0,2pi)
            %   weights (1 x 2N)
            %       weight for each sample
            %
            % Igor Gilitschenski, Gerhard Kurz, Uwe D. Hanebeck, Roland Siegwart,
            % Optimal Quantization of Circular Distributions
            % Proceedings of the 19th International Conference on Information Fusion (Fusion 2016), Heidelberg, Germany, July 2016.
            arguments
                this (1,1) BinghamDistribution
                N (1,1) {mustBeInteger, mustBePositive}
            end
            assert(this.dim == 2, 'sampleOptimalQuantization only implemented for 2d Bingham');
            
            mu = atan2(this.M(2,2), this.M(1,2));
            kappa = (this.Z(2)-this.Z(1))/2;
            
            [samples, weights] = VMDistribution(0, kappa).sampleOptimalQuantization(N);
            samples = [samples/2 (samples/2+pi)];
            samples = mod(samples+mu,2*pi);
            
            weights = [weights weights]/2;
        end
        
        function [s,w] = sampleWeighted(this, n)
            % Weighted sample generator.
            % Generates uniform (w.r.t. the Haar Measure) random samples on
            % the unit sphere and assigns each sample a weight based on the
            % pdf.
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            %
            % Returns:
            %   samples (d x n  matrix)
            %     generated samples (one sample per column)
            %   weights (1 x n vector)
            arguments
                this (1,1) BinghamDistribution
                n (1,1) {mustBeInteger, mustBePositive}
            end
            s = mvnrnd(zeros(1,this.dim), eye(this.dim), n)';
            s = s./repmat(sqrt(sum(s.^2,1)),this.dim, 1);
            
            w = this.pdf(s);
            w = w/sum(w); % normalize weights
        end
        
        function alpha = bounds(this, p)
            % Calculates the confidence interval for a given confidence p
            % Only works for d=2 so far.
            % Parameters:
            %   p (scalar)
            %       confidence to achieve, 0<p<1
            % Returns:
            %   alpha (scalar)
            %       the confidence interval ranges from mode-alpha to
            %       mode+alpha
            arguments
                this (1,1) BinghamDistribution
                p (1,1) {mustBeInRange(p,0,1,'exclusive')}
            end
            if this.dim==2
                assert(p>0)
                assert(p<1)
                f = @(phi) this.pdf([cos(phi); sin(phi)]);
                mB = this.mode();
                mBangle = atan2(mB(2), mB(1));
                g = @(alpha) integral(f,mBangle,mBangle+alpha, 'AbsTol', 0.01) - p/2;
                alpha = fsolve(g, 0.1, optimset('display', 'none'));
            else
                error('unsupported dimension');
            end
        end
        
        function P = gaussianCovariance(this, angle)
            % Calculates the covariance to obtain a comparable Gaussian.
            % Attention: this does not compute the covariance of the entire
            % Bingham but only on the half sphere!
            %
            % Returns:
            %   P (1 x 1 matrix for d=2, if angle = true (angle represenation)
            %      d x d matrix for d>=2, if angle = false (vector
            %      representation)
            arguments
                this (1,1) BinghamDistribution
                angle = false
            end
            if angle
                assert(this.dim==2)
                m = this.mode();
                mAngle = atan2(m(2),m(1));
                
                % sample-based solution:
                % samples = this.sample(1000);
                % sampleAngles = atan2(samples(:,2), samples(:,1));
                % sampleAngles = mod(sampleAngles + pi/2 + mAngle, pi);
                % P = cov(sampleAngles);
                
                % numerical-integration-based solution
                P = 2*integral(@(phi) this.pdf([cos(mAngle+phi);sin(mAngle+phi)]).*phi.^2, -pi/2, +pi/2);
            else
                if this.dim==2
                    % numerical-integration-based solution:
                    m = this.mode();
                    mAngle = atan2(m(2),m(1));
                    fx  = @(phi) 2*this.pdf([cos(phi);sin(phi)]).*(cos(phi)-m(1)).^2;
                    fxy = @(phi) 2*this.pdf([cos(phi);sin(phi)]).*(cos(phi)-m(1)).*(sin(phi)-m(2));
                    fy  = @(phi) 2*this.pdf([cos(phi);sin(phi)]).*(sin(phi)-m(2)).^2;
                    P(1,1) = integral(fx, mAngle-pi/2, mAngle+pi/2);
                    P(1,2) = integral(fxy, mAngle-pi/2, mAngle+pi/2);
                    P(2,1) = P(1,2);
                    P(2,2) = integral(fy, mAngle-pi/2, mAngle+pi/2);
                else
                    count = 1000;
                    samples = this.sample(count);
                    for i=1:count
                        if samples(:,i)'*this.mode()<0 %scalar product
                            samples(:,i)=-samples(:,i);
                        end
                    end
                    %P = cov(samples);
                    % do not use cov, because the mean of the Gaussian will
                    % be at the mode of the Bingham, not at the mean of the
                    % samples
                    samples = samples - repmat(this.mode(), 1, count);
                    P = samples*samples'/count;
                end
            end
        end
    end
    
    methods (Static)
        function F = computeF(Z, mode)
            % Compute normalization constant
            % Parameters:
            %   Z (d x d matrix)
            %       concentration matrix
            %   mode (string, optional)
            %       choses the algorithm to compute the normalization constant
            % Returns:
            %   F (scalar)
            %       the calculated normalization constant
            assert(size(Z,2) == 1);
            dim = length(Z);
            
            if nargin<2
                mode = 'default';
            end
            
            if dim == 2
                if strcmp(mode, 'default') || strcmp(mode, 'bessel')
                    % Gerhard Kurz, Igor Gilitschenski, Simon Julier, Uwe D. Hanebeck,
                    % Recursive Bingham Filter for Directional Estimation Involving 180 Degree Symmetry
                    % Journal of Advances in Information Fusion, 9(2):90 - 105, December 2014.
                    F = exp(Z(2))* BinghamDistribution.S2 * besseli(0,(Z(1)-Z(2))/2) * exp((Z(1)-Z(2))/2);
                elseif strcmp(mode, 'hypergeom')
                    % Gerhard Kurz, Igor Gilitschenski, Simon J. Julier, Uwe D. Hanebeck,
                    % Recursive Estimation of Orientation Based on the Bingham Distribution
                    % Proceedings of the 16th International Conference on Information Fusion (Fusion 2013), Istanbul, Turkey, July 2013.
                    F = exp(Z(2))* BinghamDistribution.S2 * double(hypergeom(0.5,1, vpa(Z(1)-Z(2))));
                elseif strcmp(mode, 'mhg')
                    % Koev, P. & Edelman, A.
                    % The Efficient Evaluation of the Hypergeometric Function of a Matrix Argument
                    % Mathematics of Computation., 2006, 75, 833-846
                    F = BinghamDistribution.S2 * mhg(100, 2, 0.5, dim/2, Z);
                elseif strcmp(mode, 'saddlepoint')
                    % Igor Gilitschenski, Gerhard Kurz, Simon J. Julier, Uwe D. Hanebeck,
                    % Efficient Bingham Filtering based on Saddlepoint Approximations
                    % Proceedings of the 2014 IEEE International Conference on Multisensor Fusion and Information Integration (MFI 2014), Beijing, China, September 2014.
                    F = numericalSaddlepointWithDerivatives(sort(-Z)+1)*exp(1);
                    F = F(3);
                elseif strcmp(mode, 'glover')
                    % https://code.google.com/p/bingham/
                    % and
                    %
                    % Glover, J. & Kaelbling, L. P.
                    % Tracking 3-D Rotations with the Quaternion Bingham Filter
                    % MIT, 2013
                    % http://dspace.mit.edu/handle/1721.1/78248
                    
                    F = glover(Z);
                else
                    error('unsupported mode');
                end
            else
                if strcmp(mode, 'default') || strcmp(mode, 'saddlepoint')
                    % Igor Gilitschenski, Gerhard Kurz, Simon J. Julier, Uwe D. Hanebeck,
                    % Efficient Bingham Filtering based on Saddlepoint Approximations
                    % Proceedings of the 2014 IEEE International Conference on Multisensor Fusion and Information Integration (MFI 2014), Beijing, China, September 2014.
                    F = numericalSaddlepointWithDerivatives(sort(-Z)+1)*exp(1);
                    F = F(3);
                elseif strcmp(mode, 'mhg')
                    % Koev, P. & Edelman, A.
                    % The Efficient Evaluation of the Hypergeometric Function of a Matrix Argument
                    % Mathematics of Computation., 2006, 75, 833-846
                    F = AbstractHypersphericalDistribution.computeUnitSphereSurface(dim) * mhg(100, 2, 0.5, dim/2, Z);
                elseif strcmp(mode, 'wood') && dim == 4
                    % ANDREW T.A. WOOD
                    % ESTIMATION OF THE CONCENTRATION PARAMETERS
                    % OF THE FISHER MATRIX DISTRIBUTION ON SO(3)
                    % AND THE BINGHAM DISTRIBUTION ON Sq, q>= 2
                    % Austral. J. Statist., S5(1), 1993, 69-79
                    J = @(Z,u) besseli(0, 0.5 .* abs(Z(1)-Z(2)) .* u) .* besseli(0, 0.5 .* abs(Z(3)-Z(4)) .* (1-u));
                    ifun = @(u) J(Z,u).*exp(0.5 .* (Z(1)+Z(2)).* u + 0.5.*(Z(3)+Z(4)).*(1-u));
                    F = 2*pi^2*integral(ifun,0,1);
                elseif strcmp(mode, 'glover') && dim <= 4
                    % https://code.google.com/p/bingham/
                    % and
                    %
                    % Glover, J. & Kaelbling, L. P.
                    % Tracking 3-D Rotations with the Quaternion Bingham Filter
                    % MIT, 2013
                    % http://dspace.mit.edu/handle/1721.1/78248
                    
                    if any(abs(Z(1:end-1))<1e-8)
                        warning('Glover''s method currently does not work for Z with more than one zero entry.');
                    end
                    
                    F = glover(Z);
                else
                    error('unsupported mode');
                end
            end
        end
        
        %todo: add glover?
        
        function dF = computeDF(Z, mode)
            % Partial derivatives of normalization constant
            % Parameters:
            %   Z (d x d matrix)
            %       concentration matrix
            % Returns:
            %   dF (scalar)
            %       the calculated derivative of the normalization constant
            assert(size(Z,2) == 1);
            
            dim = size(Z,1);
            dF = zeros(1,dim);
            
            if nargin<2
                mode = 'default';
            end
            
            if dim==2
                if strcmp(mode, 'default') || strcmp(mode, 'bessel')
                    % Gerhard Kurz, Igor Gilitschenski, Simon Julier, Uwe D. Hanebeck,
                    % Recursive Bingham Filter for Directional Estimation Involving 180 Degree Symmetry
                    % Journal of Advances in Information Fusion, 9(2):90 - 105, December 2014.
                    b1 = besseli(1,(Z(1)-Z(2))/2);
                    b0 = besseli(0,(Z(1)-Z(2))/2);
                    dF(1) = BinghamDistribution.S2/2 * (b1 + b0)* exp((Z(1)+Z(2))/2);
                    dF(2) = BinghamDistribution.S2/2 * (-b1 + b0 )* exp((Z(1)+Z(2))/2);
                elseif strcmp(mode, 'hypergeom')
                    % Gerhard Kurz, Igor Gilitschenski, Simon J. Julier, Uwe D. Hanebeck,
                    % Recursive Estimation of Orientation Based on the Bingham Distribution
                    % Proceedings of the 16th International Conference on Information Fusion (Fusion 2013), Istanbul, Turkey, July 2013.
                    h = double(hypergeom(1.5,2,vpa(Z(1)-Z(2))));
                    dF(1) = BinghamDistribution.S2 * exp(Z(2)) * 0.5 * h;
                    dF(2) = BinghamDistribution.S2 * exp(Z(2)) * (double(hypergeom(0.5,1, vpa(Z(1)-Z(2)))) - 0.5*h);
                elseif strcmp(mode, 'saddlepoint')
                    % Igor Gilitschenski, Gerhard Kurz, Simon J. Julier, Uwe D. Hanebeck,
                    % Efficient Bingham Filtering based on Saddlepoint Approximations
                    % Proceedings of the 2014 IEEE International Conference on Multisensor Fusion and Information Integration (MFI 2014), Beijing, China, September 2014.
                    for i=1:dim
                        ModZ = Z([1:i i i:dim]);
                        T = numericalSaddlepointWithDerivatives(sort(-ModZ)+1)*exp(1)/(2*pi);
                        dF(i) = T(3);
                    end
                elseif strncmp(mode, 'finitedifferences', 17)
                    % Approximation by finite Differences
                    % use mode='finitedifferences-method', where method
                    % is a method for calculating the normalizaton
                    % constant
                    for i=1:dim
                        epsilon=0.001;
                        dZ = [zeros(i-1,1);  epsilon; zeros(dim-i,1)];
                        F1 = BinghamDistribution.computeF(Z + dZ, mode(19:end));
                        F2 = BinghamDistribution.computeF(Z - dZ, mode(19:end));
                        dF(i) = (F1-F2)/(2*epsilon);
                    end
                else
                    error('unsupported mode');
                end
            else
                if strcmp(mode, 'default') || strcmp(mode, 'saddlepoint')
                    % Igor Gilitschenski, Gerhard Kurz, Simon J. Julier, Uwe D. Hanebeck,
                    % Efficient Bingham Filtering based on Saddlepoint Approximations
                    % Proceedings of the 2014 IEEE International Conference on Multisensor Fusion and Information Integration (MFI 2014), Beijing, China, September 2014.
                    for i=1:dim
                        ModZ = Z([1:i i i:dim]);
                        T = numericalSaddlepointWithDerivatives(sort(-ModZ)+1)*exp(1)/(2*pi);
                        dF(i) = T(3);
                    end
                elseif strncmp(mode, 'finitedifferences', 17)
                    for i=1:dim
                        % Approximation by finite differences
                        % use mode='finitedifferences-method', where method
                        % is a method for calculating the normalizaton
                        % constant
                        epsilon=0.001;
                        dZ = [zeros(i-1,1);  epsilon; zeros(dim-i,1)];
                        F1 = BinghamDistribution.computeF(Z + dZ, mode(19:end));
                        F2 = BinghamDistribution.computeF(Z - dZ, mode(19:end));
                        dF(i) = (F1-F2)/(2*epsilon);
                    end
                else
                    error('unsupported mode');
                end
            end
        end
        
        function B = fit(samples, weights, options)
            % Fits Bingham parameters to a set of samples
            % Parameters:
            %   samples (d x n matrix)
            %       matrix that contains one sample per column
            %   weights (1 x n row vector)
            %       weight for each sample
            %   options (struct)
            %       parameter to select the MLE algorithm
            % Returns:
            %   B (BinghamDistribution)
            %       the MLE estimate for a Bingham distribution given the
            %       samples
            n = size(samples,2);
            if nargin<2
                C = samples*samples'/n;
            else
                assert(abs(sum(weights)-1) < 1E-10, 'weights must sum to 1'); %check normalization
                assert(size(weights,1)==1, 'weights needs to be a row vector');
                assert(size(weights,2)==n, 'number of weights and samples needs to match');
                C = samples.*weights*samples';
            end
            
            C = (C+C')/2; % ensure symmetry
            
            if nargin<3
                B = BinghamDistribution.fitToMoment(C);
            else
                B = BinghamDistribution.fitToMoment(C, options);
            end
        end
               
        function B = fitToMoment(S, options)
            % Finds a Bingham distribution with a given second moment
            %
            % Parameters:
            %	S (d x d matrix)
            %       matrix representing second moment.
            %   options (struct)
            %       parameters to configure the MLE algorithm
            % Returns:
            %   B (BinghamDistribution)
            %       the MLE estimate for a Bingham distribution given the
            %       scatter matrix S
            
            assert(all(all(S == S')), 'S must be symmetric');
            if nargin < 2
                options.algorithm = 'default';
            end
            
            assert(isfield(options,'algorithm'), ...
                'Options need to contain an algorithm field');
            
            if strcmp(options.algorithm, 'default') || strcmp(options.algorithm, 'fsolve')
                % Gerhard Kurz, Igor Gilitschenski, Simon Julier, Uwe D. Hanebeck,
                % Recursive Bingham Filter for Directional Estimation Involving 180 Degree Symmetry
                % Journal of Advances in Information Fusion, 9(2):90 - 105, December 2014.
                [eigenvectors,omega] = eig(S);
                [omega, order] = sort(diag(omega),'ascend');  % sort eigenvalues
                M_ = eigenvectors(:,order); % swap columns to match the sorted eigenvalues
                omega = omega / sum(omega); % ensure that entries sum to one
                Z_ = BinghamDistribution.mleFsolve(omega, options);
                
                % This reordering shouldn't be necessary. However, it can
                % become necessary as a consequence of numerical errors when
                % fitting to moment matrices with almost equal eigenvalues.
                [Z_, order] = sort(Z_,'ascend');  %sort eigenvalues
                M_ = M_(:,order);
                Z_=Z_-Z_(size(Z_,1)); %subtract last entry to ensure that last entry is zero
                B = BinghamDistribution(Z_,M_);
            elseif strcmp(options.algorithm, 'fminunc')
                [eigenvectors,omega] = eig(S);
                [omega, order] = sort(diag(omega),'ascend');  % sort eigenvalues
                M_ = eigenvectors(:,order); % swap columns to match the sorted eigenvalues
                omega = omega / sum(omega); % ensure that entries sum to one
                Z_ = BinghamDistribution.mleFminunc(omega, options);
                
                % This reordering shouldn't be necessary. However, it can
                % become necessary as a consequence of numerical errors when
                % fitting to moment matrices with almost equal eigenvalues.
                [Z_, order] = sort(Z_,'ascend');  %sort eigenvalues
                M_ = M_(:,order);
                Z_=Z_-Z_(size(Z_,1)); %subtract last entry to ensure that last entry is zero
                B = BinghamDistribution(Z_,M_);
            elseif strcmp(options.algorithm, 'gaussnewton')
                % Igor Gilitschenski, Gerhard Kurz, Simon J. Julier, Uwe D. Hanebeck,
                % Efficient Bingham Filtering based on Saddlepoint Approximations
                % Proceedings of the 2014 IEEE International Conference on Multisensor Fusion and Information Integration (MFI 2014), Beijing, China, September 2014.
                [eigenvectors,omega] = eig(S);
                [omega, order] = sort(diag(omega),'ascend');  %sort eigenvalues
                M_ = eigenvectors(:,order); %swap columns to match the sorted eigenvalues
                omega = omega / sum(omega); % ensure that entries sum to one
                Z_ = numericalBinghamMLE(omega);
                
                % This reordering shouldn't be necessary. However, it can
                % become necessary as a consequence of numerical errors when
                % fitting to moment matrices with almost equal eigenvalues.
                [Z_, order] = sort(Z_,'ascend');  %sort eigenvalues
                M_ = M_(:,order);
                Z_=Z_-Z_(size(Z_,1)); %subtract last entry to ensure that last entry is zero
                B = BinghamDistribution(Z_,M_);
            else
                error('Unsupported estimation algorithm');
            end
        end
        
        function Z = mleFsolve (omega, options)
            % Calculate maximum likelihood estimate of Z.
            % Considers only the first three values of omega.
            %
            % Parameters:
            %   omega (d x 1 column vector)
            %       eigenvalues of the scatter matrix
            % Returns:
            %   Z (d x 1 column vector)
            
            if nargin < 2 || ~isfield(options, 'Fmethod')
                options.Fmethod = 'default';
            end
            
            if nargin < 2 || ~isfield(options, 'dFmethod')
                options.dFmethod = 'default';
            end
            
            function r = mleGoalFun(z, rhs)
                % objective function of MLE.
                d = size(z,1)+1;
                
                %if d == 2
                a = BinghamDistribution.computeF([z;0], options.Fmethod);
                b = BinghamDistribution.computeDF([z;0], options.dFmethod);
                %else
                %    [a,b] = numericalSaddlepointWithDerivatives([-z; 0]);
                %    a = a(3);
                %    b = b(3,:);
                %end
                
                r = zeros(d-1,1);
                for i=1:(d-1)
                    r(i) = b(i)/a - rhs(i);
                end
            end
            
            dim = size(omega,1);
            
            f = @(z) mleGoalFun(z, omega);
            Z = fsolve(f, -ones(dim-1,1), optimset('display', 'off', 'algorithm', 'levenberg-marquardt'));
            Z = [Z; 0];
        end
        
        function Z = mleFminunc (omega, options)
            % Calculate maximum likelihood estimate of Z.
            % Considers all four values of omega.
            %
            % Parameters:
            %   omega (d x 1 column vector)
            %       eigenvalues of the scatter matrix
            % Returns:
            %   Z (d x 1 column vector)
            
            if nargin < 2 || ~isfield(options, 'Fmethod')
                options.Fmethod = 'default';
            end
            
            if nargin < 2 || ~isfield(options, 'dFmethod')
                options.dFmethod = 'default';
            end
            
            function r = mleGoalFun(z, rhs)
                % objective function of MLE.
                d = size(z,1)+1;
                
                %if d == 2
                a = BinghamDistribution.computeF([z;0], options.Fmethod);
                b = BinghamDistribution.computeDF([z;0], options.dFmethod);
                %else
                %    [a,b] = numericalSaddlepointWithDerivatives([-z; 0]);
                %    a = a(3);
                %    b = b(3,:);
                %end
                
                r = zeros(d,1);
                for i=1:d
                    r(i) = b(i)/a - rhs(i);
                end
                r=norm(r);
            end
            
            dim = size(omega,1);
            
            f = @(z) mleGoalFun(z, omega);
            Z = fminunc(f, -ones(dim-1,1), optimoptions('fminunc','algorithm','quasi-newton', 'display', 'off'));
            Z = [Z; 0];
        end
    end
end

function P = acgpdf_pcs(X,z,M) %taken from libBingham
    %P = acgpdf_pcs(X,z,;) -- z and M are the sqrt(eigenvalues) and
    %eigenvectors of the covariance matrix; x's are in the rows of X

    S_inv = M*diag(1./(z.^2))*M';

    d = size(X,2);
    P = repmat(1 / (prod(z) * BinghamDistribution.computeUnitSphereSurface(d)), [size(X,1),1]);
    md = sum((X*S_inv).*X, 2);  % mahalanobis distance
    P = P .* md.^(-d/2);
end
