classdef SE2BinghamDistribution < AbstractSE2Distribution
    % Distribution on S^{1} x R^2.
    %   We assume the following parametrization of the density
    %     f(x) = exp(x^T * C * x).
    %
    % Igor Gilitschenski, Gerhard Kurz, Simon J. Julier, Uwe D. Hanebeck,
    % A New Probability Distribution for Simultaneous Representation of Uncertain Position and Orientation
    % Proceedings of the 17th International Conference on Information Fusion (Fusion 2014),
    % Salamanca, Spain, July 2014.
    
    properties
        % this is a bit redundant, we could use functions to obtain C1,C2,C3 from C or vice versa
        C (4,4) double  % Parameter matrix
        C1 (2,2) double % Parameter matrix for spherical part.
        C2 (2,2) double % Parameter matrix for dependence.
        C3 (2,2) double % Parameter matrix for unrestricted part.
        NC (1,1) double {mustBePositive} = 1 % Normalization constant
    end
    
    methods
        function this = SE2BinghamDistribution(C, C2, C3)
            % Generate an SE2BinghamDistribution object.
            %
            % Parameters:
            %   C - Full Parameter Matrix (Alternative: sperical part).
            %   C2 - Parameter matrix for dependency (optional).
            %   C3 - Parameter matrix for unrestricted part (optional).
            arguments
                C (:,:) double
                C2 (:,:) double = []
                C3 (:,:) double = []
            end
            this.linD = 2;
            this.boundD = 1;
            this.dim = 3;
            
            assert(isempty(C2)==isempty(C3), 'Wrong number of arguments');
            
            if isempty(C2)
                % Sets main parameter matrix C and adapts other matrices.
                assert(isequal(size(C),[4 4]), 'C has wrong size.');
                assert(isequal(C,C'), 'C must be symmetric.');
                
                this.C = C;
                this.C1 = C(1:2,1:2);
                this.C2 = C(1:2,3:4)';
                this.C3 = C(3:4,3:4);
            else
                % Check sizes.
                assert(isequal(size(C),[2 2]), 'C has wrong size.');
                assert(isequal(size(C2),[2 2]), 'C2 has wrong size.');
                assert(isequal(size(C3),[2 2]), 'C3 has wrong size.');
                
                % Check symmetry
                assert(isequal(C,C'), 'C must be symmetric.');
                assert(isequal(C3,C3'), 'C3 must be symmetric.');
                
                this.C = [C C2'; C2 C3];
                this.C1 = C;
                this.C2 = C2;
                this.C3 = C3;
            end
            
            this.NC = this.computeNC();
            
            assert(all(eig(this.C3)<=0), ...
                'Parameter matrix for unrestricted part must be negative definite.');
        end
        
        function p = pdf(this, xa)
            arguments
                this (1,1) SE2BinghamDistribution
                xa (:,:) double
            end
            if size(xa)==3 % Convert to quaternion if given in angle-pos form
                xa = AbstractSE2Distribution.anglePosToDualQuaternion(xa);
            end    
            assert(size(xa,1)==4);
            p = 1/this.NC * exp(sum(xa.*(this.C*xa)));
        end
        
        function C = computeCovarianceMCMC(this, n)
            % Computes covariance using n random samples
            %
            % Parameters:
            %   n (scalar)
            %       number of samples (optional)
            % Returns:
            %   C (4x4 matrix)
            %       covariance of the distribution
            
            % Set default number of samples
            arguments
                this (1,1) SE2BinghamDistribution
                n (1,1) double {mustBeInteger,mustBePositive} = 10000
            end            
            samples = this.sample(n);
            C = cov(samples',1);
        end
        
        function nc = computeNC(this)
            % Compute normalization constant of the distribution.
            % 
            % Returns:
            %   nc (scalar)
            %       normalization constant
            arguments
                this (1,1) SE2BinghamDistribution
            end  
            BM = this.C1-(this.C2' * pinv(this.C3) * this.C2);
            [~,Z] = eig(BM);
            Z = sort(diag(Z),'ascend');
            bNC = BinghamDistribution.computeF(Z);
            
            nc = 2*pi*sqrt(det(-0.5 * pinv(this.C3)))*bNC;
        end
        
        function m = mode(this)
            % Computes one of the modes of the distribution.
            %
            % Because the distribution has antipodal symmetry -m is also a
            % valid mode.
            %
            % Returns:
            %   m (4 x 1 column vector)
            %       mode of the distribution
            arguments
                this (1,1) SE2BinghamDistribution
            end 
            m = zeros(4,1);
            
            % Compute the rotational part
            binghamC = this.C1 - this.C2'*pinv(this.C3)*this.C2;
            [M, Z] = eig(binghamC);
            
            if Z(1,1) <= Z(2,2)
                m(1:2) = M(:,2);
            else
                m(1:2) = M(:,1);
            end
            
            m(3:4) = - pinv(this.C3)*this.C2*m(1:2);
        end
        
        
        function [samples, weights] = sampleDeterministic(this)
            % Generates deterministic samples.
            % Uses default sampling of the Bingham class and a naive variant
            % of the UKF for deterministic sampling.
            %
            % Returns:
            %   samples (4 x n matrix)
            %       deterministic samples, one per column
            %   weights (1 x n row vector)
            %       weight for each sample
            %
            % Igor Gilitschenski, Gerhard Kurz, Uwe D. Hanebeck,
            % A Stochastic Filter for Planar Rigid-Body Motions
            % Proceedings of the 2015 IEEE International Conference on 
            % Multisensor Fusion and Information Integration (MFI 2015), 
            % San Diego, California, USA, September 2015.
            
            % Generate Sampling of Bingham marginal
            arguments
                this (1,1) SE2BinghamDistribution
            end
            b = this.getBinghamMarginal();
            [bSamples, bWeights] = b.sampleDeterministic('uniform');
            numBinghamSamples = length(bWeights);
            
            % Generate Samples using UKF.
            % Generate Gaussian Samples using UKF.
            ukfSampling = GaussianSamplingUKF();
            [gSamples, gWeights, numGaussianSamples] = ukfSampling.getSamples(Gaussian(zeros(2,1), -0.5*pinv(this.C3)));
            if isscalar(gWeights)
                gWeights = repmat(gWeights, 1, numGaussianSamples);
            end            
            
            % Compute repositioning matrix
            A = -pinv(this.C3)*this.C2;
            
            % Initialize return paramters.
            numSamples = numBinghamSamples * numGaussianSamples;
            samples = zeros(4,numSamples);
            weights = zeros(1,numSamples);
            
            % Computation of the actual samples
            for i=1:numBinghamSamples
                for j = 1:numGaussianSamples
                    cur = numGaussianSamples*(i-1) + j;
                    samples(1:2,cur) = bSamples(:,i);
                    samples(3:4,cur) = gSamples(:,j) + A*bSamples(:,i);
                    weights(cur) = bWeights(i)*gWeights(j);
                end
            end
        end
        
        function s = sample(this,n)
            % Samples from current distribution.
            %
            % A two step sampling procedure is used. It is based on sampling a
            % Bingham density first, because it is the marginal density of the
            % first two entries. Then, this is used to sample a corresponding
            % Gaussian, which has a mean dependent on the Bingham.
            %
            %   Parameters:
            %       n - number of samples.
            %
            %   Returns:
            %       s - (4 x n) Matrix containing the s samples.
            arguments
                this (1,1) SE2BinghamDistribution
                n (1,1) double {mustBeInteger,mustBePositive} = 10000
            end
            s = zeros(4,n); % Initialize return values.
            
            % STEP 1: Sample from Bingham density
            % Compute Bingham parameters using the Eigendecomposition of
            % the Schur complement of C.
            binghamC = this.C1 - this.C2'*pinv(this.C3)*this.C2;
            [M, Z] = eig(binghamC);
            
            if Z(1,1)>Z(2,2)
                M = [M(:, 2) M(:,1)];
                Z = diag([Z(2,2) Z(1,1)]);
            end
            
            %             assert(Z(1,1)<=Z(2,2), ...
            %                 'Eigenvalues seem to be in wrong order. Implement sort.');
            z=diag(Z)-Z(2,2);
            
            % Generate the actual samples.
            b = BinghamDistribution(z,M);
            binghamSamples = b.sample(n);
            s(1:2,:) = binghamSamples;
            
            % STEP 2: Sample corresponding Gaussian
            % The mean of the gaussian given the circular part x_c
            % is -C_3^{-1} * C_2 * x_c
            % The Covariance is given by -0.5 * C_3^{-1}
            
            R = chol(-0.5*pinv(this.C3));
            means = (-pinv(this.C3)*this.C2*binghamSamples);
            s(3:4,:) =  means + (randn(n,2)*R)';
        end
        
        function h = plotState(this, scalingFactor, circleColor, angleColor, samplesForMatching)
            arguments
                this (1,1) AbstractSE2Distribution
                scalingFactor (1,1) double = 1
                circleColor (3,1) double = [0    0.4470    0.7410]
                angleColor (3,1) double = [0.8500    0.3250    0.0980]
                samplesForMatching (1,1) double {mustBeInteger,mustBePositive} = 10000
            end
            samples = this.sample(samplesForMatching);
            dist = SE2PWDDistribution(AbstractSE2Distribution.dualQuaternionToAnglePos(samples));
            h = dist.plotState(scalingFactor, circleColor, angleColor);
        end
    end
    
    methods (Access=private)
        function b = getBinghamMarginal(this)
            % Computes Bingham marginal of circular part.
            %   Caution: The normalization constant of the returned Bingham
            %   distribution object does not correspond to the true
            %   normalization constant of the rotational part.
            %
            %   Returns:
            %       b - BinghamDistribution object.
            
            BM = this.C1-(this.C2' * pinv(this.C3) * this.C2);
            [M,Z] = eig(BM);
            [Z,order] = sort(diag(Z),'ascend');
            Z=Z-Z(2);
            M = M(:,order);
            b = BinghamDistribution(Z,M);
        end
    end
    
    methods (Static)
        function res = fit(samples, weights)
            % Estimates parameters of SE2 distribution from samples.
            %
            %   Parameters:
            %       samples - (4 x n) Matrix with a sample in each column.
            %
            %   Return Value:
            %       res - resulting distribution object.
            arguments
                samples (:,:) double
                weights (1,:) = ones(1,size(samples,2))/size(samples,2)
            end
            if size(samples,1)==3
                % Given in angle-pos representation, convert to quaternions
                samples = AbstractSE2Distribution.anglePosToDualQuaternion(samples);
            end
            assert(size(samples,1)==4);
            assert(size(weights,2) == size(samples,2));
            
            % Estimate Bingham parameters first
            b = BinghamDistribution.fit(samples(1:2,:), weights);
            
            % This corresponds to C1 - C2'*C3^{-1}*C2.
            Tmp = b.M*diag(b.Z)*b.M';
            
            % Estimation of -0.5*C3^{-1} and -C3^{-1}*C2 using formulas from
            % (Anderson, 2003, p. 294 Th. 8.2.1)
            
            % todo check use of weights           
            regC = samples(3:4,:).*weights*samples(1:2,:)'*size(samples,2);
            regA = samples(1:2,:).*weights*samples(1:2,:)'*size(samples,2);
            regBeta = regC*pinv(regA);
            
            regCov = (samples(3:4,:)-regBeta*samples(1:2,:)) .* ...
                weights * (samples(3:4,:)-regBeta*samples(1:2,:))';
            
            % Transform into distribution parameters.
            C3est = pinv(-2*regCov);
            
            % Numerics break symmetry. We fix it.
            C3est = 0.5*(C3est + C3est');
            C2est = -C3est * regBeta;
            
            % Obtain estimate of C1.
            C1est = Tmp + C2est'*pinv(C3est)*C2est;
            
            % Same problem here.
            C1est = 0.5*(C1est + C1est');
            
            res = SE2BinghamDistribution(C1est, C2est, C3est);
        end
    end
end

