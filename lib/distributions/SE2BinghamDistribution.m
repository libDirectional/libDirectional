classdef SE2BinghamDistribution
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
        C  % Parameter matrix
        C1 % Parameter matrix for spherical part.
        C2 % Parameter matrix for dependence.
        C3 % Parameter matrix for unrestricted part.
        NC % Normalization constant
    end
    
    methods
        function this = SE2BinghamDistribution(C, C2, C3)
            % Generate an SE2BinghamDistribution object.
            %
            % Parameters:
            %   C - Full Parameter Matrix (Alternative: sperical part).
            %   C2 - Parameter matrix for dependency (optional).
            %   C3 - Parameter matrix for unrestricted part (optional).
            
            assert(nargin == 3 || nargin == 1, 'Wrong number of Arguments');
            
            if nargin == 1
                % Sets main parameter matrix C and adapts other matrices.
                assert(all(size(C)==[4 4]), 'C has wrong size.');
                assert(all(all(C==C')), 'C must be symmetric.');
                
                this.C = C;
                this.C1 = C(1:2,1:2);
                this.C2 = C(1:2,3:4)';
                this.C3 = C(3:4,3:4);
            else % nargin == 3
                % Check sizes.
                assert(all(size(C)==[2 2]), 'C has wrong size.');
                assert(all(size(C2)==[2 2]), 'C2 has wrong size.');
                assert(all(size(C3)==[2 2]), 'C3 has wrong size.');
                
                % Check symmetry
                assert(all(all(C==C')), 'C must be symmetric.');
                assert(all(all(C3==C3')), 'C3 must be symmetric.');
                
                this.C = [C C2'; C2 C3];
                this.C1 = C;
                this.C2 = C2;
                this.C3 = C3;
                this.NC = this.computeNC();
            end
            
            assert(all(eig(this.C3)<=0), ...
                'Parameter matrix for unrestricted part must be negative definite.');
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
            if nargin < 2
                n = 10000;
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
            b = this.getBinghamMarginal();
            [bSamples, bWeights] = b.sampleDeterministic('uniform');
            numBinghamSamples = length(bWeights);
            
            % Generate Samples using UKF.
            % Generate Gaussian Samples using UKF.
            ukfSampling = GaussianSamplingUKF();
            [gSamples, gWeights, numGaussianSamples] = ukfSampling.getSamples(Gaussian(zeros(2,1), -0.5*pinv(this.C3)));
            
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
            
            assert(isscalar(n));
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
        function res = fit(samples)
            % Estimates parameters of SE2 distribution from samples.
            %
            %   Parameters:
            %       samples - (4 x n) Matrix with a sample in each column.
            %
            %   Return Value:
            %       res - resulting distribution object.
            
            sampleSize = size(samples);
            assert(sampleSize(1) == 4, 'Wrong sample matrix dimension');
            
            % Estimate Bingham parameters first
            b = BinghamDistribution.fit(samples(1:2,:));
            
            % This corresponds to C1 - C2'*C3^{-1}*C2.
            Tmp = b.M*diag(b.Z)*b.M';
            
            % Estimation of -0.5*C3^{-1} and -C3^{-1}*C2 using formulas from
            % (Anderson, 2003, p. 294 Th. 8.2.1)
            
            regC = samples(3:4,:)*samples(1:2,:)';
            regA = samples(1:2,:)*samples(1:2,:)';
            regBeta = regC*pinv(regA);
            
            numSamples = size(samples,2);
            regCov = (samples(3:4,:)-regBeta*samples(1:2,:)) * ...
                (samples(3:4,:)-regBeta*samples(1:2,:))'/numSamples;
            
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

