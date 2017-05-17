
classdef Utils
    % This class provides various utility functions.
    %
    % Utils Methods:
    %   getMeanAndCov             - Compute sample mean and sample covariance.
    %   getMeanCovAndCrossCov     - Compute sample mean, covariance, and cross-covariance.
    %   kalmanUpdate              - Perform a Kalman update.
    %   decomposedStateUpdate     - Perform an update for a system state decomposed into two parts A and B.
    %   blockDiag                 - Create a block diagonal matrix.
    %   baseBlockDiag             - Create a block diagonal matrix.
    %   drawGaussianRndSamples    - Draw random samples from a multivariate Gaussian distribution.
    %   resampling                - Perform a simple resampling.
    %   systematicResampling      - Perform a systematic resampling.
    %   rndOrthogonalMatrix       - Creates a random orthogonal matrix of the specified dimension.
    %   getStateSamples           - Get a set of samples approximating a Gaussian distributed system state.
    %   getStateNoiseSamples      - Get a set of samples approximating a jointly Gaussian distributed system state and (system/measurement) noise.
    %   diffQuotientState         - Compute first-order and second-order difference quotients of a function at the given nominal system state.
    %   diffQuotientStateAndNoise - Compute first-order and second-order difference quotients of a function at the given nominal system state and nominal noise.
    
    % >> This function/class is part of the Nonlinear Estimation Toolbox
    %
    %    For more information, see https://bitbucket.org/nonlinearestimation/toolbox
    %
    %    Copyright (C) 2015  Jannik Steinbring <jannik.steinbring@kit.edu>
    %
    %                        Institute for Anthropomatics and Robotics
    %                        Chair for Intelligent Sensor-Actuator-Systems (ISAS)
    %                        Karlsruhe Institute of Technology (KIT), Germany
    %
    %                        http://isas.uka.de
    %
    %    This program is free software: you can redistribute it and/or modify
    %    it under the terms of the GNU General Public License as published by
    %    the Free Software Foundation, either version 3 of the License, or
    %    (at your option) any later version.
    %
    %    This program is distributed in the hope that it will be useful,
    %    but WITHOUT ANY WARRANTY; without even the implied warranty of
    %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %    GNU General Public License for more details.
    %
    %    You should have received a copy of the GNU General Public License
    %    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    methods (Static)
        function [mean, cov] = getMeanAndCov(samples, weights)
            % Compute sample mean and sample covariance.
            %
            % Parameters:
            %   >> samples (Matrix)
            %      Column-wise arranged samples.
            %
            %   >> weights (Row vector)
            %      Column-wise arranged corresponding normalized sample weights.
            %      If no weights are passed, all samples are assumed
            %      to be equally weighted.
            %
            % Returns:
            %   << mean (Column vector)
            %      The sample mean.
            %
            %   << cov (Positive definite matrix)
            %      The sample covariance matrix.
            
            if nargin < 2
                numSamples = size(samples, 2);
                
                % Compute mean
                mean = sum(samples, 2) / numSamples;
                
                % Compute covariance
                if nargout > 1
                    diffSamples = bsxfun(@minus, samples, mean);
                    
                    cov = (diffSamples * diffSamples') / numSamples;
                end
            else
                % Compute mean
                mean = samples * weights';
                
                % Compute covariance
                if nargout > 1
                    diffSamples = bsxfun(@minus, samples, mean);
                    
                    % Weights can be negative => we have to treat them separately
                    
                    % Positive weights
                    idx = weights >= 0;
                    
                    sqrtWeights         = sqrt(weights(idx));
                    weightedDiffSamples = bsxfun(@times, diffSamples(:, idx), sqrtWeights);
                    
                    cov = weightedDiffSamples * weightedDiffSamples';
                    
                    % Negative weights
                    if ~all(idx)
                        idx = ~idx;
                        
                        sqrtWeights         = sqrt(abs(weights(idx)));
                        weightedDiffSamples = bsxfun(@times, diffSamples(:, idx), sqrtWeights);
                        
                        cov = cov - weightedDiffSamples * weightedDiffSamples';
                    end
                end
            end
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov] = getMeanCovAndCrossCov(stateMean, stateSamples, ...
                                                             measSamples, weights)
            % Compute sample mean, covariance, and cross-covariance.
            %
            % Parameters:
            %   >> stateMean (Column vector)
            %      State mean.
            %
            %   >> stateSamples (Matrix)
            %      Column-wise arranged state samples.
            %
            %   >> measSamples (Matrix)
            %      Column-wise arranged measurement samples.
            %
            %   >> weights (Row vector)
            %      Column-wise arranged corresponding normalized sample weights.
            %      If no weights are passed, all samples are assumed
            %      to be equally weighted.
            %
            % Returns:
            %   << measMean (Column vector)
            %      The sample measurement mean.
            %
            %   << measCov (Positive definite matrix)
            %      The sample measurement covariance matrix.
            %
            %   << stateMeasCrossCov (Matrix)
            %      The sample state measurement cross-covariance matrix.
            
            if nargin < 4
                numSamples = size(stateSamples, 2);
                
                % Compute measurement mean
                measMean = sum(measSamples, 2) / numSamples;
                
                % Compute measurement covariance
                diffMeasSamples = bsxfun(@minus, measSamples, measMean);
                
                measCov = (diffMeasSamples * diffMeasSamples') / numSamples;
                
                % Compute state measurement cross-covariance
                diffStateSamples = bsxfun(@minus, stateSamples, stateMean);
                
                stateMeasCrossCov = (diffStateSamples * diffMeasSamples') / numSamples;
            else
                % Compute measurement mean
                measMean = measSamples * weights';
                
                % Compute measurement covariance and state measurement cross-covariance
                diffMeasSamples  = bsxfun(@minus, measSamples, measMean);
                diffStateSamples = bsxfun(@minus, stateSamples, stateMean);
                
                % Weights can be negative => we have to treat them separately
                
                % Positive weights
                idx = weights >= 0;
                
                sqrtWeights              = sqrt(weights(idx));
                weightedDiffMeasSamples  = bsxfun(@times, diffMeasSamples(:, idx), sqrtWeights);
                weightedDiffStateSamples = bsxfun(@times, diffStateSamples(:, idx), sqrtWeights);
                
                measCov = weightedDiffMeasSamples * weightedDiffMeasSamples';
                
                stateMeasCrossCov = weightedDiffStateSamples * weightedDiffMeasSamples';
                
                % Negative weights
                if ~all(idx)
                    idx = ~idx;
                    
                    sqrtWeights              = sqrt(abs(weights(idx)));
                    weightedDiffMeasSamples  = bsxfun(@times, diffMeasSamples(:, idx), sqrtWeights);
                    weightedDiffStateSamples = bsxfun(@times, diffStateSamples(:, idx), sqrtWeights);
                    
                    measCov = measCov - weightedDiffMeasSamples * weightedDiffMeasSamples';
                    
                    stateMeasCrossCov = stateMeasCrossCov - weightedDiffStateSamples * weightedDiffMeasSamples';
                end
            end
        end
        
        function [updatedStateMean, ...
                  updatedStateCov] = kalmanUpdate(stateMean, stateCov, measurement, ...
                                                  measMean, measCov, stateMeasCrossCov)
            % Perform a Kalman update.
            %
            % Parameters:
            %   >> stateMean (Column vector)
            %      Prior state mean.
            %
            %   >> stateCov (Positive definite matrix)
            %      Prior state covariance matrix.
            %
            %   >> measurement (Column vector)
            %      Measurement vector.
            %
            %   >> measMean (Column vector)
            %      Measurement mean.
            %
            %   >> measCov (Positive definite matrix)
            %      Measurement covariance matrix.
            %
            %   >> stateMeasCrossCov (Matrix)
            %      State measurement cross-covariance matrix.
            %
            % Returns:
            %   << updatedStateMean (Column vector)
            %      Posterior state mean.
            %
            %   << updatedStateCov (Positive definite matrix)
            %      Posterior state covariance matrix.
            
            [measCovSqrt, isNonPos] = chol(measCov);
            
            if isNonPos
                error('Utils:InvalidMeasurementCovariance', ...
                      'Measurement covariance matrix is not positive definite.');
            end
            
            % Compute Kalman gain
            A = stateMeasCrossCov / measCovSqrt;
            
            % Compute updated state mean
            updatedStateMean = stateMean + A * (measCovSqrt' \ (measurement - measMean));
            
            % Compute updated state covariance
            updatedStateCov = stateCov - A * A';
        end
        
        function [updatedStateMean, ...
                  updatedStateCov] = decomposedStateUpdate(stateMean, stateCov, stateCovSqrt, ...
                                                           updatedStateMeanA, updatedStateCovA, updatedStateCovASqrt)
            % Perform an update for a system state decomposed into two parts A and B.
            %
            % Parameters:
            %   >> stateMean (Column vector)
            %      Prior mean of the entire state.
            %
            %   >> stateCov (Positive definite matrix)
            %      Prior covariance matrix of the entire state.
            %
            %   >> stateCovSqrt (Square matrix)
            %      Square root of the prior covariance matrix.
            %
            %   >> updatedStateMeanA (Column vector)
            %      Already updated mean of subspace A.
            %
            %   >> updatedStateCovA (Positive definite matrix)
            %      Already updated covariance matrix of subspace A.
            %
            %   >> updatedStateCovASqrt (Square matrix)
            %      Square root of the already updated covariance matrix of subspace A.
            %
            % Returns:
            %   << updatedStateMean (Column vector)
            %      Posterior mean of the entire state.
            %
            %   << updatedStateCov (Positive definite matrix)
            %      Posterior covariance matrix of the entire state.
            
            D = size(updatedStateMeanA, 1);
            
            priorStateMeanA    = stateMean(1:D);
            priorStateMeanB    = stateMean(D+1:end);
            priorStateCovA     = stateCov(1:D, 1:D);
            priorStateCovB     = stateCov(D+1:end, D+1:end);
            priorStateCovBA    = stateCov(D+1:end, 1:D);
            priorStateCovASqrt = stateCovSqrt(1:D, 1:D);
            
            % Computed updated mean, covariance, and cross-covariance for the subspace B
            K = priorStateCovBA / priorStateCovA;
            A = K * updatedStateCovASqrt;
            B = K * priorStateCovASqrt;
            updatedStateMeanB = priorStateMeanB + K * (updatedStateMeanA - priorStateMeanA);
            updatedStateCovB  = priorStateCovB + A * A' - B * B';
            updatedStateCovBA = K * updatedStateCovA;
            
            % Construct updated state mean and covariance
            updatedStateMean = [updatedStateMeanA
                                updatedStateMeanB];
            
            updatedStateCov = [updatedStateCovA  updatedStateCovBA'
                               updatedStateCovBA updatedStateCovB  ];
        end
        
        function blockMat = blockDiag(matrix, numRepetitions)
            % Create a block diagonal matrix.
            %
            % Parameters:
            %   >> matrix (Matrix)
            %      Matrix.
            %
            %   >> numRepetitions (Positive scalar)
            %      Number of matrix repetitions along the diagonal.
            %
            % Returns:
            %   << blockMat (Matrix)
            %      Block diagonal matrix.
            
            blockMat = kron(speye(numRepetitions), matrix);
        end
        
        function blockMat = baseBlockDiag(matrixBase, matrixDiag, numRepetitions)
            % Create a block diagonal matrix.
            %
            % Parameters:
            %   >> matrixBase (Matrix)
            %      Matrix.
            %
            %   >> matrixDiag (Matrix)
            %      Matrix of same size as matrixBase.
            %
            %   >> numRepetitions (Positive scalar)
            %      Number of matrix repetitions along the diagonal.
            %
            % Returns:
            %   << blockMat (Matrix)
            %      Block diagonal matrix.
            
            blockMat = repmat(matrixBase, numRepetitions, numRepetitions) + ...
                       Utils.blockDiag(matrixDiag, numRepetitions);
        end
        
        function rndSamples = drawGaussianRndSamples(mean, covSqrt, numSamples)
            % Draw random samples from a multivariate Gaussian distribution.
            %
            % Parameters:
            %   >> mean (Column vector)
            %      Mean vector.
            %
            %   >> covSqrt (Square matrix)
            %      Square root of the covariance matrix.
            %
            %   >> numSamples (Positive scalar)
            %      Number of samples to draw from the given Gaussian distribution.
            %
            % Returns:
            %   << rndSamples (Matrix)
            %      Column-wise arranged samples drawn from the given Gaussian distribution.
            
            dim = size(mean, 1);
            
            rndSamples = covSqrt * randn(dim, numSamples);
            
            rndSamples = bsxfun(@plus, rndSamples, mean);
        end
        
        function [rndSamples, idx] = resampling(samples, cumWeights, numSamples)
            % Perform a simple resampling.
            %
            % Parameters:
            %   >> samples (Matrix)
            %      Set of column-wise arranged sample positions to resample from.
            %
            %   >> cumWeights (Vector)
            %      Vector containing the cumulative sample weights.
            %
            %   >> numSamples (Positive scalar)
            %      Number of samples to draw from the given sample distribution.
            %
            % Returns:
            %   << rndSamples (Matrix)
            %      Column-wise arranged samples drawn from the given sample distribution.
            %
            %   << idx (Row vector)
            %      Corresponding indices of the samples that were resampled from.
            
            u = rand(1, numSamples);
            
            u = sort(u);
            
            idx = zeros(1, numSamples);
            
            i = 1;
            
            for j = 1:numSamples
                while u(j) > cumWeights(i)
                    i = i + 1;
                end
                
                idx(j) = i;
            end
            
            rndSamples = samples(:, idx);
        end
        
        function [rndSamples, idx] = systematicResampling(samples, cumWeights, numSamples)
            % Perform a systematic resampling.
            %
            % Implements the systematic resampling algorithm from:
            %
            %   Branko Ristic, Sanjeev Arulampalam, and Neil Gordon,
            %   Beyond the Kalman Filter: Particle filters for Tracking Applications,
            %   Artech House Publishers, 2004,
            %   Section 3.3
            %
            % Parameters:
            %   >> samples (Matrix)
            %      Set of column-wise arranged sample positions to resample from.
            %
            %   >> cumWeights (Vector)
            %      Vector containing the cumulative sample weights.
            %
            %   >> numSamples (Positive scalar)
            %      Number of samples to draw from the given sample distribution.
            %
            % Returns:
            %   << rndSamples (Matrix)
            %      Column-wise arranged samples drawn from the given sample distribution.
            %
            %   << idx (Row vector)
            %      Corresponding indices of the samples that were resampled from.
            
            csw = cumWeights * numSamples;
            
            idx = zeros(1, numSamples);
            
            u1 = rand(1);
            
            i = 1;
            
            for j = 1:numSamples
                uj = u1 + (j - 1);
                
                while uj > csw(i)
                    i = i + 1;
                end
                
                idx(j) = i;
            end
            
            rndSamples = samples(:, idx);
        end
        
        function rndMat = rndOrthogonalMatrix(dim)
            % Creates a random orthogonal matrix of the specified dimension.
            %
            % Parameters:
            %   >> dim (Positive scalar)
            %      Dimension of the desired random orthogonal matrix.
            %
            % Returns:
            %   << rndMat (Square matrix)
            %      A random orthogonal matrix of the specified dimension.
            
            mat = randn(dim, dim);
            
            [Q, R] = qr(mat);
            
            D = diag(sign(diag(R)));
            
            rndMat = Q * D;
        end
        
        function [stateSamples, ...
                  weights, ...
                  numSamples] = getStateSamples(sampling, stateMean, stateCovSqrt)
            % Get a set of samples approximating a Gaussian distributed system state.
            %
            % The number of samples, their positions, and their weights are determined (and
            % controlled) by the respective Gaussian sampling technique.
            %
            % Parameters:
            %   >> sampling (Subclass of GaussianSampling)
            %      Gaussian sampling technique that controls the sample generation.
            %
            %   >> stateMean (Column vector)
            %      State mean.
            %
            %   >> stateCovSqrt (Square matrix)
            %      Square root of the state covariance, e.g., the lower matrix of a Cholesky decomposition.
            %
            % Returns:
            %   << stateSamples (Matrix)
            %      Column-wise arranged sample positions approximating the Gaussian system state.
            %
            %   << weights (Row vector)
            %      Column-wise arranged corresponding sample weights.
            %
            %   << numSamples (Positive scalar)
            %      Number of samples approximating the Gaussian system state.
            
            dimState = size(stateMean, 1);
            
            % Get standard normal approximation
            [stdNormalSamples, weights, numSamples] = sampling.getStdNormalSamples(dimState);
            
            % Generate state samples
            stateSamples = stateCovSqrt * stdNormalSamples;
            stateSamples = bsxfun(@plus, stateSamples, stateMean);
        end
        
        function [stateSamples, ...
                  noiseSamples, ...
                  weights, ...
                  numSamples] = getStateNoiseSamples(sampling, stateMean, stateCovSqrt, ...
                                                     noiseMean, noiseCovSqrt)
            % Get a set of samples approximating a jointly Gaussian distributed system state and (system/measurement) noise.
            %
            % It is assumed that state and noise are mutually independent.
            %
            % The number of samples, their positions, and their weights are determined (and
            % controlled) by the respective Gaussian sampling technique.
            %
            % Parameters:
            %   >> sampling (Subclass of GaussianSampling)
            %      Gaussian sampling technique that controls the sample generation.
            %
            %   >> stateMean (Column vector)
            %      State mean.
            %
            %   >> stateCovSqrt (Square matrix)
            %      Square root of the state covariance, e.g., the lower matrix of a Cholesky decomposition.
            %
            %   >> noiseMean (Column vector)
            %      Noise mean.
            %
            %   >> noiseCovSqrt (Square matrix)
            %      Square root of the noise covariance, e.g., the lower matrix of a Cholesky decomposition.
            %
            % Returns:
            %   << stateSamples (Matrix)
            %      Column-wise arranged sample positions approximating the system state.
            %
            %   << noiseSamples (Matrix)
            %      Column-wise arranged samples approximating the noise.
            %
            %   << weights (Row vector)
            %      Column-wise arranged corresponding sample weights.
            %
            %   << numSamples (Positive scalar)
            %      Number of samples approximating the Gaussian joint distribution.
            
            dimState    = size(stateMean, 1);
            dimNoise    = size(noiseMean, 1);
            dimAugState = dimState + dimNoise;
            
            % Get standard normal approximation
            [stdNormalSamples, weights, numSamples] = sampling.getStdNormalSamples(dimAugState);
            
            % Generate state samples
            stateSamples = stateCovSqrt * stdNormalSamples(1:dimState, :);
            stateSamples = bsxfun(@plus, stateSamples, stateMean);
            
            % Generate noise samples
            noiseSamples = noiseCovSqrt * stdNormalSamples(dimState+1:end, :);
            noiseSamples = bsxfun(@plus, noiseSamples, noiseMean);
        end
        
        function [stateJacobian, stateHessians] = diffQuotientState(func, nominalState, step)
            % Compute first-order and second-order difference quotients of a function at the given nominal system state.
            %
            % Parameters:
            %   >> func (Function handle)
            %      System/Measurement function.
            %
            %   >> nominalState (Column vector)
            %      Nominal system state.
            %
            %   >> step (Positive scalar)
            %      Step size for computing the difference quotients.
            %      Default: eps^(1/4)
            %
            % Returns:
            %   << stateJacobian (Square matrix)
            %      First-order difference quotients of the system state
            %      variables, i.e., an approxiamtion of the Jacobian.
            %
            %   << stateHessians (3D matrix)
            %      Set of second-order difference quotients of the system
            %      state variables, i.e., approxiamtions of the Hessians.
            
            % Default value for step
            if nargin < 3
                step = eps^(1/4);
            end
            
            dimState = size(nominalState, 1);
            
            % State Jacobian
            stateSamples = Utils.getJacobianSamples(dimState, nominalState, step);
            
            valuesState = func(stateSamples);
            
            stateJacobian = Utils.getJacobian(dimState, valuesState, step);
            
            if nargout == 2
                % State Hessians
                [stateSamples, L] = Utils.getHessiansSamples(dimState, nominalState, step);
                
                valuesState2 = func(stateSamples);
                
                stateHessians = Utils.getHessians(dimState, valuesState, valuesState2, L, step);
            end
        end
        
        function [stateJacobian, noiseJacobian, ...
                  stateHessians, noiseHessians] = diffQuotientStateAndNoise(func, nominalState, nominalNoise, step)
            % Compute first-order and second-order difference quotients of a function at the given nominal system state and nominal noise.
            %
            % Parameters:
            %   >> func (Function handle)
            %      System/Measurement function.
            %
            %   >> nominalState (Column vector)
            %      Nominal system state.
            %
            %   >> nominalNoise (Column vector)
            %      Nominal noise.
            %
            %   >> step (Positive scalar)
            %      Step size to compute the finite difference.
            %      Default: eps^(1/4)
            %
            % Returns:
            %   << stateJacobian (Square matrix)
            %      First-order difference quotients of the system state
            %      variables, i.e., an approxiamtion of the Jacobian.
            %
            %   << noiseJacobian (Square matrix)
            %      First-order difference quotients of the noise variables,
            %      i.e., an approxiamtion of the Jacobian.
            %
            %   << stateHessians (3D matrix)
            %      Set of second-order difference quotients of the system
            %      state variables, i.e., approxiamtions of the Hessians.
            %
            %   << noiseHessians (3D matrix)
            %      Set of second-order difference quotients of the noise
            %      variables, i.e., approxiamtions of the Hessians.
            
            % Default value for step
            if nargin < 4
                step = eps^(1/4);
            end
            
            dimState = size(nominalState, 1);
            dimNoise = size(nominalNoise, 1);
            
            % State Jacobian
            stateSamples = Utils.getJacobianSamples(dimState, nominalState, step);
            noiseSamples = repmat(nominalNoise, 1, 2 * dimState + 1);
            
            valuesState = func(stateSamples, noiseSamples);
            
            stateJacobian = Utils.getJacobian(dimState, valuesState, step);
            
            % Noise Jacobian
            noiseSamples = Utils.getJacobianSamples(dimNoise, nominalNoise, step);
            stateSamples = repmat(nominalState, 1, 2 * dimNoise + 1);
            
            valuesNoise = func(stateSamples, noiseSamples);
            
            noiseJacobian = Utils.getJacobian(dimNoise, valuesNoise, step);
            
            if nargout == 4
                % State Hessians
                [stateSamples, L] = Utils.getHessiansSamples(dimState, nominalState, step);
                noiseSamples = repmat(nominalNoise, 1, 2 * L);
                
                valuesState2 = func(stateSamples, noiseSamples);
                
                stateHessians = Utils.getHessians(dimState, valuesState, valuesState2, L, step);
                
                % Noise Hessians
                [noiseSamples, L] = Utils.getHessiansSamples(dimNoise, nominalNoise, step);
                stateSamples = repmat(nominalState, 1, 2 * L);
                
                valuesNoise2 = func(stateSamples, noiseSamples);
                
                noiseHessians = Utils.getHessians(dimNoise, valuesNoise, valuesNoise2, L, step);
            end
        end
    end
    
    methods (Static, Access = 'private')
        function samples = getJacobianSamples(dim, nominalVec, step)
            samples = bsxfun(@plus, [step*eye(dim) -step*eye(dim) zeros(dim, 1)], nominalVec);
        end
        
        function jacobian = getJacobian(dim, values, step)
            idx      = 1:dim;
            jacobian = values(:, idx) - values(:, dim + idx);
            jacobian = jacobian / (2 * step);
        end
        
        function [samples, L] = getHessiansSamples(dim, nominalVec, step)
            L     = (dim * (dim + 1)) * 0.5 - dim;
            steps = zeros(dim, L);
            
            a = 1;
            b = dim - 1;
            for i = 1:dim - 1
                d = dim - i;
                
                steps(i,         a:b) = step;
                steps(i + 1:end, a:b) = step * eye(d);
                
                a = b + 1;
                b = a + d - 2;
            end
            
            samples = bsxfun(@plus, [steps -steps], nominalVec);
        end
        
        function hessians = getHessians(dim, values, values2, L, step)
            idx = 1:dim;
            a   = 2 * values(:, end);
            b   = bsxfun(@plus, values2(:, 1:L) + values2(:, L + 1:end), a);
            c   = values(:, idx) + values(:, dim + idx);
            d   = bsxfun(@minus, c, a);
            
            dimFunc  = size(values, 1);
            hessians = nan(dim, dim, dimFunc);
            
            k = 1;
            for i = 1:dim
                hessians(i, i, :) = d(:, i);
                
                for j = (i + 1):dim
                    vec = (b(:, k) - c(:, i) - c(:, j)) * 0.5;
                    
                    hessians(i, j, :) = vec;
                    hessians(j, i, :) = vec;
                    
                    k = k + 1;
                end
            end
            
            hessians = hessians / (step * step);
        end
    end
end
