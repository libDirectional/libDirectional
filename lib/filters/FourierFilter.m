classdef FourierFilter < AbstractCircularFilter
    % Fourier filters.  Depending on the parameter in the initialization,
    % either Fourier identity or Fourier square root filter is used.
    %
    % see Florian Pfaff, Gerhard Kurz, and Uwe D. Hanebeck,
    % Multimodal Circular Filtering Using Fourier Series
    % Proceedings of the 18th International Conference on Information Fusion (Fusion 2015),
    % Washington, D.C., USA, July 2015.
    %
    % Florian Pfaff, Gerhard Kurz, Uwe D. Hanebeck,
    % Nonlinear Prediction for Circular Filtering Using Fourier Series
    % Proceedings of the 19th International Conference on Information Fusion (Fusion 2016),
    % Heidelberg, Germany, July 2016.
    
    properties
        fd FourierDistribution
    end
    
    methods
        function this = FourierFilter(noOfCoefficients, transformation)
            % Constructor
            %
            % Parameters:
            %   noOfCoefficients (odd integer > 0)
            %       number of Fourier coefficients to use
            arguments
                noOfCoefficients (1,1) double
                transformation char {mustBeMember(transformation,{'sqrt','identity'})} = 'sqrt'
            end
            this.fd = FourierDistribution.fromDistribution(CircularUniformDistribution(), noOfCoefficients, transformation);
        end
        
        function setState(this, fd_)
            % Sets the current system state
            %
            % Parameters:
            %   fd_ (FourierDistribution)
            %       new state
            arguments
                this (1,1) FourierFilter
                fd_ (1,1) AbstractCircularDistribution
            end
            if ~(isa(fd_, 'FourierDistribution'))
                warning('setState:nonFourier', 'fd_ is not a FourierDistribution. Transforming with a number of coefficients that is equal to that of the filter.');
                fd_ = FourierDistribution.fromDistribution(fd_, 2*length(this.fd.a)-1, this.fd.transformation);
            else
                if ~strcmp(this.fd.transformation, fd_.transformation)
                    warning('setState:transDiffer', 'New density is transformed differently.');
                end
                if length(fd_.a) ~= length(this.fd.a)
                    warning('setState:noOfCoeffsDiffer', 'New density has different number of coefficients.')
                end
            end
            this.fd = fd_;
        end
        
        function fd = getEstimate(this)
            % Return current estimate
            %
            % Returns:
            %   fd (FourierDistribution)
            %       current estimate
            arguments
                this (1,1) FourierFilter
            end
            fd = this.fd;
        end
        
        function mean = getEstimateMean(this)
            arguments
                this (1,1) FourierFilter
            end
            mean = this.fd.circularMean;
        end
        
        function predictIdentity(this, dSys)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by dSys.
            %
            % Parameters:
            %   dSys (AbstractCircularDistribution or samples on grid)
            %       distribution of additive noise or samples on a grid. If
            %       values on a grid are given, take care that they are
            %       not rooted yet. Take care that if
            %       linspace(0,2*pi,noCoeffs) is used, that the last grid
            %       point has to be discarded. For example
            %       xvals=0:2*pi/noCoeffs:2*pi-2*pi/noCoeffs;
            if isa(dSys, 'AbstractCircularDistribution')
                % If a distribution on a circle is given, use Hadamard
                % product approach
                if ~isa(dSys, 'FourierDistribution')
                    warning('Predict:automaticConversion', ...
                        'dSys is not a FourierDistribution. Transforming with a number of coefficients that is equal to that of the filter. For non-varying noises, transforming once is much more efficient and should be preferred.');
                    dSys = FourierDistribution.fromDistribution(dSys, 2*length(this.fd.a)-1, this.fd.transformation);
                end
                this.fd = this.fd.convolve(dSys);
            elseif isnumeric(dSys)
                % Use cconv-based approach
                noCoeffs = 2 * length(this.fd.a) - 1;
                assert(strcmp(this.fd.transformation, 'sqrt'), 'Only sqrt transformation currently supported');
                assert(numel(dSys) == noCoeffs, 'Assume that as many grid points are used as there are coefficients.');
                fdvals = (ifft(ifftshift(this.fd.c), 'symmetric')).^2;
                fPredVals = cconv(fdvals, dSys, noCoeffs) * noCoeffs * 2 * pi;
                this.fd = FourierDistribution.fromComplex(fftshift(fft(sqrt(fPredVals)))/noCoeffs, 'sqrt');
            else
                error('Input format of dSys is not supported');
            end
            
        end
        
        function updateIdentity(this, dMeas, z)
            % Updates assuming identity measurement model, i.e.,
            % z(k) = x(k) + v(k)    mod 2pi,
            % where v(k) is additive noise given by dMeas.
            %
            % Parameters:
            %   dMeas (AbstractCircularDistribution)
            %       distribution of additive noise
            %   z (scalar)
            %       measurement in [0, 2pi)
            arguments
                this (1,1) FourierFilter
                dMeas (1,1) AbstractCircularDistribution
                z (1,1) double = 0
            end
            if ~isa(dMeas, 'FourierDistribution')
                warning('Update:automaticConversion', ...
                    'dMeas is not a FourierDistribution. Transforming with a number of coefficients that is equal to that of the filter. For non-varying noises, transforming once is much more efficient and should be preferred.');
                dMeas = FourierDistribution.fromDistribution(dMeas, 2*length(this.fd.a)-1, this.fd.transformation);
            end
            if z~=0
                dMeasShifted = dMeas.shift(z);
            else
                dMeasShifted = dMeas;
            end
            this.fd = this.fd.multiply(dMeasShifted, 2*length(this.fd.a)-1);
        end
        
        function predictNonlinear(this, f, noiseDistribution, truncateJointSqrt)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k)    mod 2pi,
            % where w(k) is additive noise given by noiseDistribution.
            % Using predictNonlinearViaTransitionDensity with a
            % pretransformed fTrans is to be preferred as the function is
            % always performed anew when using predictNonlinear.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi) to [0,2pi).
            %       fuction MUST support row vectors of different angles
            %       and output a row vector!
            %   noiseDistribution (AbstractCircularDistribution)
            %       distribution of additive noise
            arguments
                this (1,1) FourierFilter
                f (1,1) function_handle
                noiseDistribution (1,1) AbstractCircularDistribution
                truncateJointSqrt (1,1) logical = true
            end
            fTrans = @(xkk, xk)reshape(noiseDistribution.pdf(xkk(:)'-f(xk(:)')), size(xk));
            this.predictNonlinearViaTransitionDensity(fTrans, truncateJointSqrt);
            
        end
        
        function predictNonlinearViaTransitionDensity(this, fTrans, truncateJointSqrt)
            % Predicts assuming a nonlinear system model using a
            % probabilistic model. See also HypertoroidalFourierFilter that
            % implements this for arbitary dimensions in a complex Fourier
            % series representation.
            %
            % Parameters:
            %   fTrans
            %       transition density f(x(k+1)|x(k))
            %       Either given as ToroidalFourierDistribution or as
            %       function.
            %       If fTrans is a ToroidalFourierDistribution:
            %       First dimension is x(k+1), second is x(k). If given as
            %       a function, the first input argument is for x(k+1) and
            %       the second for x(k).
            %
            % To reuse code, fTrans is treated as a distribution and is
            % thus normalized to one although it should integrate to
            % 2*pi. This is because fTrans is not a bivariate density but a
            % conditional density that, when integrated over x_k+1
            % always yields 1. By also integrating over x_k, we get
            % 2*pi.
            arguments
                this (1,1) FourierFilter
                fTrans
                truncateJointSqrt (1,1) logical = true
            end
            
            % As explained above, fTrans acutally integrates to 2*pi.
            % This would result in a warning due to the lack of
            % normalization. Warning is therefore disabled.
            noOfCoeffs = 2 * length(this.fd.a) - 1;
            warnStruct = warning('off', 'Normalization:notNormalized');
            if isa(fTrans, 'function_handle')
                assert(nargin(fTrans) == 2);
                fTrans = ToroidalFourierDistribution.fromFunction(fTrans, noOfCoeffs*[1, 1], this.fd.transformation);
            else
                assert(strcmp(this.fd.transformation, fTrans.transformation));
            end
            
            if strcmp(this.fd.transformation, 'identity')
                % Multiply and marginalize in one step. See test case in
                % FourierFilter test that illustrates the validity of this
                % approach. On top of the 2*pi from the marginalization,
                % we also have to respect that we have fTrans/(2*pi).
                % Thus, we multiply by (2*pi)^2.
                % Note: fd.c does NOT need to be transposed because x_k is
                % already a row vector by default.
                cPredictedId = (2 * pi)^2 * conv2(fTrans.C, this.fd.c, 'valid');
            elseif strcmp(this.fd.transformation, 'sqrt')
                % There are two ways to ensure nonnegativity of the result.
                % 1) Calculating the full coefficient vectors/tensors for
                % the identity transformation for both ('full' is needed as
                % coefficint vectors must not be truncated) and then by
                % convolving them with 'valid' is sufficient to obtain the
                % result of the marginalization. This can be implemented as
                % conv2(conv2(fTrans.C,fTrans.C,'full'),conv2(this.fd.c,this.fd.c,'full'),'valid')
                % 2) Calculating the (possibly truncated) result for the
                % coefficient tensor for the square root of the joint and
                % then convolve the result with itself. We use this
                % approach as it is more flexible and featured a better run
                % time performance (even if we chose to not truncate the
                % joint density).
                if ~truncateJointSqrt
                    % Since we are still in square root representation, we can
                    % scale fTrans with sqrt(2*pi) to ensure normalization
                    % Note: fd.c does NOT need to be transposed because it x_k is
                    % already a row vector by default.
                    cJointSqrt = conv2(sqrt(2*pi)*fTrans.C, this.fd.c, 'full');
                else
                    % See commments above. Due to the additional
                    % truncation, normalization is not ensured even
                    % though we scale it by sqrt(2*pi)
                    cJointSqrt = conv2(sqrt(2*pi)*fTrans.C, this.fd.c, 'same');
                end
                % Do NOT simply convolve with 'same' to prevent truncation
                % along the preserved dimension. We pad along this
                % dimension to preserve the relevant entries. Any truncation along
                % this dimension in the identity representation could
                % otherwise lead to negative function values. The entries
                % along the other dimension are to be discarded as the
                % dimension is marginalized out. An additional factor
                % 2*pi is required for preserving the normalization in the
                % the marginalization.
                % It is faster to use this padding and then call conv2 with
                % 'valid' than to calculate 'full' and then truncate.
                additionalColumns = 2 * length(this.fd.b); % If we set this value to length(this.fd.b) instead of 2*length(this.fd.b), we would only get the desired number of coefficients but we could no longer guarantee that the function is always nonnegative.
                cPredictedId = 2 * pi * conv2( ...
                    [zeros(additionalColumns, size(cJointSqrt, 2)); cJointSqrt; zeros(additionalColumns, size(cJointSqrt, 2))], ...
                    cJointSqrt, 'valid');
            end
            if strcmp(this.fd.transformation, 'identity') || ~truncateJointSqrt
                % If identity transformation is used or no truncation has
                % been performed, no normalization issues are expected
                % (except numerical issues)
                warning(warnStruct); % Reenable warnings first as should be no problem with identity
                % No need to transpose becaues fromComplex expects rowVector (unlike HypertoroidalFourierDistribution in the 1-D case)
                this.fd = FourierDistribution.fromComplex(cPredictedId.', 'identity');
            else % Normalization is not ensured. Disable warning afterward
                this.fd = FourierDistribution.fromComplex(cPredictedId.', 'identity');
                warning(warnStruct);
            end
            
            if strcmp(fTrans.transformation, 'sqrt')
                this.fd = this.fd.transformViaFFT('sqrt', noOfCoeffs);
            end
        end
        
        function updateNonlinear(this, likelihood, z) %measurement z, likelihood(z,x)=P(Z|X)
            arguments
                this (1,1) FourierFilter
                likelihood (1,1) function_handle
                z (:,1) double
            end
            fdMeas = FourierDistribution.fromFunction( ...
                ... % Assume likelihood can use implicit expansion (for scalars also possible in older Matlab versions)
                @(x)reshape(likelihood(z, x(:)'), size(x)), ...
                2*numel(this.fd.a)-1, this.fd.transformation);
            this.updateIdentity(fdMeas, zeros(size(z)));
        end
        
        function updateNonlinearViaIFFT(this, likelihood, z) %measurement z, likelihood(z,x)=P(Z|X)
            % Transform to state space via IFFT, then multiply with
            % likelihood and transform back. This leads to (more) aliasing
            % than the other update approach if no padding is used. In
            % FourierFilterTest, we show that we can obtain the same result
            % as using the other update approach when we pad sufficiently.
            arguments
                this (1,1) FourierFilter
                likelihood (1,1) function_handle
                z (:,1) double
            end
            cCurr = this.fd.c;
            priorVals = ifft(ifftshift(cCurr), 'symmetric') * numel(cCurr);
            xVals = linspace(0, 2*pi, numel(cCurr)+1);
            if strcmp(this.fd.transformation, 'identity')
                posteriorVals = priorVals .* likelihood(z, xVals(1:end-1));
            elseif strcmp(this.fd.transformation, 'sqrt')
                % Nonnegativity of pdf is not voided because we stay
                % in sqrt space.
                posteriorVals = priorVals .* sqrt(likelihood(z, xVals(1:end-1)));
            else
                error('Transforamtion currently not supported');
            end
            warnStruct = warning('off', 'Normalization:notNormalized');
            this.fd = FourierDistribution.fromComplex(fftshift(fft(posteriorVals)), this.fd.transformation);
            warning(warnStruct);
        end
        
        function likelihoodVal = associationLikelihood(this, likelihood)
            % see Florian Pfaff, Kailai Li, and Uwe D. Hanebeck,
            % Association Likelihoods for Directional Estimation
            % Proceedings of the 2019 IEEE International Conference on
            % Multisensor Fusion and Integration for Intelligent Systems (MFI 2019),
            % Taipei, Republic of China, May, 2019.
            arguments
                this (1,1) FourierFilter
                likelihood (1,1) FourierDistribution
            end
            assert(numel(this.getEstimate.a) == numel(likelihood.a));
            assert(numel(this.getEstimate.transformation) == numel(likelihood.transformation));
            
            if strcmp(this.getEstimate.transformation, 'identity')
                likelihoodVal = 2 * pi * real(this.getEstimate.c*likelihood.c');
            elseif strcmp(this.getEstimate.transformation, 'sqrt')
                likelihoodVal = 2 * pi * norm(conv(this.getEstimate.c, likelihood.c)).^2;
            else
                error('Transformation not supported');
            end
        end
        
    end
    
end
