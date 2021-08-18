classdef SphericalHarmonicsFilter < AbstractHypersphericalFilter
    % Filter using spherical harmonics
    
    % see Florian Pfaff, Gerhard Kurz, and Uwe D. Hanebeck,
    % Filtering on the Unit Sphere Using Spherical Harmonics
    % Proceedings of the 2017 IEEE International Conference 
    % on Multisensor Fusion and Integration for Intelligent Systems (MFI 2017),
    % Daegu, Korea, November 2017.
    properties
        state SphericalHarmonicsDistributionComplex
    end
    
    methods
        function this = SphericalHarmonicsFilter(degree, transformation)
            arguments
                degree (1,1) double
                transformation char = 'identity'
            end
            coeffMat = zeros(degree+1, 2*degree+1);
            switch transformation
                case 'identity'
                    coeffMat(1, 1) = 1 / sqrt(4*pi);
                case 'sqrt'
                    coeffMat(1, 1) = 1;
            end
            this.state = SphericalHarmonicsDistributionComplex(coeffMat, transformation);
        end
        
        function setState(this, state)
            assert(isa(state, 'SphericalHarmonicsDistributionComplex'));
            if ~strcmp(this.state.transformation, state.transformation)
                warning('setState:transDiffer', 'New density is transformed differently.');
            end
            if ~isequal(size(this.state.coeffMat), size(state.coeffMat))
                warning('setState:noOfCoeffsDiffer', 'New density has different number of coefficients.')
            end
            this.state = state;
        end
        
        function est = getEstimate(this)
            est = this.state;
        end
        
        function mean = getEstimateMean(this)
            mean = this.state.meanDirection;
        end
        
        function predictIdentity(this, sysNoise)
            % sysNoise has to be zonal.
            assert(isa(sysNoise, 'SphericalHarmonicsDistributionComplex'));
            if strcmp(this.state.transformation, 'sqrt') && strcmp(sysNoise.transformation, 'identity')
                assert(2*(size(this.state.coeffMat, 1) - 1) == size(sysNoise.coeffMat, 1)-1, ...
                    'If sqrt variant is used and sysNoise is given in identity form, sysNoise should normally have a higher degree.');
            end
            this.setState(this.state.convolve(sysNoise));
        end
        
        function updateIdentity(this, measNoise, z)
            % Splitting the noise into a noise term and a measurement is not a
            % universal solution. This function was added for interface
            % compatibility with the VMF Filter. If z is used, measNoise needs to be
            % zonal with rotational symmetry around the z-axis and have a
            % mean of [0;0;1].
            % It is possible to use 0 or [0;0;1] as measurement if no rotation should be
            % performed. Please note that it is usually more efficient to encode it in
            % measNoise than to perform a rotation.
            
            assert(isa(measNoise, 'SphericalHarmonicsDistributionComplex'));
            
            if (norm(z) > 1E-6) && norm(z-[0; 0; 1]) > 1E-6
                warning('SphericalHarmonicsFilter:rotationRequired', ...
                    'Performance may be low for z~=[0;0;1]. Using updateNonlinear may yield faster results.');
                [phi, theta] = cart2sph(z(1), z(2), z(3));
                theta = -theta + pi / 2;
                measNoise = measNoise.rotate(0, theta, phi);
            end
            this.setState(this.state.multiply(measNoise));
        end
        
        function updateNonlinear(this, likelihood, z) %measurement z, likelihood(z,x)=P(Z|X)
            % No zonality is required. Can use array valued likelihoods and
            % matrix using cell array of likelihoods and measurements. The
            % update step used here is similar to the updateStepUsingIFFT
            % in the FourierDistribution class. The convolution-like update
            % step is not used because it is not very efficient.
            degree = size(this.state.coeffMat, 1) - 1;
            lat = linspace(0, 2*pi, 2*degree+2);
            lon = linspace(pi/2, -pi/2, degree+2);
            [latMesh, lonMesh] = meshgrid(lat, lon);
            [x1, x2, x3] = sph2cart(latMesh, lonMesh, 1);
            fvalCurr = plm2xyz(this.state.plmWithDegreeAndOrder);
            
            % Likelihood is expected to evaluate 3 x samplePoints matrices
            % (just as .pdf does)
            if ~iscell(likelihood)
                likelihoodVals = reshape(likelihood(z, [x1(:)'; x2(:)'; x3(:)']), size(latMesh));
            else % If likelihood is a cell array, multiply likelihoods for all
                likelihoodVals = likelihood{1}(z{1}, [x1(:)'; x2(:)'; x3(:)']);
                for i = 2:numel(likelihood)
                    likelihoodVals = likelihoodVals .* likelihood{i}(z{i}, [x1(:)'; x2(:)'; x3(:)']);
                end
                likelihoodVals = reshape(likelihoodVals, size(latMesh));
            end
            
            switch this.state.transformation
                % Also multiply with 2^numel(likelihood) to avoid that an
                % error will occur because the density is so far from being
                % normalized
                case 'identity'
                    fvalNew = 2^numel(likelihood) * fvalCurr .* likelihoodVals;
                case 'sqrt'
                    fvalNew = 2^numel(likelihood) * fvalCurr .* sqrt(likelihoodVals);
                otherwise
                    error('Transformation unsupported');
            end
            
            warningSettings = warning('off', 'Normalization:notNormalized');
            this.state = SphericalHarmonicsDistributionComplex.fromPlmWithDegreeAndOrder(xyz2plm(fvalNew), this.state.transformation);
            warning(warningSettings);
        end
    end
    
end