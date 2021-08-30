classdef (Abstract) AbstractSphericalHarmonicsDistribution < AbstractHypersphericalDistribution
    
    % see Florian Pfaff, Gerhard Kurz, and Uwe D. Hanebeck,
    % Filtering on the Unit Sphere Using Spherical Harmonics
    % Proceedings of the 2017 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI 2017),
    % Daegu, Korea, November 2017.
    properties %#ok<*PROP>
        coeffMat;
        transformation;
    end
    
    methods
        function this = AbstractSphericalHarmonicsDistribution(coeffMat, transformation)
            assert(size(coeffMat, 2) == size(coeffMat, 1)*2-1, 'CoefficientMatrix:Size', 'Dimensions of coefficient Matrix are incompatible.');
            this.dim = 3;
            % Ignore irrelevant entries of coeffMat and set to NaN
            this.coeffMat = coeffMat + [[zeros(size(coeffMat, 1)-1, 1), ...
                kron(triu(NaN(size(coeffMat, 1)-1)), [1, 1])]; ...
                zeros(1, 2*size(coeffMat, 1)-1)];
            this.transformation = transformation;
            this = this.normalize;
        end
        function shd = convolve(this, other)
            % For the convolution, we only permit zonal densities as
            % "other". In our library, a density is seen as zonal if it is
            % zonal around the z-axis. Changes along
            % the z-axis are described by the coefficients of order m=0.
            % (This is a quite common convention).
            % Can use zonality test identically for sqrt transformation. If
            % density is zonal, its square root is as well.
            zonalTestMat = other.coeffMat;
            % The entries for m=0 don't matter for the zonality test, just
            % set them to NaN
            zonalTestMat(sub2ind(size(zonalTestMat), 1:size(zonalTestMat, 1), 1:size(zonalTestMat, 1))) = NaN;
            % If all other entires are NaN or almost zero (allow for numerical imprecision), then the
            % function is zonal.
            assert(all((zonalTestMat <= 1E-5) | isnan(zonalTestMat),[1,2]), 'Other is not zonal.');
            
            % We can truncate the harmonics to the smaller degree because
            % all entries for higher degrees would become zero
            if size(other.coeffMat, 1) < size(this.coeffMat, 1)
                this = this.truncate(size(other.coeffMat, 1)-1);
            elseif size(this.coeffMat, 1) < size(other.coeffMat, 1)
                other = other.truncate(size(this.coeffMat, 1)-1);
            end
            switch this.transformation
                case 'identity'
                    newCoeffMat = this.coeffMat .* ... % Get coefficients of other of order 0. Use implicit expansion.
                        (other.coeffMat(sub2ind(size(other.coeffMat), 1:size(other.coeffMat, 1), 1:size(other.coeffMat, 1))') ...
                        .* sqrt(4*pi./(2 * (0:size(other.coeffMat, 1) - 1) + 1))');% Times sqrt(4*pi/(l+1))
                    shd = this; % Do not use constructor to allow for inheritability
                    shd.coeffMat = newCoeffMat;
                case 'sqrt'
                    % First pad to degree 2L, so multiplication in state
                    % space does not void nonnegativity of the
                    % function with corresponding spherical harmonics
                    % coefficients. Currently only supported for
                    % SphericalHarmonicsDistributionComplex
                    assert(isa(this, 'SphericalHarmonicsDistributionComplex'));
                    degree = size(this.coeffMat, 1) - 1;
                    thisId = this.transformViaCoefficients('square', 2*degree);
                    switch other.transformation
                        case 'sqrt'
                            otherId = other.transformViaCoefficients('square', 2*degree);
                        case 'identity'
                            % We rely on the user to have ensured that this
                            % density does not attain negative function
                            % values
                            otherId = other;
                        otherwise
                            error('Transformation not supported');
                    end
                    resultIdFromSqrt = thisId.convolve(otherId);
                    % Truncate again in transformation to sqrt
                    shd = resultIdFromSqrt.transformViaCoefficients('sqrt', degree);
                otherwise
                    error('Transformation not supported');
            end
            
        end
        
        function shd = normalize(this)
            
            %%
            % Do not use constructor to allow for inheritance and avoid
            % another normalization iteration
            intVal = integralAnalytical(this);
            shd = this;
            if intVal < 0
                warning('Normalization:negative', 'Coefficient for first degree is negative. This can either be caused by a user error or due to negativity caused by non-square rooted version');
            elseif abs(intVal) < 1e-12
                error('Normalization:almostZero', 'Coefficient for first degree is too close to zero, this usually points to a user error');
            elseif abs(intVal-1) > 1e-5 % Warning massively impedes speed. Comment out for best performance.
                warning('Normalization:notNormalized', 'Coefficients apparently do not belong to normalized density. Normalizing...');
            else
                return % Normalized, return original density
            end
            if strcmp(this.transformation, 'identity')
                shd.coeffMat = this.coeffMat / intVal;
            elseif strcmp(this.transformation, 'sqrt')
                shd.coeffMat = this.coeffMat / sqrt(intVal);
            else
                warning('Currently cannot normalize');
            end
        end
        function intVal = integralAnalytical(this)
            switch this.transformation
                case 'identity'
                    intVal = this.coeffMat(1, 1) * sqrt(4*pi);
                case 'sqrt'
                    intVal = norm(this.coeffMat(~isnan(this.coeffMat)))^2;
                otherwise
                    error('No analytical formula for normalization available');
            end
            
        end
        
        function shd = truncate(this, degree)
            shd = this;
            if size(this.coeffMat, 1) - 1 > degree
                shd.coeffMat = this.coeffMat(1:degree+1, 1:2*degree+1);
            elseif size(this.coeffMat, 1) - 1 < degree
                warning('Truncate:TooFewCoefficients', 'Less coefficients than desired, filling up with zeros')
                % Do not use constructor to allow for inheritance and avoid
                % unnecessary checks (normalization is always preserved)
                shd.coeffMat = zeros(degree+1, 2*degree+1);
                shd.coeffMat(1:size(this.coeffMat, 1), 1:2*size(this.coeffMat, 1)-1) = this.coeffMat;
            end
        end
        function vals = pdf(this, xa)
            if strcmp(this.transformation, 'identity')
                vals = this.value(xa);
            elseif strcmp(this.transformation, 'sqrt')
                vals = this.value(xa).^2;
            else
                error('Transformation currently not supported')
            end
        end
    end
    methods(Abstract)
        vals = value(this, xa)
    end
    
end
