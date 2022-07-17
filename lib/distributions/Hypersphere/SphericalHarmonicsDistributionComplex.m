classdef SphericalHarmonicsDistributionComplex < AbstractSphericalHarmonicsDistribution
    % Class for representing a density using a complex spherical harmonic.
    % See the value function for the convention used.
    
    % see Florian Pfaff, Gerhard Kurz, and Uwe D. Hanebeck,
    % Filtering on the Unit Sphere Using Spherical Harmonics
    % Proceedings of the 2017 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI 2017),
    % Daegu, Korea, November 2017.
    methods
        function this = SphericalHarmonicsDistributionComplex(coeffMat, transformation)
            % coeffMat is assumed to be a L+1 x 2*L+1 matrix. The
            % coefficients need to be given as a lower triangular matrix
            % (including the diagonal). The coefficients of degree l are
            % given in the l+1th row in column 1 to 2*l+1. For better
            % legibility, all functions of this class set all upper
            % triangular parts to NaN or zero.
            
            if nargin == 1, transformation = 'identity';end
            assert(size(coeffMat, 2) == (2 * size(coeffMat, 1) - 1), 'Size of coeffMat is invalid.');
            assert(~anynan(tril(coeffMat)), 'Lower triangular part must not be NaN.');
            % Complex spherical harmonics can describe complex functions.
            % Therefore, we have to verify that the coefficients describe a
            % real function.
            for l = 0:size(coeffMat, 1) - 1
                if ~(max(abs(coeffMat(l+1, 2*l+1:-1:l+1)-(-1).^(l:-1:0).*conj(coeffMat(l+1, 1:l+1)))) < 1e-3)
                    warning('CoefficientMatrix:ConjugateIncorrect', ...
                        'Entries of the coefficient matrix do not describe a real function. This may happen sometimes for low numbers of coefficients. However, if it happens regularly, there may be a bug. Now ensuring it describes a real function again.');
                    coeffMat(l+1, 2*l+1:-1:l+1) = (-1).^(l:-1:0) .* conj(coeffMat(l+1, 1:l+1));
                end
            end
            this@AbstractSphericalHarmonicsDistribution(coeffMat, transformation);
        end
        function vals = value(this, xa)
            [phi, theta] = cart2sph(xa(1, :), xa(2, :), xa(3, :));
            theta = -theta + pi / 2; % Convert to different convention for spherical coordinates
            % We follow the convention that includes the Condon-Shortley
            % phase. See, e.g., "Mathematical Methods for Physicists, 7th
            % ed." by Arfken, Weber and Harris, Ch. 15.5.
            %
            % A table for the basis functions is available in the
            % aforementioned book in Table 15.4. The convention used
            % is in accordance with the convention used in
            % https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
            
            
            function ylmRes = genylm(l, theta, phi)
                % The associated legendre polynomials as obtained by simply
                % calling legendre in Matlab include the Condon-Shortley
                % phase. However, we use the Schmidt Seminormalized
                % Associated Legendre Functions as this way,
                % we do not have to address the factorials. However, in the
                % Schmidt Seminormalized form, the Condon-Shortley phase is
                % elminiate and an unwanted factor sqrt(2) is used on top.
                legendreFactors = legendre(l, cos(theta), 'sch')';
                
                % First, generate the negative spherical harmonics as no
                % Condon-Shortley phase applies to them. Also elminate
                % unwanted sqrt(2). Remaining factors of the normalization
                % are applied later on.
                ylmNeg = 1 / sqrt(2) * legendreFactors(:, end:-1:2) .* exp(1i*phi'*(-l:-1));
                % Use identity for negative spherical harmonics as given in
                % "Evaluation of the rotation matrices in the basis of real
                % spherical harmonics" by Blanco et al and "Hilbert Space
                % Methods in Signal Processing" Sec 7.3.3 by Kennedy et al.
                % This reintroduces the Condon-Shortley phase.
                ylmPos = (-1).^(1:l) .* conj(fliplr(ylmNeg));
                % Apply remaining factors and add spherical harmonic with order 0.
                ylmRes = sqrt((2 * l + 1)/(4 * pi)) * [ylmNeg, legendreFactors(:, 1), ylmPos];
            end
            
            vals = zeros(1, size(xa, 2));
            for l = 0:size(this.coeffMat, 1) - 1
                vals = vals + real((this.coeffMat(l+1, 1:2*l+1))*genylm(l, theta, phi).');
            end
        end
        function shd = toSphericalHarmonicsDistributionReal(this)
            if ~strcmp(this.transformation, 'identity')
                error('Transformation currently not supported')
            end
            coeffMat = NaN(size(this.coeffMat));
            for l = 0:size(this.coeffMat, 1) - 1
                coeffMat(l+1, l+1+(-l:-1)) = sqrt(2) * (-1).^((-l:-1) + 1) .* imag(this.coeffMat(l+1, l+1-(-l:-1)));
                coeffMat(l+1, l+1) = this.coeffMat(l+1, l+1);
                coeffMat(l+1, l+1+(1:l)) = sqrt(2) * (-1).^(1:l) .* real(this.coeffMat(l+1, l+1+(1:l)));
            end
            assert(abs(max(imag(coeffMat(:))))<1e-10)
            shd = SphericalHarmonicsDistributionReal(real(coeffMat), this.transformation);
        end
        function f = transformViaCoefficients(this, desiredTransformation, degree)
            if degree > size(this.coeffMat, 1) - 1
                % Do not use truncate to skip unnecessary operations such
                % as the normalization and also avoid warning for too few c
                % coeffs.
                coeffMatExtended = zeros(degree+1, 2*degree+1);
                coeffMatExtended(1:size(this.coeffMat, 1), 1:size(this.coeffMat, 2)) = this.coeffMat;
                this.coeffMat = coeffMatExtended;
            end
            if strcmp(this.transformation, 'identity') && strcmp(desiredTransformation, 'sqrt')
                % Calculate values on a grid, calculate square root, transform to coefficients
                plmNew = xyz2plm(sqrt(plm2xyz(this.plmWithDegreeAndOrder)));
                f = SphericalHarmonicsDistributionComplex.fromPlmWithDegreeAndOrder( ...
                    ... % Truncate plm before creating distribution
                    plmNew(1:min((degree + 1)*(degree + 2)/2, size(plmNew, 1)), :), ...
                    'sqrt');
            elseif strcmp(this.transformation, 'sqrt') && strcmp(desiredTransformation, 'square')
                % Calculate values on a grid, square values, transform to
                % coefficients.
                plmNew = xyz2plm((plm2xyz(this.plmWithDegreeAndOrder)).^2);
                f = SphericalHarmonicsDistributionComplex.fromPlmWithDegreeAndOrder( ...
                    ... % Truncate plm before creating distribution
                    plmNew(1:min((degree + 1)*(degree + 2)/2, size(plmNew, 1)), :), ...
                    'identity');
            else
                error('Transformation currently not supported');
            end
        end
        function plm = plmWithDegreeAndOrder(this)
            % Returns plm with order and degree (as required by slepian)
            L = size(this.coeffMat, 1) - 1;
            [mVec, lVec] = addmon(L);
            
            coeffs = NaN(numel(mVec), 2);
            index = 1;
            for l = 0:L % Do for all m of one l at once
                coeffs(index:index+l, 1) = real(this.coeffMat(l+1, l+1:-1:1));
                coeffs(index:index+l, 2) = imag(this.coeffMat(l+1, l+1:-1:1));
                index = index + l + 1;
            end
            coeffs([1, 1 + cumsum(1:L)], 1) = coeffs([1, 1 + cumsum(1:L)], 1) / sqrt(2); % Entries for m=0 have to be adjusted
            
            plm = [lVec, mVec, coeffs / sqrt(2*pi)];
        end
        function shd = multiply(this, other, degree)
            if nargin == 2
                degree = size(this.coeffMat, 1) - 1;
            end
            assert(strcmp(this.transformation, other.transformation), 'Transformation has to be identical to perform multiplication');
            if (degree + 1) > size(this.coeffMat, 1) + size(other.coeffMat, 1) - 1
                warning('Multiplication:degreeTooHigh', 'A degree for which no entries could possible be nonzero was requested on the call of multiply. Filling up with zeros...');
            end
            
            % Pad to respect higher frequencies (this could be made
            % variable)
            warningSettings = warning('off', 'Truncate:TooFewCoefficients');
            thistrun = this.truncate(size(this.coeffMat, 1)+size(other.coeffMat, 1)-2);
            othertrunc = other.truncate(size(this.coeffMat, 1)+size(other.coeffMat, 1)-2);
            warning(warningSettings);
            % plmSlepianToLibD is defined below the class definition
            plm = xyz2plm(plm2xyz(thistrun.plmWithDegreeAndOrder).*plm2xyz(othertrunc.plmWithDegreeAndOrder));
            % Normalization is left up to the constructor but is expected
            warningSettings = warning('off', 'Normalization:notNormalized');
            shd = SphericalHarmonicsDistributionComplex.fromPlmWithDegreeAndOrder(plm, this.transformation);
            shd = shd.truncate(degree);
            warning(warningSettings);
        end
        function shd = multiplyViaCoefficients(this, other, degree)
            if nargin == 2
                degree = size(this.coeffMat, 1) - 1;
            end
            assert(strcmp(this.transformation, other.transformation), 'Transformation has to be identical to perform multiplication');
            if (degree + 1) > size(this.coeffMat, 1) + size(other.coeffMat, 1) - 1
                warning('Multiplication:degreeTooHigh', 'A degree for which no entries could possible be nonzero was requested on the call of multiply. Filling up with zeros...');
            end
            
            coeffMat = zeros(degree+1, 2*degree+1);
            % For the multiplication formula,
            % see "Computing Fourier Transforms
            % and Convolutions on the 2-Sphere"
            % by Driscoll and Healy. For the
            % Wiegner Coefficients, the
            % convention from the book "Angular
            % momentum in quantum physics" by
            % Biedernharn and Louck is used.
            % This convention, given in Ch.
            % 3.12, is slightly different from
            % Wiegner3j.m in the expected sign
            % of the last m entry and a factor
            % of sqrt(2*L+1). Since two Wigner
            % coefficients are multiplied, we
            % have to multiply the formula by
            % (2*L+1)
            function Lm1m2 = LCellFill(m1, m2) % To determine all relevant L for a kombination of l1, l2, m1, m2
                L = max(abs(l1-l2), abs(m1+m2)):min(l1+l2, degree);
                Lm1m2 = {[L; m1 * ones(1, numel(L)); m2 * ones(1, numel(L))]};
            end
            [l1Mesh, l2Mesh, LMesh] = ndgrid(0:size(this.coeffMat, 1)-1, 0:size(other.coeffMat, 1)-1, 0:min(size(this.coeffMat, 1)+size(other.coeffMat, 1), degree));
            wignerl1l2L = reshape(threej(l1Mesh(:), l2Mesh(:), LMesh(:), 0, 0, 0), size(l1Mesh));
            for l1 = 0:size(this.coeffMat, 1) - 1
                for l2 = 0:size(other.coeffMat, 1) - 1
                    % First determin all necessary L to get the Wigner
                    % coefficients for all m1xm2xL simultaneously
                    LCell = cell(2*l1+1, 2*l2+1);
                    
                    % Filter out entries for which both
                    % relevant coefficients are <1e-8
                    m1bigger = abs(this.coeffMat(l1+1, l1+1+(-l1:l1)))' > 1e-8;
                    m2bigger = abs(other.coeffMat(l2+1, l2+1+(-l2:l2))) > 1e-8;
                    [m1Mesh, m2Mesh] = meshgrid(-l1:l1, -l2:l2);
                    m1Mesh = m1Mesh';
                    m2Mesh = m2Mesh';
                    relevant = repmat(m1bigger, 1, numel(m2bigger)) & repmat(m2bigger, numel(m1bigger), 1);
                    
                    if ~any(relevant(:))
                        continue
                    end
                    % Determine for which there is any L to consider
                    m1MeshRel = m1Mesh(relevant);
                    m2MeshRel = m2Mesh(relevant);
                    LCell(relevant) = arrayfun(@LCellFill, m1MeshRel, m2MeshRel);
                    relevant = relevant & ~cellfun(@isempty, LCell);
                    if ~any(relevant(:))
                        continue
                    end
                    LAll = [LCell{:}];
                    % Get all necessary Wigner coefficients for this
                    % combination of l1 and l2
                    wignerCoeffs = threej(l1*ones(1, size(LAll, 2)), l2*ones(1, size(LAll, 2)), LAll(1, :), LAll(2, :), LAll(3, :), -LAll(2, :)-LAll(3, :));
                    wignertmp = wignerl1l2L(l1+1, l2+1, LAll(1, :)+1);
                    wignerCoeffs = wignerCoeffs .* wignertmp(:)';
                    
                    LIndex = 1;
                    for i = 1:numel(m1MeshRel)
                        m1 = m1MeshRel(i);
                        m2 = m2MeshRel(i);
                        LRange = max(abs(l1-l2), abs(m1+m2)):min(l1+l2, degree); % Could also be obtained from LAll(1,lIndex:lIndex+numel(LRange)-1))
                        indices = sub2ind(size(coeffMat), LRange+1, LRange+1+m1+m2);
                        coeffMat(indices) = coeffMat(indices) + (-1)^(m1 + m2) * sqrt((2 * l1 + 1)*(2 * l2 + 1)*(2 * LRange + 1)/(4 * pi)) ...
                            .* wignerCoeffs(LIndex:LIndex+numel(LRange)-1) ...
                            * (this.coeffMat(l1+1, l1+1+m1) .* other.coeffMat(l2+1, l2+1+m2));
                        LIndex = LIndex + numel(LRange);
                    end
                end
            end
            % Normalization is left up to the constructor but is expected
            warningSettings = warning('off', 'Normalization:notNormalized');
            shd = SphericalHarmonicsDistributionComplex(coeffMat, this.transformation);
            warning(warningSettings);
        end
        function shd = rotate(this, alpha, beta, gamma)
            % The function is rotated using the z-y-z
            % convention by alpha, beta, and gamma.
            % First, the Wigner (small) d Matrix is obtained using the naive
            % closed form formula given in eq. 37 of "Evaluation of the rotation
            % matrices in the basis of real spherical harmonics"
            % by Miguel A. Blanco, M. Flórez, and M. Bermejo.
            dCell = arrayfun(@(l){zeros(l+1, 2*l+1)}, 0:size(this.coeffMat, 1)-1); % Initizalize
            factTable = factorial(0:2*size(this.coeffMat, 1)); % precompute factorials
            for l = 0:size(this.coeffMat, 1) - 1
                for mNew = -l:l
                    factTemp1 = factTable(l+mNew+1);
                    factTemp2 = factTable(l-mNew+1);
                    for mOrig = -l:l
                        s = max(0, -(mNew - mOrig)):min(l-mNew, l+mOrig);
                        dCell{l+1}(l + 1 + mNew, l + 1 + mOrig) = (-1)^(mNew - mOrig) * sqrt(factTemp1*factTemp2*factTable(l+mOrig+1)*factTable(l-mOrig+1)) ...
                            * sum((-1).^s.*cos(beta/2).^(2 * (l - s) + mOrig - mNew).*sin(beta/2).^(2 * s - mOrig + mNew) ...
                            ./(factTable(l+mOrig-s+1) .* factTable(s+1) .* factTable(mNew-mOrig+s+1) .* factTable(l-mNew-s+1)));
                    end
                end
            end
            D = cell(size(dCell));
            % Calculcate Wigner (capital) D from above result. See eq. 36
            % of the above mentioned paper.
            for l = 0:length(D) - 1
                expmatmOrig = repmat((-l:l), [size(dCell{l+1}, 1), 1]);
                expmatmNew = repmat((-l:l)', [1, size(dCell{l+1}, 2)]);
                D{l+1} = exp(-1i*expmatmOrig*alpha-1i*expmatmNew*gamma) .* dCell{l+1};
            end
            coeffMatRot = NaN(size(this.coeffMat), 'like', 1i);
            for l = 0:size(this.coeffMat, 1) - 1
                for mNew = -l:l
                    coeffMatRot(l+1, l+1+mNew) = sum(this.coeffMat(l+1, 1:2*l+1).*D{l+1}(l + 1 + mNew, :));
                end
            end
            shd = this;
            shd.coeffMat = coeffMatRot; % Do not use constructor to skip normalization test
        end
        function mu = meanDirection(this)
            if ~(size(this.coeffMat) > 1)
                error('Too few coefficients available to calculate the mean');
            end
            y = imag(this.coeffMat(2, 1)+this.coeffMat(2, 3)) / sqrt(2);
            x = real(this.coeffMat(2, 1)-this.coeffMat(2, 3)) / sqrt(2);
            z = real(this.coeffMat(2, 2));
            if norm([x; y; z]) < 1E-9
                error('Coefficients of degree 1 are almost zero. Therefore, no meaningful mean is available');
            end
            mu = [x; y; z] ./ norm([x; y; z]);
        end
    end
    methods(Static)
        function shd = fromDistributionViaIntegral(dist, degree, transformation)
            if nargin == 2 % Default to identity
                transformation = 'identity';
            end
            % This is just a convenience function to convert distributions
            % to complex spherical harmonics. No closed form solutions are
            % available (unlike in the FourierDistribution case)
            assert(isa(dist, 'AbstractHypersphericalDistribution') && dist.dim == 3, 'dist must be a distribution on the sphere.');
            shd = SphericalHarmonicsDistributionComplex.fromFunctionViaIntegralCart(@(x)dist.pdf(x), degree, transformation);
        end
        function shd = fromDistribution(dist, degree, transformation)
            if nargin == 2 % Default to identity
                transformation = 'identity';
            end
            if isa(dist,'SphericalGridDistribution')
                shd = SphericalHarmonicsDistributionComplex.fromGrid(this.fvals, this.grid, transformation);
            else
                shd = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(dist, degree, transformation);
            end
        end
        function shd = fromDistributionNumericalFast(dist, degree, transformation)
            % Just a convenience function to call from function
            if nargin == 2 % Default to identity
                transformation = 'identity';
            end
            assert(isa(dist, 'AbstractHypersphericalDistribution') && dist.dim == 3, 'dist must be a distribution on the sphere.');
            shd = SphericalHarmonicsDistributionComplex.fromFunctionFast(@dist.pdf, degree, transformation);
        end
        function shd = fromFunctionFast(fun, degree, transformation)
            if nargin == 2 % Default to identity
                transformation = 'identity';
            end
            % Converts a real-valued function depending on cartesian
            % coordinates to a complex spherical harmonics representation
            % The function needs to be able to handle a matrix with three
            % rows (x, y and z coordindates) as the input argument.
            % This is the case for all .pdf functions on the sphere.
            assert(degree >= 1);
            lat = linspace(0, 2*pi, 2*degree+2);
            lon = linspace(pi/2, -pi/2, degree+2);
            [latMesh, lonMesh] = meshgrid(lat, lon);
            fval = funSph(latMesh, lonMesh);
            switch transformation
                case 'identity' % Do not need to transform
                case 'sqrt'
                    fval = sqrt(fval);
                otherwise
                    error('Currently, only identity transformation is supported');
            end
            shd = SphericalHarmonicsDistributionComplex.fromPlmWithDegreeAndOrder(xyz2plm(fval), transformation);
            function vals = funSph(theta, phi)
                [x, y, z] = sph2cart(theta, phi, 1);
                vals = reshape(fun([x(:)'; y(:)'; z(:)']), size(theta));
            end
        end
        function shd = fromFunctionViaIntegralCart(fun, degree, transformation)
            if nargin == 2 % Default to identity
                transformation = 'identity';
            end
            % Converts a real-valued function depending on cartesian
            % coordinates to a complex spherical harmonics representation
            % The function needs to be able to handle a matrix with three
            % rows (x, y and z coordindates) as the input argument.
            % This is the case for all .pdf functions on the sphere.
            shd = SphericalHarmonicsDistributionComplex.fromFunctionViaIntegralSph(@funSph, degree, transformation);
            function vals = funSph(theta, phi)
                theta = -theta + pi / 2; % Convert to different convention for spherical coordinates
                [x, y, z] = sph2cart(phi, theta, 1);
                vals = fun([x; y; z]);
            end
        end
        function shd = fromFunctionViaIntegralSph(fun, degree, transformation)
            if nargin == 2 % Default to identity
                transformation = 'identity';
            end
            % Converts a real-valued function depending on spherical
            % coordinates to a complex spherical harmonics representation
            % Expects a function taking theta (elevation from 0 to pi) and
            % phi (azimuth) as two separate row vectors
            
            % This functionality is EXPERIMENTAL and may require further
            % testing.
            
            if strcmp(transformation, 'sqrt')
                fun = @(theta, phi)sqrt(fun(theta, phi));
            elseif ~strcmp(transformation, 'identity')
                error('Transformation not supported')
            end
            
            % The unwanted sqrt(2) is directly undone here.
            ylmTrunc = @(l, m, theta, phi)((m ~= 0) / sqrt(2) + (m == 0)) * sqrt((2 * l + 1)/(4 * pi)) * exp(1i*phi'*m) ...
                ... % Eliminate all except for the necessary one at abs(m)
                .* (legendre(l, cos(theta), 'sch')' * [zeros(abs(m), 1); 1; zeros(l-abs(m), 1)]);
            
            coeffMat = NaN(degree+1, 2*degree+1, 'like', 1i);
            for l = 0:degree
                for m = -l:l
                    % Attention: the ' transposes and performs the
                    % conjugate (this behavior is explicitly desired)
                    coeffMat(l+1, m+l+1) = ((m <= 0) + (m > 0) * (-1)^m) * integral2(@(phi, theta)reshape(fun(theta(:)', phi(:)') ...
                        .*ylmTrunc(l, m, theta(:)', phi(:)')'.*sin(theta(:)'), size(theta)), 0, 2*pi, 0, pi);
                end
            end
            shd = SphericalHarmonicsDistributionComplex(coeffMat, transformation);
        end
        function shd = fromGrid(gridValues, grid, transformation, degree)
            % Values are assumed to be given without any transformation.
            % If no grid is given (i.e., only one input or an empty grid),
            % a grid that is directly compatible with
            % the spherical harmonics transform is assumed. For other
            % regular grids, directly provide grid (it is treated as
            % irregular grid then).
            arguments
                gridValues (:,:) double
                grid (3,:) double = zeros(3,0)
                transformation char = 'identity';
                degree (1,1) double = (-6+sqrt(36-8*(4-numel(gridValues))))/4;
            end
            if strcmp(transformation, 'sqrt')
                gridValues=sqrt(gridValues);
            end
            if nargin==1 || isempty(grid)
                assert(size(gridValues,2)>1, 'For regular grids, provide values as a matrix.');
                assert(nargin==3 || floor(degree)==degree, ...
                    'Based on the number of values, this grid is definitely not directly compatible with spherical harmonics transform.');
                shd = SphericalHarmonicsDistributionComplex.fromPlmWithDegreeAndOrder(xyz2plm(reshape(gridValues,degree+2,2*degree+2)), transformation);
            else
                assert(size(gridValues,1)==numel(gridValues) && size(gridValues,1)==size(grid,2));
                [lon,lat]=cart2sph(grid(1,:),grid(2,:),grid(3,:));
                lon=rad2deg(lon);
                lat=rad2deg(lat);
                shd = SphericalHarmonicsDistributionComplex.fromPlmWithDegreeAndOrder(...
                    xyz2plm(gridValues(:),ceil(degree),'irr',lat(:),lon(:)), transformation);
            end
        end
        function shd = fromPlmWithDegreeAndOrder(plm, transformation)
            coeffMat = NaN(plm(end, 1)+1, 2*plm(end, 1)+1, 'like', 1i);
            coeffMat(sub2ind(size(coeffMat), plm(:, 1)+1, plm(:, 1)+1-plm(:, 2))) = plm(:, 3) .* sqrt(2).^(plm(:, 2) == 0) + plm(:, 4) * 1i;
            coeffMat(sub2ind(size(coeffMat), plm(:, 1)+1, plm(:, 1)+1+plm(:, 2))) = (-1).^plm(:, 2) .* plm(:, 3) .* sqrt(2).^(plm(:, 2) == 0) ...
                +(-1).^(plm(:, 2) + 1) .* plm(:, 4) * 1i;
            coeffMat = sqrt(2*pi) * coeffMat;
            shd = SphericalHarmonicsDistributionComplex(coeffMat, transformation);
        end
    end
end
