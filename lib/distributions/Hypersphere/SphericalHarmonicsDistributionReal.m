classdef SphericalHarmonicsDistributionReal < AbstractSphericalHarmonicsDistribution
    methods
        function this = SphericalHarmonicsDistributionReal(coeffMat, transformation)
            arguments
                coeffMat {mustBeReal}
                transformation char = 'identity'
            end
            this@AbstractSphericalHarmonicsDistribution(coeffMat, transformation);
        end
        function vals = value(this, xa)
            [phi, theta] = cart2sph(xa(1, :), xa(2, :), xa(3, :));
            theta = -theta + pi / 2; % Convert to different convention for spherical coordinates
            
            function ylmRes = genylm(l, theta, phi)
                % Use 'sch' to take care of the factorials, sqrt(2) and (-1)^m
                % (which is used to obtain signless spherical harmonics). See
                % "Evaluation of the rotation matrices in the basis of real
                % spherical harmonics" by Miguel A. Blanco, M. Flórez, and
                % M. Bermejo and Chapter 4 of "Group Theoretical Techniques in
                % Quantum Chemistry" by C. D. H. Chrisholm for more details.
                legendreFactors = legendre(l, cos(theta), 'sch')';
                % We use 'sch' to take care of the factorials and (-1)^m.
                % additional factors depending only on l are dealt with
                % directly in the loop later
                ylmRes = [sin(phi'*abs(-l:1:-1)), ones(length(phi), 1), cos(phi'*(1:l))] ...
                    .* legendreFactors(:, [end:-1:1, 2:end]);% Legendre functions for -m to m
            end
            
            vals = zeros(1, size(xa, 2));
            for i = 0:size(this.coeffMat, 1) - 1
                vals = vals + sqrt((2 * i + 1)/(4 * pi)) * this.coeffMat(i+1, 1:2*i+1) * genylm(i, theta, phi)';
            end
        end
        function shd = toSphericalHarmonicsDistributionComplex(this)
            if ~strcmp(this.transformation, 'identity')
                error('Transformation currently not supported')
            end
            coeffMat = NaN(size(this.coeffMat));
            for l = 0:size(this.coeffMat, 1) - 1
                coeffMat(l+1, l+1+(-l:-1)) = (1i * this.coeffMat(l+1, l+1+(-l:-1))) / sqrt(2) + this.coeffMat(l+1, l+1-(-l:-1)) / sqrt(2);
                coeffMat(l+1, l+1) = this.coeffMat(l+1, l+1);
                coeffMat(l+1, l+1+(1:l)) = (-1).^(1:l) .* ((-1i * this.coeffMat(l+1, l+1-(1:l))) + this.coeffMat(l+1, l+1+(1:l))) / sqrt(2);
            end
            shd = SphericalHarmonicsDistributionComplex(coeffMat, this.transformation);
        end
        function shd = rotate(this, alpha, beta, gamma)
            % Currently, a conversion to complex is used to use the
            % rotation implemented there
            shdComplex = this.toSphericalHarmonicsDistributionComplex;
            shdComplexRot = shdComplex.rotate(alpha, beta, gamma);
            shd = shdComplexRot.toSphericalHarmonicsDistributionReal;
        end
        function mu = meanDirection(this)
            if ~strcmp(this.transformation, 'identity')
                error('Transformation currently not supported')
            end
            if ~(size(this.coeffMat) > 1)
                error('Too few coefficients available to calculate the mean');
            elseif norm(this.coeffMat(2, [3, 1, 2])) < 1E-8
                error('Coefficients of degree 1 are almost zero. Therefore, no meaningful mean is available');
            else
                mu = this.coeffMat(2, [3, 1, 2])' ./ norm(this.coeffMat(2, [3, 1, 2]));
            end
        end
    end
    
end