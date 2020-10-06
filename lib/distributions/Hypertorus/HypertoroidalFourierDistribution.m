classdef HypertoroidalFourierDistribution < AbstractHypertoroidalDistribution
    % Distribution on the Hypertorus using Fourier series
    %
    % Florian Pfaff, Gerhard Kurz, Uwe D. Hanebeck,
    % Multivariate Angular Filtering Using Fourier Series
    % Journal of Advances in Information Fusion, 11(2):206-226, December 2016.
    
    properties
        % Saving the whole matrix may not be numerically the smartest way to
        % handle it but is really easy. Currently, the symmetric argument in
        % ifft ensures that the two parts do not diverge indefinitely
        C double
        transformation char
    end
    
    methods
        function this = HypertoroidalFourierDistribution(C, transformation)
            % We assume this is at least twodimensional. Use
            % FourierDistribution for the one dimensional case
            if isa(C, 'AbstractHypertoroidalDistribution')
                error('fourierCoefficients:invalidCoefficientMatrix', 'You gave a distribution as the first argument. To convert distributions to a distribution in Fourier representation, use .fromDistribution');
            elseif numel(C) == 1
                warning('fourierCoefficients:singleCoefficient', 'Fourier series only has one element, assuming dimension 1.');
            end
            this.dim = ndims(C) - (ismatrix(C) && size(C, 2) == 1); % Correction term for 1D
            if nargin == 1 % Square root of density is standard case
                this.transformation = 'sqrt';
            else
                this.transformation = transformation;
            end
            this.C = C;
            % Check if normalized. If not: Normalize!
            this = this.normalize;
        end
        
        function p = value(this, xa)
            assert(size(xa, 1) == this.dim);
            maxk = (size(this.C) - 1) / 2;
            kRanges = arrayfun(@(currMaxk){-currMaxk:currMaxk}, maxk);
            individualIndexMatrices = cell(size(kRanges));
            p = NaN(1, size(xa, 2), 'like', 1i); % Intermediate results will be complex due to numerical inaccuracy
            [individualIndexMatrices{:}] = ndgrid(kRanges{:});
            % Alternative implementation to the one below that are usually
            % not faster because of higher memory use.
            % 1) With implicit expansion
            % exponents=cat(size(xa,1)+1,individualIndexMatrices{:}).*reshape(xa,[ones(1,size(xa,1)),size(xa)]);
            % multWithC=this.C.*exp(sum(exponents,size(xa,1)+1)*1i);
            % p=sum(reshape(multWithC,[],size(xa,2)),1);
            % 2) With bsxfun
            % exponents=bsxfun(@times,cat(size(xa,1)+1,individualIndexMatrices{:}),reshape(xa,[ones(1,size(xa,1)),size(xa)]));
            % multWithC=bsxfun(@times,exp(sum(exponents,size(xa,1)+1)*1i),this.C);
            % p=sum(reshape(multWithC,[],size(xa,2)),1);
            % 3) Implementation using loops that is actually used
            for i = 1:size(xa, 2)
                matCurr = individualIndexMatrices{1} * xa(1, i);
                for currDim = 2:this.dim
                    matCurr = matCurr + individualIndexMatrices{currDim} * xa(currDim, i);
                end
                matCurr = this.C .* exp(matCurr*1i);
                p(i) = sum(matCurr(:));
            end
            assert(all(imag(p) < 0.1));
            p = real(p);
        end
        
        function p = pdf(this, xa)
            % Evaluate pdf at each column of xa
            %
            % Parameters:
            %   xa (2 x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            val = value(this, xa);
            switch this.transformation
                case 'sqrt'
                    assert(all(imag(val) < 0.0001));
                    p = real(val).^2;
                case 'identity'
                    p = val;
                case 'log'
                    warning('pdf:mayNotBeNormalized', 'Density may not be normalized');
                    p = exp(val);
                otherwise
                    error('transformation:unrecognizedTransformation', 'Transformation not recognized or unsupported');
            end
        end
        
        function hfd = truncate(this, noOfCoefficients, forceNormalization)
            % Truncates Fourier series. Fills up if there are less coefficients
            % Expects number of complex coefficients (or sum of number of real
            % coefficients)
            if nargin < 3
                forceNormalization = false;
            end
            if numel(noOfCoefficients) == 1 && this.dim ~= 1
                noOfCoefficients = ones(1, this.dim) * noOfCoefficients; % If only one value given, truncate evenly
            elseif numel(noOfCoefficients) == 1 && this.dim == 1
                noOfCoefficients = [noOfCoefficients, 1]; % To handle 1D case
            end
            assert(all(mod(noOfCoefficients-1, 2) == 0) && (all((noOfCoefficients > 1)) || this.dim == 1), ...
                'Invalid number of coefficients, numbers for all dimensions have to be odd and greater than 1.');
            if isequal(size(this.C), noOfCoefficients) % Do not need to truncate as already correct size
                if ~forceNormalization
                    hfd = this;
                else
                    warnstruct = warning('off', 'Normalization:notNormalized'); % We want to normalize, do not throw warning
                    hfd = this.normalize;
                    warning(warnstruct);
                end
                return
            elseif any(size(this.C) < noOfCoefficients)
                warning('Truncate:TooFewCoefficients', 'At least in one dimension, truncate has to fill up due to too few coefficients.')
            end
            
            maxksold = (size(this.C) - 1) / 2;
            maxksnew = (noOfCoefficients - 1) / 2;
            Cnew = zeros(noOfCoefficients, 'like', 1i);
            indicesNewMat = arrayfun(@(i){max(maxksnew(i)-maxksold(i)+1, 1):min(maxksnew(i)+maxksold(i)+1, size(Cnew, i))}, 1:this.dim);
            indicesOldMat = arrayfun(@(i){max(maxksold(i)-maxksnew(i)+1, 1):min(maxksold(i)+maxksnew(i)+1, size(this.C, i))}, 1:this.dim);
            Cnew(indicesNewMat{:}) = this.C(indicesOldMat{:});
            hfd = this;
            hfd.C = Cnew;
            % Truncation can void normalization if transformation is not
            % identity. Normalize only if enforced by the uesr or if
            % truncation can void normalization.
            if forceNormalization || (~strcmp(this.transformation, 'identity')) && any(noOfCoefficients < size(this.C))
                % Disable warning as we expect normalization to be
                % necessary
                warnStruct = warning('off', 'Normalization:notNormalized');
                hfd = hfd.normalize;
                warning(warnStruct);
            end
        end
        
        function f = multiply(this, f2, noOfCoefficients)
            assert(strcmp(this.transformation, f2.transformation));
            if nargin == 2
                if ismatrix(this.C)&&size(this.C,2)==1 % 1-D case
                    noOfCoefficients=numel(this.C);
                else
                    noOfCoefficients = size(this.C);
                end
            elseif numel(noOfCoefficients) == 1
                noOfCoefficients = ones(1, this.dim) * noOfCoefficients;
            end
            if strcmp(this.transformation, 'log')
                f = this;
                if ~isequal(noOfCoefficients, size(f2.C))
                    this = this.truncate(noOfCoefficients);
                end
                if ~isequal(noOfCoefficients, size(f2.C))
                    f2 = f2.truncate(noOfCoefficients);
                end
                f.C = this.C + f2.C;
                warning('Multiply:NotNormalizing', 'Not performing normalization when using log transformation.');
            elseif strcmp(this.transformation, 'identity') || strcmp(this.transformation, 'sqrt')
                f = this;
                if all(noOfCoefficients <= size(this.C)) % Do not calculate additional coefficients if they are not needed
                    f.C = convnc(this.C, f2.C, 'same');
                else
                    f.C = convnc(this.C, f2.C, 'full');
                end
            else
                error('Multiply:unsupportedTransformation', 'Transformation not recognized or unsupported.');
            end
            % Normalization is intended.
            f = f.truncate(noOfCoefficients, true);
        end
        
        function hfd = convolve(this, f2, noOfCoefficients)
            assert(strcmp(this.transformation, f2.transformation));
            if nargin == 2
                if ismatrix(this.C)&&size(this.C,2)==1 % 1-D case
                    noOfCoefficients=numel(this.C);
                else
                    noOfCoefficients = size(this.C);
                end
            elseif numel(noOfCoefficients) == 1
                noOfCoefficients = ones(1, this.dim) * noOfCoefficients;
            end
            if ~isequal(size(this.C), size(f2.C)) % Adjust coefficient matricies if they are not equal
                warnStruct = warning('off', 'Truncate:TooFewCoefficients');
                hfdtmp = this.truncate(max([size(this.C); size(f2.C)]));
                f2 = f2.truncate(max([size(this.C); size(f2.C)]));
                warning(warnStruct);
            else
                hfdtmp = this;
            end
            switch this.transformation
                case 'sqrt'
                    % May not be truncated as it would void nonnegativity.
                    % Also note that the convolution is not distributive in
                    % regard to the multiplication.
                    cConv = (2 * pi)^(hfdtmp.dim) * convnc(hfdtmp.C, hfdtmp.C, 'full') .* convnc(f2.C, f2.C, 'full');
                    hfdtmp.C = cConv;
                    hfdtmp.transformation = 'identity';
                    hfd = hfdtmp.transformViaFFT('sqrt', noOfCoefficients);
                    hfd = hfd.truncate(noOfCoefficients); % Normalization is performed in truncation
                case 'identity'
                    cConv = (2 * pi)^(this.dim) * hfdtmp.C .* f2.C;
                    hfdtmp.C = cConv;
                    hfdtmp.transformation = 'identity';
                    hfd = hfdtmp.truncate(noOfCoefficients, true); % Enforce normalization just to be safe in case numerical issues occurred
                otherwise
                    error('transformation:unrecognizedTransformation', 'Transformation not recognized or unsupported');
            end
        end
        
        function hfd = normalize(this)
            % Normalize Fourier density while taking its type into account
            tol = 1e-4;
            switch this.transformation
                case 'sqrt'
                    c00 = norm(this.C(:))^2; % Square root calculated later to use norm and not squared norm
                    factorForId = c00 * (2 * pi)^(this.dim);
                    normalizationFactor = sqrt(factorForId);
                case 'identity'
                    % This will always get the most central element that
                    % corresponds to c00.
                    c00 = this.C((numel(this.C) + 1)/2);
                    factorForId = c00 * (2 * pi)^(this.dim);
                    normalizationFactor = factorForId;
                otherwise
                    warning('Normalization:cannotTest', 'Unable to test if normalized');
                    hfd = this;
                    return
            end
            hfd = this;
            if c00 < 0
                warning('Normalization:negative', 'C00 is negative. This can either be caused by a user error or due to negativity caused by non-square rooted version');
            elseif abs(c00) < 1e-200 % Tolerance has to be that low to avoid unnecessary errors on multiply
                error('Normalization:almostZero', 'C00 is too close to zero, this usually points to a user error');
            elseif abs(factorForId-1) > tol
                warning('Normalization:notNormalized', 'Coefficients apparently do not belong to normalized density. Normalizing...');
            else
                return % Normalized, return original density
            end
            hfd.C = this.C / normalizationFactor;
        end
        
        function hfd = marginalizeOut(this, dimension)
            assert(strcmp(this.transformation, 'identity'), 'Marginals are currently only supported for the identity transformation')
            assert(dimension >= 1 || dimension <= this.dim);
            hfd = this;
            
            indices = repmat({':'}, 1, this.dim);
            indices{dimension} = (size(this.C, dimension) + 1) / 2;
            coeffs = size(this.C);
            coeffs(dimension) = [];
            if numel(coeffs) == 1 % Return a FourierDistribution if only one dimension remains
                hfd = FourierDistribution.fromComplex(reshape(2*pi*this.C(indices{:}), [1, coeffs]), 'identity');
                return
            end
            hfd.C = reshape(2*pi*this.C(indices{:}), coeffs);
            hfd.dim = this.dim - 1;
        end
        
        function fd = marginalizeTo1D(this, dimension)
            assert(strcmp(this.transformation, 'identity'), 'Marginals are currently only supported for the identity transformation')
            assert(dimension >= 1 || dimension <= this.dim);
            
            indices = num2cell((size(this.C) + 1)/2);
            indices{dimension} = ':';
            fd = FourierDistribution.fromComplex(reshape((2 * pi)^(this.dim - 1)*squeeze(this.C(indices{:})), 1, size(this.C, dimension)), 'identity');
        end
        
        function hfd = transformViaFFT(this, desiredTransformation, noOfCoefficients)
            % Calculates transformation of Fourier series via FFT
            % Expects number of complex coefficients (or sum of number of real
            % coefficients)
            assert(strcmp(this.transformation, 'identity'), 'Cannot transform via FFT if already transformed')
            if nargin == 2
                if ismatrix(this.C)&&size(this.C,2)==1 % 1-D case
                    noOfCoefficients=numel(this.C);
                else
                    noOfCoefficients = size(this.C);
                end
            end
            fvals = ifftn(ifftshift(this.C), 'symmetric') * numel(this.C); %Calculate function values via IFFT
            ftmp = HypertoroidalFourierDistribution.fromFunctionValues(fvals, noOfCoefficients, desiredTransformation);
            hfd = this; % To allow for inheritance
            hfd.C = ftmp.C;
            hfd.transformation = ftmp.transformation;
        end
        
        function hfd = transformViaCoefficients(this, desiredTransformation, noOfCoefficients)
            % Calculates transformations using Fourier coefficients
            if nargin == 2, noOfCoefficients = size(this.C);end
            switch desiredTransformation
                case 'identity'
                    hfd = this;
                case 'square'
                    switch this.transformation
                        case 'sqrt'
                            newTrans = 'identity';
                        case 'identity'
                            newTrans = 'square';
                        otherwise
                            newTrans = 'multiple';
                    end
                    hfd = this;
                    % Do not use .multiply as it has different semantic and
                    % would normalize according to sqrt for sqrt
                    % transformation.
                    if all(noOfCoefficients <= size(this.C))
                        hfd.C = convnc(this.C, this.C, 'same');
                    else
                        hfd.C = convnc(this.C, this.C, 'full');
                    end
                    hfd.transformation = newTrans;
                otherwise
                    error('Desired transformation not supported via coefficients');
            end
            hfd = hfd.truncate(noOfCoefficients, true); % Enforce normalization
        end
        
        function hfd = shift(this, shiftAngles)
            % Shift distribution by shiftAngles.
            assert(isequal(numel(shiftAngles), this.dim));
            if ~(any(shiftAngles)) % All angles are zero, do not change anything.
                hfd = this;
                return;
            end
            maxk = (size(this.C) - 1) / 2;
            kRanges = arrayfun(@(currMaxk){-currMaxk:currMaxk}, maxk);
            expFactorMatrices = cell(size(kRanges));
            [expFactorMatrices{:}] = ndgrid(kRanges{:});
            hfd = this;
            % Perform the shift operation (using implicit expansion)
            hfd.C = hfd.C .* exp(-sum(cat(this.dim+1, expFactorMatrices{:}).*reshape(shiftAngles, [ones(1, this.dim), this.dim]), this.dim+1)*1i);
            % The above is equivalent to (use this if no implicit expansion
            % is available)
            %             for i=1:this.dim
            %                 hfd.C=hfd.C.*exp(-expFactorMatrices{i}*shiftAngles(i)*1i);
            %             end
        end
        
        function int = integral(this, l, r)
            % Calculates the integral of the pdf from l to r
            %
            % Parameters:
            %   l (dim x 1 column vector)
            %       left bound of integral in each dimension, default 0
            %   r (dim x 1 column vector)
            %       right bound of integral in each dimension, default 2*pi
            %borders=repmat([0;2*pi],[1,this.dim]);
            %borders(1:length(varargin))=cell2mat(varargin);
            if nargin < 2
                int = 1;
                return;
            end
            if nargin < 3
                r = 2 * pi * ones(this.dim, 1);
            end
            assert(all(size(l) == [this.dim, 1]));
            assert(all(size(r) == [this.dim, 1]));
            
            borders = [l'; r'];
            
            switch this.transformation
                case 'sqrt'
                    hfd = this.transformViaCoefficients('square', 2*size(this.C)-1);
                case 'identity'
                    hfd = this;
                otherwise
                    error('transformation:unrecognizedTransformation', 'Transformation not recognized or unsupported');
            end
            
            maxk = (size(hfd.C) - 1) / 2;
            kRanges = arrayfun(@(currMaxk){-currMaxk:currMaxk}, maxk);
            factorVectors = cell(size(kRanges));
            
            for d = 1:this.dim
                factorVectors{d} = -1i ./ kRanges{d} .* (exp(1i*kRanges{d}*borders(2, d)) - exp(1i*kRanges{d}*borders(1, d)));
                factorVectors{d}(maxk(d) + 1) = borders(2, d) - borders(1, d);
            end
            
            factorMatrices = cell(size(kRanges));
            [factorMatrices{:}] = ndgrid(factorVectors{:});
            allSummands = hfd.C;
            for d = 1:this.dim
                allSummands = allSummands .* factorMatrices{d};
            end
            int = real(sum(allSummands(:)));
        end
        
        function Cov2dimD = covariance2dimD(this)
            % Calculates covariance of [cos(x1), sin(x1), cos(x2), sin(x2), ...,cos(xd), sin(xd)]
            %
            % Returns:
            %   CovdimD (2*dim x 2*dim)
            %       covariance matrix of [cos(x1), sin(x1), cos(x2), sin(x2), ...,cos(xd), sin(xd)]
            switch this.transformation
                case 'sqrt'
                    hfd = this.transformViaCoefficients('square', 5*ones(1, this.dim));
                case 'identity'
                    hfd = this.truncate(5*ones(1, this.dim)); % Resize to convenient size and prevent index out of bounds
                otherwise
                    error('transformation:unrecognizedTransformation', 'Transformation not recognized or unsupported');
            end
            m = hfd.mean2dimD();
            Cov2dimD = zeros(2*this.dim);
            % Formulae are given in the appendix of the paper. The
            % (2*pi)^(this.dim-1) or (2*pi)^(this.dim-2) stem from the
            % notational trick used.
            for d = 1:hfd.dim
                % We know that 3 is center as we truncated to 5.
                indices = num2cell(3*ones(1, this.dim));
                indices(d) = {':'};
                Ccurr = hfd.C(indices{:});
                % Calculate cos(x_d)cos(x_d) terms
                Cov2dimD(2*d-1, 2*d-1) = (2 * pi)^(this.dim - 1) * (pi * real(Ccurr(5)) - 4 * pi * m(2*d-1) * real(Ccurr(4))) + 0.5 * (1 + 2 * m(2*d-1)^2);
                % Calculate sin(x_d)sin(x_d) terms
                Cov2dimD(2*d, 2*d) = (2 * pi)^(this.dim - 1) * (-pi * real(Ccurr(5)) + 4 * pi * imag(Ccurr(4)) * m(2*d)) + 0.5 * (1 + 2 * m(2*d)^2);
                % Calculate cos(x_d)sin(x_d) terms
                Cov2dimD(2*d-1, 2*d) = (2 * pi)^(this.dim - 1) * (-pi * imag(Ccurr(5)) + 2 * pi * imag(Ccurr(4)) * m(2*d-1) - 2 * pi * real(Ccurr(4)) * m(2*d)) + m(2*d-1) * m(2*d);
            end
            for d1 = 1:hfd.dim
                for d2 = d1 + 1:hfd.dim
                    subIndices = repmat({3 * ones(1, 4)}, [1, this.dim]);
                    subIndices(d1) = {[4, 4, 3, 4]};
                    subIndices(d2) = {[2, 3, 4, 4]};
                    Ccurr = hfd.C(sub2ind(size(hfd.C), subIndices{:}));
                    % Calculate cos(x_d1)cos(x_d2) term
                    Cov2dimD(2*d1-1, 2*d2-1) = (2 * pi)^(this.dim - 2) * (2 * pi^2 * real(Ccurr(4)) + 2 * pi^2 * real(Ccurr(1))) ...
                        +(2 * pi)^(this.dim - 1) * (-2 * pi * real(Ccurr(3)) * m(2*d1-1) - 2 * pi * real(Ccurr(2)) * m(2*d2-1)) ...
                        +m(2*d1-1) * m(2*d2-1);
                    % Calculate cos(x_d1)sin(x_d2) term
                    Cov2dimD(2*d1-1, 2*d2) = (2 * pi)^(this.dim - 2) * (2 * pi^2 * imag(-Ccurr(4)+Ccurr(1))) ...
                        +(2 * pi)^(this.dim - 1) * (-2 * pi * real(Ccurr(2)) * m(2*d2) + 2 * pi * imag(Ccurr(3)) * m(2*d1-1)) ...
                        +m(2*d1-1) * m(2*d2);
                    % Calculate sin(x_d1)cos(x_d2) term
                    Cov2dimD(2*d1, 2*d2-1) = (2 * pi)^(this.dim - 2) * (2 * pi^2 * imag(-Ccurr(4)-Ccurr(1))) ...
                        +(2 * pi)^(this.dim - 1) * (-2 * pi * real(Ccurr(3)) * m(2*d1) + 2 * pi * imag(Ccurr(2)) * m(2*d2-1)) ...
                        +m(2*d1) * m(2*d2-1);
                    % Calculate sin(x_d1)sin(x_d2) term
                    Cov2dimD(2*d1, 2*d2) = (2 * pi)^(this.dim - 2) * 2 * pi^2 * real(-Ccurr(4)+Ccurr(1)) ...
                        +(2 * pi)^(this.dim - 1) * (2 * pi * imag(Ccurr(3)) * m(2*d1) + 2 * pi * imag(Ccurr(2)) * m(2*d2)) ...
                        +m(2*d1) * m(2*d2);
                end
            end
            Cov2dimD = Cov2dimD + triu(Cov2dimD, 1)';
        end
        
        function m = trigonometricMoment(this, n)
            % Calculate n-th angular moment from coefficients
            if n == 0
                m = ones(this.dim, 1);
            elseif n < 0
                m = conj(trigonometricMoment(this, -n));
            else
                switch this.transformation
                    case 'sqrt'
                        tfdtmp = this.transformViaCoefficients('square', 2*n*ones(1, this.dim)+1);
                    case 'identity'
                        tfdtmp = this.truncate(2*n*ones(1, this.dim)+1);
                    otherwise
                        error('transformation:unrecognizedTransformation', 'Transformation not recognized or unsupported');
                end
                indicesCell = num2cell((n + 1)*ones(this.dim)-n*eye(this.dim), 2);
                m = (2 * pi)^this.dim * tfdtmp.C(sub2ind(size(tfdtmp.C), indicesCell{:})).';
            end
        end
        
        function p = plot(this, varargin)
            if this.dim ~= 2 % For other dimensions, just fall back to regular plotting
                p = plot@AbstractHypertoroidalDistribution(this, varargin{:});
                return
            end
            noPoints = 102;
            factors = ceil(size(this.C)/(noPoints - 1));
            factorsOdd = factors + ~mod(factors, 2); % Ensure that all factors are odd
            warnStruct = warning('off', 'Truncate:TooFewCoefficients');
            tfdtmp = this.truncate(factorsOdd*(noPoints - 1));
            warning(warnStruct);
            fvals = ifftn(ifftshift(tfdtmp.C), 'symmetric') * numel(tfdtmp.C);
            [alpha, beta] = meshgrid(linspace(0, 2*pi, noPoints));
            p = surf(alpha, beta, fvals([1:factorsOdd(1):end, 1], [1:factorsOdd(2):end, 1]).');
        end
    end
    
    methods(Static)
        function hfd = fromFunction(fun, noOfCoefficients, desiredTransformation)
            % Creates Fourier distribution from function
            % Function must be able to take vector arguments
            % Dimension of noOfCoefficients has to be in accordance to
            % dimensionality of the function to approximate.
            if nargin == 2, desiredTransformation = 'sqrt';end
            assert((length(noOfCoefficients) == nargin(fun)) || (nargin(fun) == -1), 'noOfCoefficient has to match dimensionality (in form of number of input arguments) of the function.');
            gridIndividualAxis = arrayfun(@(currNo){0:2 * pi / currNo:2 * pi - 2 * pi / currNo}, noOfCoefficients);
            gridCell = cell(1, length(noOfCoefficients));
            [gridCell{:}] = ndgrid(gridIndividualAxis{:});
            fvals = fun(gridCell{:}); %Assume functions are vectorized!
            assert(numel(fvals) == prod(noOfCoefficients), 'Size of output of pdf is incorrect. Please ensure that pdf returns only one scalar per dim-tupel.');
            hfd = HypertoroidalFourierDistribution.fromFunctionValues(fvals, noOfCoefficients, desiredTransformation);
        end
        
        function hfd = fromFunctionValues(fvals, noOfCoefficients, desiredTransformation)
            % Creates Fourier distribution from function values
            % Assumes fvals are not yet transformed, use custom if they already
            % are transformed
            assert(all(noOfCoefficients > 1) && all(mod(noOfCoefficients-1, 2) == 0),...
                'Invalid number of coefficients, numbers for all dimensions have to be odd and greater than 1.');
            % Cannot directly compare the size of fvals with noOfCoefficients because we
            % allow truncation afterward. But we ensure there is no 1 x n
            % matrix by accident, since it has to be n x 1.
            assert(all(size(fvals)>1)||ismatrix(fvals)&&size(fvals,2)==1,...% Later condition ensures [n,1] matrices work
                'Some dimension has only one entry along it. Fix this.');
            % Ensure dimensions match
            assert(numel(noOfCoefficients)==ndims(fvals) || numel(noOfCoefficients)==1);
            switch desiredTransformation
                case 'sqrt'
                    fvals = sqrt(fvals);
                case 'log'
                    fvals = log(fvals);
                case 'identity' %keep them unchanged
                case 'custom' %already transformed
                otherwise
                    error('fromFunctionValues:unrecognizedTranformation', 'Transformation not recognized or unsupported by transformation via FFT');
            end
            hfd = HypertoroidalFourierDistribution(fftshift(fftn(fvals)/numel(fvals)), desiredTransformation);
            hfd = hfd.truncate(noOfCoefficients);
        end
        
        function hfd = fromDistribution(distribution, noOfCoefficients, desiredTransformation)
            % Creates Fourier distribution from a different distribution
            % Always uses FFT as no formulae for pupular distributions
            % available
            if numel(noOfCoefficients) == 1
                noOfCoefficients = ones(1, distribution.dim) * noOfCoefficients;
            end
            assert(all(mod(noOfCoefficients,2)==1),'noOfCoefficients must be odd');
            assert(isa(distribution, 'AbstractToroidalDistribution') && (size(noOfCoefficients, 2) == 2) || ...
                isa(distribution, 'AbstractHypertoroidalDistribution') && (size(noOfCoefficients, 2) == distribution.dim), ...
                'fromDistribution:invalidObject', 'First argument has to be a (hyper)toroidal distribution with appropriate dimensionality.');
            if nargin == 2
                desiredTransformation = 'sqrt';
            end
            if isa(distribution, 'AbstractCircularDistribution') % Convert via FourierDistribution for circular distributions
                fd = FourierDistribution.fromDistribution(distribution, noOfCoefficients, desiredTransformation);
                hfd = HypertoroidalFourierDistribution(fd.c.', desiredTransformation);
            elseif isa(distribution, 'HypertoroidalUniformDistribution')
                C = zeros([noOfCoefficients, 1]); % Extra 1 is for 1D case
                indexc00 = num2cell((size(C) + 1)/2);
                switch desiredTransformation
                    case 'sqrt'
                        C(indexc00{:}) = 1 / sqrt((2 * pi)^distribution.dim);
                    case 'identity'
                        C(indexc00{:}) = 1 / (2 * pi)^distribution.dim;
                    otherwise
                        error('Transformation not recognized or unsupported');
                end
                hfd = HypertoroidalFourierDistribution(C, desiredTransformation);
            elseif isa(distribution, 'HypertoroidalWDDistribution')
                coeffGenerator = @(k) sum(exp(-1i*(k*distribution.d)).*distribution.w);
                gridCell = cell(1, distribution.dim);
                maxDeg = (noOfCoefficients - 1) / 2;
                [gridCell{:}] = ndgrid(-maxDeg:maxDeg);
                C = arrayfun(@(varargin)coeffGenerator(cat(2, varargin{:})), gridCell{:});
                hfdtmp = HypertoroidalFourierDistribution(1/(2*pi)^distribution.dim * C, 'identity');
                switch desiredTransformation
                    case 'identity'
                        hfd = hfdtmp;
                    case 'sqrt'
                        warning('Conversion:WDMatching', 'Approximating WD by matching moments');
                        hfd = hfdtmp.transformViaFFT('sqrt', noOfCoefficients);
                    case 'log'
                        warning('Conversion:WDMatching', 'Approximating WD by matching moments');
                        hfd = hfdtmp.transformViaFFT('log', noOfCoefficients);
                    otherwise
                        error('Transformation not recognized or unsupported');
                end
            elseif strcmp(desiredTransformation, 'identity') && isa(distribution, 'HypertoroidalWNDistribution')
                % Could also handle 1-D explicitly
                if (distribution.dim == 2) % Solution for 2D that is more efficient in both computation and memory
                    maxk = (noOfCoefficients - 1) / 2;
                    C = 1 / (2 * pi)^2 * (exp(-1i*(-maxk(1):maxk(1))'*distribution.mu(1)) * exp(-1i*(-maxk(2):maxk(2))*distribution.mu(2)) ...
                        ./ (exp(0.5*(-maxk(1):maxk(1))'.^2*distribution.C(1, 1)) * exp(0.5*(-maxk(2):maxk(2)).^2*distribution.C(2, 2))) ...
                        ./ exp((-maxk(1):maxk(1))'*(-maxk(2):maxk(2))*distribution.C(1, 2)));
                else
                    kRanges = arrayfun(@(currMaxk){-currMaxk:currMaxk}, (noOfCoefficients - 1)/2);
                    individualIndexMatrices = cell(size(kRanges))';
                    [individualIndexMatrices{:}] = ndgrid(kRanges{:});
                    indAsVec = cell2mat(cellfun(@(mat){mat(:)'}, individualIndexMatrices));
                    C = reshape(exp(-1i*indAsVec'*distribution.mu- ...
                        0.5*pdist2(indAsVec', zeros(1, distribution.dim), 'mahalanobis', inv(distribution.C)).^2) ...
                        /(2 * pi)^distribution.dim, size(individualIndexMatrices{1}));
                end
                if ~any(isnan(C)) % The exponential function can lead to Inf/NaN values. If this occurs, use FFT instead
                    hfd = HypertoroidalFourierDistribution(C, desiredTransformation);
                else
                    pdfFun = @(varargin)reshape(distribution.pdf(cell2mat(cellfun(@(currCell){currCell(:)}, varargin))'), size(varargin{1}));
                    hfd = HypertoroidalFourierDistribution.fromFunction(pdfFun, noOfCoefficients, desiredTransformation);
                end
            else
                pdfFun = @(varargin)reshape(distribution.pdf(cell2mat(cellfun(@(currCell){currCell(:)}, varargin))'), size(varargin{1}));
                hfd = HypertoroidalFourierDistribution.fromFunction(pdfFun, noOfCoefficients, desiredTransformation);
            end
        end
        
        function hfd = fromSamples(samples, noOfCoefficients, desiredTransformation)
            if nargin == 2
                desiredTransformation = 'sqrt';
            end
            wd = HypertoroidalWDDistribution(samples);
            hfd = HypertoroidalFourierDistribution.fromDistribution(wd,...
                noOfCoefficients, desiredTransformation);
        end
    end
end