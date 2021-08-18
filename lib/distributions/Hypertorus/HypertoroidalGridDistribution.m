classdef HypertoroidalGridDistribution < AbstractHypertoroidalDistribution & AbstractGridDistribution
	properties
        noOfGridPoints double = []
	end
    methods
        function this = HypertoroidalGridDistribution(grid_, gridValues_, enforcePdfNonnegative_, dim)
            % Constructor
            arguments % Use to set default value
                grid_ double {mustBeGreaterThanOrEqual(grid_,0)}
                gridValues_ (:,1) double {mustBeNonnegative}
                enforcePdfNonnegative_ (1,1) logical = true
                dim (1,1) double {mustBePositive,mustBeInteger} = max(size(grid_,1),1); % Can enforce dimension when no grid is given
            end
            % Verify grid size (skip if empty)
            assert(size(grid_,1)==dim||isempty(grid_));
            assert(size(grid_,2)==size(gridValues_,1)||isempty(grid_));
            if size(grid_,1)>size(grid_,2)
                warning('Dimension is higher than number of grid points. Verify that this is really intended.');
            end
            this.dim = dim;
            this.grid = grid_;
            this.gridValues = gridValues_;
            this.noOfGridPoints = numel(gridValues_);
            this.enforcePdfNonnegative = enforcePdfNonnegative_;
            % Check if normalized. If not: Normalize!
            this = this.normalize();
        end
        
        function int = integral(this, l, r)
            arguments
                this (1,1) HypertoroidalGridDistribution
                l double = []
                r double = []
            end
            int = integral@AbstractGridDistribution(this, l, r);
        end
        
        function f = normalize(this, opt)
            arguments
                this (1,1) HypertoroidalGridDistribution
                opt.tol (1,1) double = 1e-2
                opt.warnUnnorm (1,1) logical = true
            end
            f = normalize@AbstractGridDistribution(this,tol=opt.tol,warnUnnorm=opt.warnUnnorm);
        end
                
        function hgd = multiply(this, other)
            arguments
                this HypertoroidalGridDistribution
                other HypertoroidalGridDistribution
            end
            assert(isequal(this.grid,other.grid), 'Multiply:IncompatibleGrid','Can only multiply for equal grids.');
            hgd = multiply@AbstractGridDistribution(this, other);
        end   
        
        function p = pdf(this, xa)
            arguments
                this (1,1) HypertoroidalGridDistribution
                xa (:,:) double
            end
            assert(strcmp(this.gridType,'CartesianProd'),'pdf is not defined. Can only interpolate for certain grid types.');
            % Using Fourier-based interpolation
            warning('PDF:UseInterpolated', 'pdf is not defined. Using interpolation with Fourier series.')
            if this.enforcePdfNonnegative
                transformation = 'sqrt';
            else
                transformation = 'identity';
            end
            if isempty(this.noOfGridPoints)
                warning('Cannot interpolate if number of grid points are not specified. Assuming equidistant grid');
                sizes = repmat(size(this.grid,2)^(1/this.dim),[1,this.dim]);
                assert(max(abs(sizes-round(sizes)))<1e-10);
                sizes=round(sizes);
                this.noOfGridPoints=sizes;
            end
            
            fd = HypertoroidalFourierDistribution.fromFunctionValues(...
                reshape(this.gridValues,[this.noOfGridPoints,1]),... % Append 1 for 1-D case
                this.noOfGridPoints+double(mod(this.noOfGridPoints,2)==0), transformation);
            p = fd.pdf(xa);
        end
        
        function grid = getGrid(this)
            if ~isempty(this.grid)
                grid = this.grid;
            elseif strcmp(this.gridType,'CartesianProd')
                warning('Grid:GenerateDuringRunTime','Generating grid anew on call to .getGrid(). If you require the grid frequently, store it in the class.')
                grid = HypertoroidalGridDistribution.generateCartesianProductGrid(this.noOfGridPoints);
            else
                error('Grid:UnknownGrid','Grid was not provided and is thus unavailable');
            end
        end
        
        function p = pdfUnnormalized(this, xs)
            % Use this when wanting to evaluate likelihoods without any
            % normalization
            arguments
                this (1,1) HypertoroidalGridDistribution
                xs (:,:) double
            end
            assert(strcmp(this.gridType,'CartesianProd'),'pdf is not defined. Can only interpolate for certain grid types.');
            warningStatus = warning('off','Normalization:notNormalized'); % We know it is unnormalized and that is okay
            p = this.integral()*this.pdf(xs);
            warning(warningStatus);
        end
        
        function gd = shift(this, angles)
            arguments
                this (1,1) HypertoroidalGridDistribution
                angles (:,1) double
            end
            if norm(angles)==0
                gd=this;
                return
            end
            assert(strcmp(this.gridType,'CartesianProd'));
            % Go via HypertoroidalFourierDistribution but do not use
            % constructor to avoid normalization. Therefore initialize with
            % some uniform distribution and then replace coefficient matrix
            coeffMatTmp = zeros(repmat(3,[1,this.dim]));
            coeffMatTmp(1) = 1 / (2 * pi)^this.dim;
            hfd = HypertoroidalFourierDistribution(fftshift(coeffMatTmp),'identity');
            hfd.C =  fftshift(fftn(reshape(this.gridValues,this.noOfGridPoints)));
            hfdShifted = hfd.shift(angles);
            gd = this;
            gd.gridValues = reshape(ifftn(ifftshift(hfdShifted.C), 'symmetric'),size(this.gridValues));
        end
        
        function p = valueOfClosest(this, xa)
            arguments
                this (1,1) HypertoroidalGridDistribution
                xa (:,:) double
            end
            % Could make more efficient with dsearchn, but have to consider
            % all possible wrap arounds
            p = NaN(1,size(xa,2));
            for i=1:size(xa,2)
                dists = sum(min((this.grid-xa(:,i)).^2,(2*pi-(this.grid-xa(:,i))).^2));
                [~,minInd] = min(dists);
                p(i) = this.gridValues(minInd);
            end
        end
        
        function gd = convolve(this, other)
            arguments
                this (1,1) HypertoroidalGridDistribution
                other (1,1) HypertoroidalGridDistribution
            end
            assert(strcmp(this.gridType,'CartesianProd')&&strcmp(other.gridType,'CartesianProd'));
            assert(isequal(this.noOfGridPoints,other.noOfGridPoints));
            assert(isequal(this.grid,other.grid)); % Should be equal due to the above two lines, but we validate it anyways
            
            thisTensor = reshape(this.gridValues,this.noOfGridPoints);
            otherTensor = reshape(other.gridValues,other.noOfGridPoints);
            
            resTensor = (2*pi)^this.dim / size(this.gridValues,1) * ifftn(fftn(thisTensor).*fftn(otherTensor));
            gd = this;
            gd.gridValues = reshape(resTensor,size(this.gridValues));
        end
        
        function h = plot(this,varargin)
            arguments
                this (1,1) HypertoroidalGridDistribution
            end
            arguments (Repeating)
                varargin
            end
            % Initially equally weighted and then overwrite to prevent a
            % normalization in the constructor the WD distribution
            hdd = HypertoroidalWDDistribution(this.grid, 1/numel(this.gridValues)*ones(1,numel(this.gridValues)));
            hdd.w = this.gridValues';
            h = hdd.plot(varargin{:});
        end
        
        function h = plotInterpolated(this,varargin)
            arguments
                this (1,1) HypertoroidalGridDistribution
            end
            arguments (Repeating)
                varargin
            end
            if this.dim>2  
                error('Can currently only plot for T1 and T2 torus.')
            end
            if this.enforcePdfNonnegative
                transformation = 'sqrt';
            else
                transformation = 'identity';
            end
            sizes = repmat(numel(this.gridValues)^(1/this.dim),[1,this.dim]);
            fd = HypertoroidalFourierDistribution.fromFunctionValues(reshape(this.gridValues,[sizes,1]),...
                sizes, transformation);
            h = fd.plot(varargin{:});
        end
        
        function m = trigonometricMoment(this, n)
            hwd = HypertoroidalWDDistribution(this.grid,this.gridValues'/sum(this.gridValues));
            m = hwd.trigonometricMoment(n);
        end
        
        function hgd = sliceAt(this, dims, val, useFFTN)
            arguments
                this (1,1) HypertoroidalGridDistribution
                dims (1,:) double {mustBeInteger,mustBePositive,mustBeNonempty}
                val (:,1) double {mustBeNonempty}
                useFFTN (1,1) logical = false % transform all dimensions instead of only the ones in question
            end
            assert(numel(dims)==numel(val));
            assert(strcmp(this.gridType,'CartesianProd'),'This operation is only supported for grids generated based on a Cartesian product.');
            assert(all(dims<=this.dim),'Cannot perform this operation for a dimension that is higher than the dimensionality of the distribution.');
            
            
            fvalsOnGrid = reshape(this.gridValues,this.noOfGridPoints);
            if isequal(val,zeros(size(val)))
                % No shifting whatsoever required. We do not need to go to
                % the Fourier domain
                gridShiftedFull = fvalsOnGrid;
            elseif useFFTN
                % Use HypertoroidalFourierDistribution, which uses FFTN
                shiftVec = zeros(this.dim,1);
                shiftVec(dims) = -val;
                hfdSqrt = HypertoroidalFourierDistribution.fromFunctionValues(fvalsOnGrid,this.noOfGridPoints,'sqrt');
                hfdSqrtShifted = hfdSqrt.shift(shiftVec);
                gridShiftedFull = hfdSqrtShifted.pdfOnGrid();
            else
                hybridVals = sqrt(fvalsOnGrid);
                for d=dims
                    hybridVals = fft(hybridVals,[],d);
                end
                shiftFactors = 1;
                % This loop will expand it into a tensor with full size for
                % every dimension that is shifted and 1 for all others.
                % Thus, it may be signficiantly smaller than the full
                % matrix.
                for i=1:numel(dims)
                    targetSize = ones(1,this.dim);
                    targetSize(dims(i)) = size(hybridVals,dims(i));
                    kd = reshape([0:(size(hybridVals,dims(i))-1)/2,-(size(hybridVals,dims(i))-1)/2:-1],targetSize);
                    % *No* minus because there is a minus in the formula
                    % and we need to use - to shift it to zero.
                    shiftFactors = shiftFactors.*exp(1i*kd*val(i));
                end
                gridShiftedFull = hybridVals.*shiftFactors;
                for d=dims
                    if d==dims(end)
                        gridShiftedFull = ifft(gridShiftedFull,[],d,'symmetric');
                    else
                        gridShiftedFull = ifft(gridShiftedFull,[],d);
                    end
                end
                assert(all(imag(gridShiftedFull)<0.001,1:this.dim));
                gridShiftedFull = real(gridShiftedFull).^2;
            end
            indArray = repmat({':'},1,this.dim);
            indArray(dims) = {1};
            sizes = this.noOfGridPoints;
            sizes(dims) = [];
            gridValsCond = reshape(gridShiftedFull(indArray{:}),[sizes,1]); % Including 1 to have at least 2 dims

            hgd = this;
            hgd.gridValues = gridValsCond(:);
            hgd.noOfGridPoints(dims) = [];
            hgd.dim = this.dim-numel(dims);
            if ~isempty(this.grid) % Generate new grid (only if it was not empty before)
                hgd.grid = HypertoroidalGridDistribution.generateCartesianProductGrid(hgd.noOfGridPoints);
            end
        end
        
        function hgd = conditionOn(this, dims, val, useFFTN)
            arguments
                this (1,1) HypertoroidalGridDistribution
                dims (1,:) double {mustBeInteger,mustBePositive,mustBeNonempty}
                val (:,1) double {mustBeNonempty}
                useFFTN (1,1) logical = false
            end
            assert(all(dims<=this.dim),'Cannot perform this operation for a dimension that is higher than the dimensionality of the distribution.');
            assert(numel(dims)==numel(val));
            hgd = this.sliceAt(dims, val, useFFTN);
            hgd = hgd.normalize(warnUnnorm=false);
        end
        
        function hgdMarginalized = marginalizeOut(this, dimensions)
            arguments
                this (1,1) HypertoroidalGridDistribution
                dimensions (1,:) double {mustBeInteger,mustBePositive}
            end
            assert(all(dimensions<=this.dim),'Cannot perform this operation for a dimension that is higher than the dimensionality of the distribution.');
            assert(strcmp(this.gridType,'CartesianProd'),'This operation is only supported for grids generated based on a Cartesian product.');
            
            hgdMarginalized = this;
            hgdMarginalized.dim = this.dim - numel(dimensions);
            hgdMarginalized.noOfGridPoints(dimensions) = [];
            
            newGridValuesTensor = (2*pi)^(numel(dimensions))*mean(reshape(this.gridValues,this.noOfGridPoints),dimensions);
            hgdMarginalized.gridValues = newGridValuesTensor(:);
            if ~isempty(this.grid)
                hgdMarginalized.grid = HypertoroidalGridDistribution.generateCartesianProductGrid(hgdMarginalized.noOfGridPoints);
            end
        end
        
        function ghdLikelihood = likelihood(this, dims, val, useFFTN)
            arguments
                this (1,1) HypertoroidalGridDistribution
                dims (1,:) double {mustBeInteger,mustBePositive}
                val (:,1) double
                useFFTN (1,1) logical = false
            end
            assert(numel(dims)==numel(val));
            ghdSliced = this.sliceAt(dims, val, useFFTN);
            ghdMarginalized = this.marginalizeOut(dims);
            ghdLikelihood = ghdSliced;
            ghdLikelihood.gridValues = ghdSliced.gridValues./ghdMarginalized.gridValues;
        end
    end
    
    methods (Static)
        function hgd = fromDistribution(dist, noOfGridPoints, gridType)
            arguments
                dist AbstractHypertoroidalDistribution
                noOfGridPoints {mustBeInteger,mustBePositive}
                gridType char = 'CartesianProd'
            end
            if strcmp(gridType,'CartesianProd')&&(isa(dist,'HypertoroidalFourierDistribution')&&...
                    (isequal(size(dist.C),noOfGridPoints)||isscalar(noOfGridPoints)&&noOfGridPoints==size(dist.C,1))) % Latter is for 1-D case
                gridValues = dist.pdfOnGrid();
                grid = HypertoroidalGridDistribution.generateCartesianProductGrid(noOfGridPoints);
                hgd = HypertoroidalGridDistribution(grid, gridValues(:)');
                hgd.gridType = gridType;
                hgd.noOfGridPoints = noOfGridPoints;
            else
                hgd = HypertoroidalGridDistribution.fromFunction(@(x)dist.pdf(x), noOfGridPoints, dist.dim, gridType);
            end
        end
        function sgd = fromFunction(fun, noOfGridPoints, dim, gridType)
            arguments
                fun function_handle
                noOfGridPoints {mustBeInteger,mustBePositive}
                dim {mustBeInteger,mustBeGreaterThanOrEqual(dim,1)}
                gridType char = 'CartesianProd'
            end
            assert(nargin(fun)==1);
            switch gridType
                case 'CartesianProd'
                    if numel(noOfGridPoints)==1
                        noOfGridPoints = repmat(noOfGridPoints,1,dim);% assuming equal resolution
                    else
                        assert(numel(noOfGridPoints)==dim);
                    end
                    grid = HypertoroidalGridDistribution.generateCartesianProductGrid(noOfGridPoints);
                otherwise
                    error('Grid scheme not recognized');
            end
            gridValues = fun(grid)';
            sgd = HypertoroidalGridDistribution(grid, gridValues);
            sgd.gridType = gridType;
            sgd.noOfGridPoints = noOfGridPoints;
        end
        function grid = generateCartesianProductGrid(noOfGridPoints)
            arguments
                noOfGridPoints (1,:) double {mustBePositive,mustBeInteger,mustBeNonempty}
            end
            gridIndividualAxis = arrayfun(@(currNo){0:2 * pi / currNo:2 * pi - 2 * pi / currNo}, noOfGridPoints);
            gridCell = cell(1, length(noOfGridPoints));
            [gridCell{:}] = ndgrid(gridIndividualAxis{:});

            dim = numel(noOfGridPoints);
            gridMatConcat=cat(dim+1,gridCell{:});
            gridMatForIteration=permute(gridMatConcat,[dim+1,1:dim]);
            grid=reshape(gridMatForIteration,[dim,prod(noOfGridPoints)]);
        end
    end
end

