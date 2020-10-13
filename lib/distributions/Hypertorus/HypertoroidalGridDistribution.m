classdef HypertoroidalGridDistribution < AbstractHypertoroidalDistribution & AbstractGridDistribution
	properties
        grid double {mustBeGreaterThanOrEqual(grid,0)}
        gridType char = 'unknown'
        noOfGridPoints double = []
	end
    methods
        function this = HypertoroidalGridDistribution(grid_, gridValues_, enforcePdfNonnegative_)
            % Constructor
            arguments % Use to set default value
                grid_ double {mustBeGreaterThanOrEqual(grid_,0)}
                gridValues_ (:,1) double {mustBeNonnegative}
                enforcePdfNonnegative_ logical = true
            end
            assert(size(grid_,2)==size(gridValues_,1));
            if size(grid_,1)>size(grid_,2)
                warning('Dimension is higher than number of grid points. Verify that is really intended.');
            end
            this.dim = size(grid_,1);
            this.grid = grid_;
            this.gridValues = gridValues_;
            this.enforcePdfNonnegative = enforcePdfNonnegative_;
            % Check if normalized. If not: Normalize!
            this = this.normalize;
        end      
        
        function f = normalize(this)
            tol = 1e-2;
            f = normalize@AbstractGridDistribution(this,tol);
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
                this.noOfGridPoints, transformation);
            p = fd.pdf(xa);
        end
        
        function gd = shift(this, angle)
            if norm(angle)==0
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
            hfdShifted = hfd.shift(angle);
            gd = this;
            gd.gridValues = reshape(ifftn(ifftshift(hfdShifted.C), 'symmetric'),size(this.gridValues));
        end
        
        function p = valueOfClosest(this, xa)
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
            assert(strcmp(this.gridType,'CartesianProd')&&strcmp(other.gridType,'CartesianProd'));
            assert(isequal(this.noOfGridPoints,other.noOfGridPoints));
            assert(isequal(this.grid,other.grid)); % Should be equal due to the above two lines, but we validate it anyways
            
            thisTensor = reshape(this.gridValues,this.noOfGridPoints);
            otherTensor = reshape(other.gridValues,other.noOfGridPoints);
            
            resTensor = (2*pi)^this.dim / size(this.gridValues,1) * ifftn(fftn(thisTensor).*fftn(otherTensor));
            gd = this;
            gd.gridValues = reshape(resTensor,size(this.gridValues));
        end
        
        function h = plot(this)
            hdd = HypertoroidalWDDistribution(this.grid, this.gridValues');
            h = hdd.plot;
        end
        
        function h = plotInterpolated(this,varargin)
            if this.dim>2  
                error('Can currently only plot for T1 and T2 torus.')
            end
            if this.enforcePdfNonnegative
                transformation = 'sqrt';
            else
                transformation = 'identity';
            end
            sizes = repmat(size(this.grid,2)^(1/this.dim),[1,this.dim]);
            if this.dim==1
                reshapeTo = [sizes,1];
            else
                reshapeTo = sizes;
            end
            fd = HypertoroidalFourierDistribution.fromFunctionValues(reshape(this.gridValues,reshapeTo),...
                sizes, transformation);
            h = fd.plot(varargin{:});
        end
        
        function m = trigonometricMoment(this, n)
            hwd = HypertoroidalWDDistribution(this.grid,this.gridValues'/sum(this.gridValues));
            m = hwd.trigonometricMoment(n);
        end
        
    end
    
    methods (Static)
        function sgd = fromDistribution(dist, noOfGridPoints, gridType)
            arguments
                dist AbstractHypertoroidalDistribution
                noOfGridPoints {mustBeInteger,mustBePositive}
                gridType char = 'CartesianProd'
            end
            sgd = HypertoroidalGridDistribution.fromFunction(@(x)dist.pdf(x), noOfGridPoints, dist.dim, gridType);
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
                    gridIndividualAxis = arrayfun(@(currNo){0:2 * pi / currNo:2 * pi - 2 * pi / currNo}, noOfGridPoints);
                    gridCell = cell(1, length(noOfGridPoints));
                    [gridCell{:}] = ndgrid(gridIndividualAxis{:});
                    
                    gridMatConcat=cat(dim+1,gridCell{:});
                    gridMatForIteration=permute(gridMatConcat,[dim+1,1:dim]);
                    grid=reshape(gridMatForIteration,[dim,prod(noOfGridPoints)]);
                otherwise
                    error('Grid scheme not recognized');
            end
            gridValues = fun(grid)';
            sgd = HypertoroidalGridDistribution(grid, gridValues);
            sgd.gridType = gridType;
            sgd.noOfGridPoints = noOfGridPoints;
        end
    end
end

