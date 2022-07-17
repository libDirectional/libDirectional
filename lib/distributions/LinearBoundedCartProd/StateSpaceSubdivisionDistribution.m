classdef StateSpaceSubdivisionDistribution < AbstractLinPeriodicDistribution
	properties
        gd
        linearDistributions
	end
    methods
        function this = StateSpaceSubdivisionDistribution(gd_, linearDistributions_)
            % Rao Blackwallized Distributions with aribtrary distributions
            % on the linear part. Use StateSpaceSubdivisionGaussian for Gaussian
            % distributions.
            arguments
                gd_ (1,1) AbstractGridDistribution
                linearDistributions_ (:,1) AbstractLinearDistribution
            end
            assert(isequal(numel(linearDistributions_),numel(gd_.gridValues)),...
                'The number of linear distributions and grid points must be equal.');
            this.linD = linearDistributions_(1).dim;
            this.boundD = gd_.dim; 
            this.dim = this.linD + this.boundD;
            this.gd = gd_;
            this.linearDistributions = linearDistributions_;
            % Normalize
            this = this.normalize;
        end      
        
        function s = sample(this, n)
            arguments
                this (1,1) StateSpaceSubdivisionDistribution
                n (1,1) {mustBePositive,mustBeInteger}
            end
            samplesBounded = this.gd.sample(n);
            % Draw from linear distriution of nearest neighbor
            [~,indices] = this.gd.getClosestPoint(samplesBounded);
            relevantGridIds = unique(indices);
            samplesLinear = NaN(this.linD, n);
            for currGridInd = relevantGridIds
                currXaIndices = indices==currGridInd;
                samplesLinear(:,currXaIndices) = this.linearDistributions(currGridInd).sample(sum(currXaIndices));
            end
            s = [samplesBounded; samplesLinear];
        end

        function d = marginalizeLinear(this)
            d = this.gd;
        end

        function d = marginalizePeriodic(this)
            d = LinearMixture(mat2cell(this.linearDistributions,ones(1,numel(this.linearDistributions)),1), this.gd.gridValues/sum(this.gd.gridValues)); %#ok<MMTC> 
        end
        
        function p = pdf(this, xa, linDistInterpolationMethod, boundDistInterpolationmethod)
            arguments
                this (1,1) StateSpaceSubdivisionDistribution
                xa (:,:) double {mustBeNonempty,mustBeNonNan}
                linDistInterpolationMethod ...
                    char {mustBeMember(linDistInterpolationMethod,{'','nearestNeighbor','mixture','convexCombinationGaussian'})} = ''
                boundDistInterpolationmethod char = 'gridDefault' % This is usually significantly better
            end
            assert(size(xa,1)==this.dim, 'Dimension of xa does not match the density''s dimension.');
            % Set reasonable defaults for the interpolation methods
            if isempty(linDistInterpolationMethod)&&isa(this.gd,'HypertoroidalGridDistribution')&&this.boundD==1
                linDistInterpolationMethod = 'mixture'; % May not be optimal to use this as default
            elseif isempty(linDistInterpolationMethod)
                linDistInterpolationMethod = 'nearestNeighbor';
            end
            xaBound = xa(1:this.boundD,:);
            xaLin = xa(this.boundD+1:end,:);
            switch linDistInterpolationMethod
                case 'nearestNeighbor'
                    % Find nearest neighbor
                    % This automatically uses modulo arithmetics for
                    % FIGDistribution
                    [~, indices] = this.gd.getClosestPoint(xaBound);
                    % Under the assumption that the evaluation is less expensive if
                    % you input all, we iterate over all grid points. However, do
                    % not do this if only one grid point is requested
                    if size(xa,2)==1
                        flincondbound = this.linearDistributions(indices).pdf(xaLin(:,1));
                    else
                        relevantGridIds = unique(indices);
                        flincondbound = NaN(1,size(xa,2));
                        for currGridInd = relevantGridIds
                            currXaIndices = indices==currGridInd;
                            flincondbound(currXaIndices) = this.linearDistributions(currGridInd).pdf(xaLin(:,currXaIndices));
                        end
                        assert(~anynan(flincondbound));
                    end
                case 'mixture'
                    % Use a mixture of the densities according to the
                    % weights to the closest points
                    assert(isa(this.gd,'HypertoroidalGridDistribution')&&this.boundD==1,'Currently only supporting circle.')
                    xaBound = mod(xaBound,2*pi);
                    posOnGridScaled = xaBound/(2*pi/this.gd.noOfGridPoints);
                    lowerIndices = floor(posOnGridScaled)+1;
                    % lowerIndex + 1 (except in the case that we leave the
                    % area)
                    higherIndices = lowerIndices+1 - (lowerIndices==this.gd.noOfGridPoints)*this.gd.noOfGridPoints;
                    % We cannot trivially it for upper and lower index
                    % simultaenous as we would get a 1x2n or 2xn vector
                    % coding for the quality of the currGridInd and the
                    % full list of grid points. This would not 
                    % for 2xn for higher and 
                    lowerAndHigherIndices = [lowerIndices,higherIndices]; % (1x2n) to handle it optimally below
                    flincondboundLowerAndUpper = NaN(1,2*size(xa,2));
                    xaLin = [xaLin,xaLin]; % Duplicate that to use the boolean indexing for the 1x2n vector
                    for currGridInd = unique(lowerAndHigherIndices)
                        currXaIndices = lowerAndHigherIndices==currGridInd;
                        flincondboundLowerAndUpper(currXaIndices) = ...
                            this.linearDistributions(currGridInd).pdf(xaLin(:,currXaIndices));
                    end
                    % Do the mixture combination with the pdf value
                    % proportional to the density
                    weightings = [-(posOnGridScaled-lowerIndices);(posOnGridScaled-lowerIndices)+1];
                    flincondbound = sum(weightings.*reshape(flincondboundLowerAndUpper,[],2)');
                case 'convexCombinationGaussian'
                    % Use a convex combination of the parameters of the
                    % Gaussians
                    assert(isa(this.gd,'HypertoroidalGridDistribution')&&this.boundD==1,'Currently only supporting circle.')
                    assert(isa(this,'StateSpaceSubdivisionGaussianDistribution')||all(arrayfun(@(d)isa(d,'GaussianDistribution'),this.linearDistributions)),'This combination method is only supported for StateSpaceSubdivisionGaussianDistributions and derived classes.');
                    xaBound = mod(xaBound,2*pi);
                    posOnGridScaled = xaBound/(2*pi/this.gd.noOfGridPoints);
                    lowerIndices = floor(posOnGridScaled)+1;
                    higherIndices = lowerIndices+1 - (lowerIndices==this.gd.noOfGridPoints)*this.gd.noOfGridPoints;
                    weightings = [-(posOnGridScaled-lowerIndices);(posOnGridScaled-lowerIndices)+1];
                    flincondbound = NaN(size(lowerIndices));
                    for i=1:size(xa,2)
                        currMu = weightings(1,i)*this.linearDistributions(lowerIndices(i)).mu...
                            + weightings(2,i)*this.linearDistributions(higherIndices(i)).mu;
                        currC = weightings(1,i)*this.linearDistributions(lowerIndices(i)).C...
                            + weightings(2,i)*this.linearDistributions(higherIndices(i)).C;
                        flincondbound(i) = GaussianDistribution(currMu,currC).pdf(xaLin(:,i));
                    end
                    
                otherwise
                    error('Interpolation of linear density unsupported');
            end
            % Use the interpolation method of the grid distribution here
            % for the periodic part (!)
            switch boundDistInterpolationmethod
                case 'nearestNeighbor'
                    if ~strcmp(linDistInterpolationMethod,'nearestNeighbor')
                        indices = mod(round(posOnGridScaled),this.gd.noOfGridPoints)+1;
                    end
                    fBound = this.gd.gridValues(indices)';
                case 'gridDefault'
                    fBound = this.gd.pdf(xaBound);
                otherwise
                    error('Interpolation of bounded density unsupported');
            end
            p = fBound.*flincondbound;
        end
        
        function dist = normalize(this)
            arguments
                this (1,1) StateSpaceSubdivisionDistribution
            end
            dist = this;
            % Normalize only grid part. Conditional part (linearDistributions) are
            % always normalized
            dist.gd = this.gd.normalize;
        end
        
        function res = multiply(~, ~)  %#ok<STOUT> keep output to ensure the call does not fail with "too many output arguments" instead of our (better) message below
            error('Not supported for arbitrary densities. Use StateSpaceSubdivisionGaussianDistribution.');
        end   
        
        function m = mode(this)
            % We have to check all the modes from the linear distributions
            % and multiply them by the weights of the grid points. The
            % highest value is then the max (not interpolated).
            arguments
                this (1,1) StateSpaceSubdivisionDistribution
            end
            pdfsAtGridPoints = NaN(size(this.gd.gridValues));
            
            mLinModesPos = NaN(size(this.gd.gridValues));
            for i=1:size(this.gd.gridValues)
                mLinModesPos(i) = this.linearDistributions(i).mode();
                pdfsAtGridPoints(i) = this.gd.gridValues(i)*this.linearDistributions(i).pdf(mLinModesPos(i));
            end
            [~, ind] = max(pdfsAtGridPoints);
            m = [this.gd.getGridPoint(ind);mLinModesPos(ind)];
        end
        
        function h = plot(this,plotCircleAsInterval)
            arguments 
                this (1,1) StateSpaceSubdivisionDistribution
                plotCircleAsInterval (1,1) logical = false
            end
            if this.linD>2 || this.boundD>3
                error('Cannot plot for this dimension.')
            end
            if ~isa(this.gd,'AbstractCircularDistribution')&&plotCircleAsInterval
                warning('StateSpaceSubdivisionDistribution:invalidPlotArgument','Provided option to plot circle as interval but the GridDistribution is not on the circle.');
            end
            if isa(this.gd,'AbstractCircularDistribution')||isa(this.gd,'AbstractHypertoroidalDistribution')&&this.gd.dim==1
                % Initialize empty and overwrite to bypass normalization
                distBounded = WDDistribution(1,1);
                distBounded.d = this.gd.getGrid();
                distBounded.w = this.gd.gridValues';
            elseif isa(this.gd,'AbstractHypersphericalDistribution')||isa(this.gd,'AbstractHyperhemisphericalDistribution')
                distBounded = HypersphericalDiracDistribution(repmat([1;zeros(this.boundD-1,1)],[1,this.gd.dim]),[1/2,1/4,1/4]); % Initialize with correct dimension
                distBounded.d = this.gd.getGrid();
                distBounded.w = this.gd.gridValues';
            else
                error('This type of distribution is not supported for this dimension.');
            end
            allAxOfFig = findobj(gcf,'Type','Axes');
            if numel(allAxOfFig)>=2 % If there are two axes, use existing
                firstAx = allAxOfFig(end);
                secondAx = allAxOfFig(end-1);
            else % Otherwise, split up figure
                firstAx = subplot(1,2,1);
                secondAx = subplot(1,2,2);
            end
            set(gcf,'CurrentAxes',firstAx);
            cla;
            if isa(this.gd,'AbstractCircularDistribution')||isa(this.gd,'AbstractHypertoroidalDistribution')&&this.gd.dim==1
                if plotCircleAsInterval
                    initialLeftSide = distBounded.plot;
                    [v1,v2] = view;
                    if abs(v1)>1 || abs(v2-90)>1
                        warning('Perspective was not good to see 2-D plot, changing it');
                        view(0, 90);
                    end
                else
                    initialLeftSide = distBounded.plot3d;
                    [v1,v2] = view;
                    if v1==0 && v2==90
                        warning('Perspective was not good to see 3-D plot, changing it');
                        view(17.4,30.4);
                    end
                end
            elseif isa(this.gd,'AbstractHyperhemisphericalDistribution')
                AbstractHyperhemisphericalDistribution.plotHemisphere;
                initialLeftSide = distBounded.plot;
                
            elseif isa(this.gd,'AbstractHypersphericalDistribution')
                AbstractHypersphericalDistribution.plotSphere;
                initialLeftSide = distBounded.plot;
            else
                error('This type of distribution is not supported for this dimesnion.');
            end
            hold on
            
            set(gcf,'CurrentAxes',secondAx);
            cla;
            plotHighlightLeft = [];
            plotPdfRight = [];
            
            % Find all 95 % regegions around the mean for all
            % distributions. Use their union as plotting region so the axis
            % do not change in between.
            scaling = sqrt(chi2inv(0.95,this.dim));
            allMu = cell2mat(arrayfun(@(i){this.linearDistributions(i).mean},1:numel(this.linearDistributions)));
            allCCell = arrayfun(@(i){this.linearDistributions(i).covariance},1:numel(this.linearDistributions));
            allC = cat(3,allCCell{:});
            range = [min(allMu-scaling*reshape(vecnorm(allC,2,1),this.linD,[]),[],2),...
                max(allMu+scaling*reshape(vecnorm(allC,2,1),this.linD,[]),[],2)];
            maxFunVal = max(mvnpdf(allMu',allMu',allC));
            
            for i=1:numel(this.gd.gridValues)
                if ~isempty([plotPdfRight, plotHighlightLeft])
                    delete([plotPdfRight, plotHighlightLeft]);
                end
                set(gcf,'CurrentAxes',firstAx);
                if isa(this.gd,'AbstractCircularDistribution')||isa(this.gd,'AbstractHypertoroidalDistribution')&&this.gd.dim==1
                    grid = this.gd.getGrid(); % Needed because .grid is a function for FIGDistribution
                    if plotCircleAsInterval
                        plotHighlightLeft = stem(grid(i),this.gd.gridValues(i),'r','filled');
                    else
                        plotHighlightLeft = stem3(cos(grid(i)),sin(grid(i)),this.gd.gridValues(i),'Color','r','MarkerFaceColor','r');
                    end
                elseif isa(this.gd,'AbstractHypersphericalDistribution')||isa(this.gd,'AbstractHyperhemisphericalDistribution')
                    currPoint = this.gd.getGridPoint(i);
                    plotHighlightLeft = scatter3(currPoint(1),currPoint(2),currPoint(3),initialLeftSide.SizeData(i),[0.8500, 0.3250, 0.0980],'filled');
                else
                    error('Grid Distribution not supported');
                end
                set(gcf,'CurrentAxes',secondAx);
                
                if this.linD==1
                    plotPdfRight = this.linearDistributions(i).plot(range);
                    ylim([0,maxFunVal]);
                else
                    fun = @(x,y)reshape(this.linearDistributions(i).pdf([x(:)';y(:)']),size(x));
                    plotPdfRight = fsurf(fun,[range(1,:),range(2,:)]);
                    zlim([0,maxFunVal]);
                end
                hold off
                drawnow;
            end
            h = [initialLeftSide, plotPdfRight];
        end
    end
    
    methods (Static)
        function hcrbd = fromDistribution(dist, noOfGridPoints, boundedManifoldType, gridType)
            arguments
                dist AbstractLinBoundedDistribution
                noOfGridPoints {mustBeInteger,mustBePositive}
                boundedManifoldType char = ''
                gridType char = 'default'
            end
            % Do not pass any arg if none is given
            if isempty(boundedManifoldType)&&isa(dist,'AbstractHypercylindricalDistribution')
                boundedManifoldType = 'hypertorus';
            elseif isa(dist,'AbstractHypercylindricalDistribution')
                assert(strcmp(boundedManifoldType,'hypertorus')||strcmp(boundedManifoldType,'circle'))
            elseif isempty(boundedManifoldType)
                error('Cannot auto detect manifold type.');
            end
            hcrbd = StateSpaceSubdivisionDistribution.fromFunction(@(x)dist.pdf(x),...
                noOfGridPoints, dist.boundD, dist.linD, boundedManifoldType, gridType);
        end
        function hcrbd = fromFunction(fun, noOfGridPoints, dimBound, dimLin, boundedManifoldType, gridType)
            arguments
                fun (1,1) function_handle
                noOfGridPoints (1,1) {mustBeInteger,mustBePositive}
                dimBound (1,1) {mustBeInteger,mustBePositive}
                dimLin (1,1) {mustBeInteger,mustBePositive}
                boundedManifoldType char
                gridType char = 'default'
            end
            assert(nargin(fun) == 1, 'Need to be given in the format used for .pdf in densities.');

            % Do not pass any arg if none is given
            if strcmp(gridType,'default')
                typeArg = {};
            else
                typeArg = {gridType};
            end
            if strcmpi(boundedManifoldType,'Hypertorus')&&dimBound==1 || strcmpi(boundedManifoldType,'Circle')
                gd = FIGDistribution.fromDistribution(CircularUniformDistribution, noOfGridPoints);
            elseif strcmpi(boundedManifoldType,'Hypertorus')
                gd = HypertoroidalGridDistribution.fromDistribution(HypertoroidalUniformDistribution(dimBound), noOfGridPoints, typeArg{:});
            elseif strcmpi(boundedManifoldType,'Hypersphere')&&dimBound==3 || strcmpi(boundedManifoldType,'Sphere')
                gd = SphericalGridDistribution.fromDistribution(SphericalUniformDistribution(), noOfGridPoints, typeArg{:});
            elseif strcmpi(boundedManifoldType,'Hypersphere')
                gd = HypersphericalGridDistribution.fromDistribution(HypersphericalUniformDistribution(dimBound), noOfGridPoints, typeArg{:});
            elseif strcmpi(boundedManifoldType,'Hyperhemisphere')&&dimBound==3 || strcmpi(boundedManifoldType,'Hemisphere')
                gd = HemisphericalGridDistribution.fromDistribution(HemisphericalUniformDistribution(), noOfGridPoints, typeArg{:});
            elseif strcmpi(boundedManifoldType,'Hyperhemisphere')
                gd = HyperhemisphericalGridDistribution.fromDistribution(HyperhemisphericalUniformDistribution(dimBound), noOfGridPoints, typeArg{:});
            else
                error('BoundedManifoldType not supported');
            end
            grid = gd.getGrid();
            % Preallocate
            cds = repmat(CustomLinearDistribution(@(x)x,1),[1,noOfGridPoints]);
            
            for i=1:noOfGridPoints
                funCurr = @(y)fun([repmat(grid(:,i),[1,size(y,2)]); y]);
                cds(i) = CustomLinearDistribution(@(x)funCurr(x),dimLin);
                gd.gridValues(i) = cds(i).integral();
                cds(i).scaleBy = 1/gd.gridValues(i);
            end
            
            hcrbd = StateSpaceSubdivisionDistribution(gd,cds);
        end
    end
end

