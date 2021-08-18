classdef TdCondTdGridDistribution < AbstractConditionalDistribution & AbstractGridDistribution
    % In this class, a conditional distribution on the hypertorus is
    % described by grid values. For the conditional distribution f(a|b), we
    % allow both a and by to vary. Thus, it is a function of the Cartesian
    % product of two hypertori. It obviously should only integrate to 1 for a
    % fixed b. To provide a grid for the Cartesian product of two hypertori, we
    % generate the Cartesian product of the grids of the individual hypertori.
    methods
        function this = TdCondTdGridDistribution(grid_, gridValues_, enforcePdfNonnegative_)
            % Conditional distribution! First dim conditioned on second
            % dim.
            % Provide grid on the sphere, Cartesian product
            % will be grid on Td x Td
            arguments % Use to set default value
                grid_ double {mustBeGreaterThanOrEqual(grid_,0)}
                gridValues_ double {mustBeNonnegative}
                enforcePdfNonnegative_ logical = true
            end
            this.dim = 2*size(grid_,1); % Setting it to twice the dimension of the hypertorus
            assert(size(gridValues_,1)==size(gridValues_,2));
            assert(size(grid_,2)==size(gridValues_,1));
            this.grid = grid_;
            this.gridValues = gridValues_;
            this.enforcePdfNonnegative = enforcePdfNonnegative_;
            this.normalize; % Not this = because we only test and do not actually normalize
        end      
        
        function this = normalize(this)
            % Do not normalize. Returning identical object for
            % compatibility.
            tol = 0.01;
            ints = mean(this.gridValues,1)*(2*pi)^(this.dim/2);
            if any(abs(ints-1)>tol)
                if all(abs(ints-1)<=tol)
                    error('Normalization:maybeWrongOrder','Not normalized but would be normalized if order of the spheres were swapped. Check input.');
                else
                    warning('Normalization:unnormalized',...
                        'When conditioning values for first torus on second, normalization is not ensured. Check input or increase tolerance.');
                    warning('Not normalizing for TdCondTdDistribution. Do this manually.');
                end
            end
            % Could normalize via this.gridValues = this.gridValues./(ints);
        end

        function sdg = multiply(this, other)
            assert(isequal(this.grid,other.grid), 'Multiply:IncompatibleGrid','Can only multiply for equal grids.');
            warning('Multiply:UnnormalizedResult','Multiplication does not yield normalized result.')
            sdg = this;
            sdg.fvals = sdg.fvals.*other.fvals;
        end
        
        function sgd = marginalizeOut(this, firstOrSecond)
            arguments
                this TdCondTdGridDistribution
                firstOrSecond (1,1) {mustBeMember(firstOrSecond,[1,2])}
            end
            switch firstOrSecond
                case 1
                    gridValuesSgd = sum(this.gridValues,1)';
                case 2
                    gridValuesSgd = sum(this.gridValues,2);
            end
            sgd = HypertoroidalGridDistribution(this.grid, gridValuesSgd);
        end
        
        function sgd = fixDim(this, firstOrSecond, point)
            % Returns slice at point for first or second dimension. Point
            % must be on the grid.
            arguments
                this TdCondTdGridDistribution
                firstOrSecond (1,1) {mustBeMember(firstOrSecond,[1,2])}
                point (:,1) double
            end
            assert(size(point,1)==this.dim/2);
            [lia,locb] = ismember(point',this.grid','rows');
            if ~lia
                error('Cannot fix value at this point because it is not on the grid');
            end
            switch firstOrSecond
                case 1
                    gridValuesSlice = this.gridValues(locb,:)';
                case 2
                    gridValuesSlice = this.gridValues(:,locb);
            end
            sgd = HypertoroidalGridDistribution(this.grid,gridValuesSlice);
        end
        
        function plot(this)
            if this.dim>6
                error('Can currently only plot for T1, T2, and T3 torus.')
            end
            hddEqual = HypertoroidalWDDistribution(this.grid, 1/size(this.grid,2)*ones(1,size(this.grid,2)));
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
            hold on
            hddEqual.plot;
            hold on
            set(gcf,'CurrentAxes',secondAx);
            cla;
            plotHighlightLeft=[];
            plotPdfRight=[];
            for i=1:size(this.grid,2)
                if ~isempty(plotPdfRight)
                    delete([plotPdfRight, plotHighlightLeft]);
                end
                set(gcf,'CurrentAxes',firstAx);
                if this.dim==4
                    plotHighlightLeft = scatter(this.grid(1,i),this.grid(2,i),100,[0.8500, 0.3250, 0.0980],'filled');
                elseif this.dim==6     
                    plotHighlightLeft = scatter3(this.grid(1,i),this.grid(2,i),this.grid(3,i),100,[0.8500, 0.3250, 0.0980],'filled');
                end
                hddCurr = HypertoroidalWDDistribution(this.grid, this.gridValues(i,:));
                set(gcf,'CurrentAxes',secondAx);
                plotPdfRight = hddCurr.plot('r');hold on
                drawnow;
            end
        end
        
        function plotInterpolated(this)
            if this.dim>6
                error('Can currently only plot for T1, T2, and T3 torus.')
            end
            hddEqual = HypertoroidalWDDistribution(this.grid, 1/size(this.grid,2)*ones(1,size(this.grid,2)));
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
            hold on
            hddEqual.plot;
            set(gcf,'CurrentAxes',secondAx);
            cla;
            plotHighlightLeft=[];
            plotPdfRight=[];
            for i=1:size(this.grid,2)
                if ~isempty(plotPdfRight)
                    delete([plotPdfRight, plotHighlightLeft]);
                end
                set(gcf,'CurrentAxes',firstAx);
                if this.dim==2
%                     plotHighlightLeft = stem(this.grid(i), this.gridValues
                elseif this.dim==4
                    plotHighlightLeft = scatter(this.grid(1,i),this.grid(2,i),100,[0.8500, 0.3250, 0.0980],'filled');
                elseif this.dim==6     
                    plotHighlightLeft = scatter3(this.grid(1,i),this.grid(2,i),this.grid(3,i),100,[0.8500, 0.3250, 0.0980],'filled');
                end
                sgdCurr = HypertoroidalGridDistribution(this.grid, this.gridValues(:,i));
                set(gcf,'CurrentAxes',secondAx);
                plotPdfRight = sgdCurr.plotInterpolated;
                drawnow;
            end
        end
        
        function getManifoldSize(~)
            error('Not defined for conditional distributions because interpretation may not be 100% obvious.');
        end
    end
    
    methods (Static)
        function sgd = fromFunction(fun, noOfGridPoints, funDoesCartesianProduct, gridType, dim)
            arguments
                fun function_handle
                noOfGridPoints (1,:) {mustBeInteger}
                % State if function does Cartesian product itself. Use this
                % if keeping one argument constant should be done by the
                % function, e.g., because it can be realized more
                % efficiently.
                funDoesCartesianProduct logical
                gridType char
                % Dim as last argument to avoid confusion due to different
                % order of arguments than for S2CondS2distribution.
                dim (1,1) {mustBeInteger,mustBePositive}
            end
            if nargin<5
                warning('Dimension was not explicitly specified. Defaulting to Cartesian product of S2 spheres.');
            end
            % x_{k+1} first argment, x_k second
            assert(nargin(fun)==2);
            assert(mod(dim,2)==0, 'The dimension must be for xk+1 and xk, so it must be even.')
            dimGrid=dim/2;
            switch gridType
                case {'CartesianProd','CartesianProduct'}
                    if numel(noOfGridPoints)==1 % Could considering giving total number of grid points instead
                        noOfGridPoints = repmat(noOfGridPoints,1,dimGrid);% assuming equal resolution
                    end
                    gridIndividualAxis = arrayfun(@(currNo){0:2 * pi / currNo:2 * pi - 2 * pi / currNo}, noOfGridPoints);
                    gridCell = cell(1, length(noOfGridPoints));
                    [gridCell{:}] = ndgrid(gridIndividualAxis{:});
                    
                    gridMatConcat=cat(dimGrid+1,gridCell{:});
                    gridMatForIteration=permute(gridMatConcat,[dimGrid+1,1:dimGrid]);
                    grid=reshape(gridMatForIteration,[dimGrid,prod(noOfGridPoints)]);
                    
                    noOfGridPoints=prod(noOfGridPoints);
                otherwise
                    error('Grid scheme not recognized');
            end
            if ~funDoesCartesianProduct
                [x,y] = meshgrid(1:noOfGridPoints,1:noOfGridPoints);
                % Transpose x and y because the after reshaping, the values
                % should iterature through first argument first (first 
                % argument 1 to n, 1 to n, 1 to n, ...). This is just the
                % other way around for meshgrid.
                fvals = fun(grid(:,x'),grid(:,y'));
            else
                fvals = fun(grid,grid);
            end
            if ~funDoesCartesianProduct && isequal(size(fvals),[noOfGridPoints^dim,noOfGridPoints^dim])
                error('Function apparently does Cartesian product! Please set funDoesCartesianProduct to true.')
            else
                sgd = TdCondTdGridDistribution(grid, reshape(fvals,[noOfGridPoints,noOfGridPoints]));
            end
        end
    end
end

