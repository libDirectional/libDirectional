classdef SdHalfCondSdHalfGridDistribution < AbstractConditionalDistribution & AbstractGridDistribution
    % In this class, a conditional distribution on the hyperhemisphere is
    % described by grid values. For the conditional distribution f(a|b), we
    % allow both a and by to vary. Thus, it is a function of the Cartesian
    % product of two hyperhemispheres. It obviously should only integrate to 1 for a
    % fixed b. To provide a grid for the Cartesian product of two hyperhemispheres, we
    % generate the Cartesian product of the grids of the individual hyperhemispheres.   
    methods
        function this = SdHalfCondSdHalfGridDistribution(grid_, gridValues_, enforcePdfNonnegative_)
            % Conditional distribution! First dim conditioned on second
            % dim.
            % Provide grid on the sphere, Cartesian product
            % will be grid on Sd x Sd
            arguments % Use to set default value
                grid_ double {mustBeLessThanOrEqual(grid_,1),mustBeGreaterThanOrEqual(grid_,-1)}
                gridValues_ double {mustBeNonnegative}
                enforcePdfNonnegative_ logical = true
            end
            assert(all(grid_(end,:)>=0), 'Always using upper hemisphere (along last dimension).');
            this.dim = 2*size(grid_,1); % Setting it to twice the dimension of sphere, which is the dimension of the space it is embedded in (convention in libDirectional)
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
            ints = mean(this.gridValues,1)*0.5*AbstractHypersphericalDistribution.computeUnitSphereSurface(this.dim/2);
            if any(abs(ints-1)>tol)
                if all(abs(ints-1)<=tol)
                    error('Normalization:maybeWrongOrder','Not normalized but would be normalized if order of the spheres were swapped. Check input.');
                else
                    warning('Normalization:unnormalized',...
                        'When conditioning values for first sphere on second, normalization is not ensured. One reason may be that you are approximating a density on the entire sphere that is not symmetrical. Otherwise, input or increase tolerance.');
                    warning('Not normalizing for SdHalfCondSdHalfDistribution. Do this manually.');
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
                this SdHalfCondSdHalfGridDistribution
                firstOrSecond (1,1) {mustBeMember(firstOrSecond,[1,2])}
            end
            switch firstOrSecond
                case 1
                    gridValuesSgd = sum(this.gridValues,1)';
                case 2
                    gridValuesSgd = sum(this.gridValues,2);
            end
            sgd = HyperhemisphericalGridDistribution(this.grid, gridValuesSgd);
        end
        
        function sgd = fixDim(this, firstOrSecond, point)
            % Returns slice at point for first or second dimension. Point
            % must be on the grid.
            arguments
                this SdHalfCondSdHalfGridDistribution
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
            sgd = HyperhemisphericalGridDistribution(this.grid,gridValuesSlice);
        end
        
        function plotFullSphere(this)
            if this.dim~=6
                error('Can currently only plot for S2 sphere.')
            end
            hddEqual = HypersphericalDiracDistribution(this.grid, 1/size(this.grid,2)*ones(1,size(this.grid,2)));
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
            AbstractHypersphericalDistribution.plotSphere;
            hold on
            hddEqual.plot;
            set(gcf,'CurrentAxes',secondAx);
            cla;
            AbstractHypersphericalDistribution.plotSphere;
            hold on
            plotHighlightLeft=[];
            plotPdfRight=[];
            for i=1:size(this.grid,2)
                if ~isempty(plotPdfRight)
                    delete([plotPdfRight, plotHighlightLeft]);
                end
                set(gcf,'CurrentAxes',firstAx);
                plotHighlightLeft = scatter3(this.grid(1,i),this.grid(2,i),this.grid(3,i),100,'filled');
                hgdCurr = HyperhemisphericalGridDistribution(this.grid, this.gridValues(:,i));
                hdgdCurr = hgdCurr.toFullSphere;
                hddCurr = HypersphericalDiracDistribution(hdgdCurr.grid, hdgdCurr.gridValues'/sum(hdgdCurr.gridValues));
                set(gcf,'CurrentAxes',secondAx);
                plotPdfRight = hddCurr.plot('r');
                drawnow;
            end
        end
        
        function plotInterpolated(this)
            if this.dim~=6
                error('Can currently only plot for S2 sphere.')
            end
            s2conds2 = S2HalfCondS2HalfGridDistribution(this.grid, this.gridValues, this.enforcePdfNonnegative);
            s2conds2.plotInterpolated();
        end
        
        function plotInterpolatedFullSphere(this)
            if this.dim~=6
                error('Can currently only plot for S2 sphere.')
            end
            s2conds2 = S2HalfCondS2HalfGridDistribution(this.grid, this.gridValues, this.enforcePdfNonnegative);
            s2conds2.plotInterpolatedFullSphere();
        end
        
        function getManifoldSize(~)
            error('Not defined for conditional distributions because interpretation may not be 100% obvious.');
        end
    end
    
    methods (Static)
        function sgd = fromFunction(fun, noGridPoints, funDoesCartesianProduct, gridType, dim)
            arguments
                fun function_handle
                noGridPoints (1,1) {mustBeInteger}
                % State if function does Cartesian product itself. Use this
                % if keeping one argument constant should be done by the
                % function, e.g., because it can be realized more
                % efficiently.
                funDoesCartesianProduct logical = false
                gridType char = 'eq_point_set_symm'
                % Dim as last argument to avoid confusion with
                % S2CondS2distribution.
                dim (1,1) {mustBeInteger,mustBePositive} = 6
            end
            if nargin<5
                warning('Dimension was not explicitly specified. Defaulting to Cartesian product of S2 spheres.');
            end
            % x_{k+1} first argment, x_k second
            assert(nargin(fun)==2);
            switch gridType
                case {'eq_point_set','eq_point_set_symm', 'eq_point_set_symmetric', 'eq_point_set_symm_plane'}
                    grid = eq_point_set_symm(dim/2-1, 2*noGridPoints, true);
                otherwise
                    error('Grid scheme not recognized');
            end
            if ~funDoesCartesianProduct
                [x,y] = meshgrid(1:noGridPoints,1:noGridPoints);
                % Transpose x and y because the after reshaping, the values
                % should iterature through first argument first (first 
                % argument 1 to n, 1 to n, 1 to n, ...). This is just the
                % other way around for meshgrid.
                fvals = fun(grid(:,x'),grid(:,y'));
            else
                fvals = fun(grid,grid);
            end
            if ~funDoesCartesianProduct && isequal(size(fvals),[noGridPoints^dim,noGridPoints^dim])
                error('Function apparently does Cartesian product! Please set funDoesCartesianProduct to true.')
            else
                sgd = SdHalfCondSdHalfGridDistribution(grid, reshape(fvals,[noGridPoints,noGridPoints]));
            end
        end
    end
end

