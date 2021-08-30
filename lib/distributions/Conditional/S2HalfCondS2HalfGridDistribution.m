classdef S2HalfCondS2HalfGridDistribution < SdHalfCondSdHalfGridDistribution
    methods
        function this = S2HalfCondS2HalfGridDistribution(grid_, gridValues_, enforcePdfNonnegative_)
            % Conditional distribution! First dim conditioned on second
            % dim.
            % Provide grid on the sphere, Cartesian product
            % will be grid on Sd x Sd
            arguments % Use to set default value
                grid_ (3,:) double {mustBeLessThanOrEqual(grid_,1),mustBeGreaterThanOrEqual(grid_,-1)}
                gridValues_ double {mustBeNonnegative}
                enforcePdfNonnegative_ logical = true
            end
            this@SdHalfCondSdHalfGridDistribution(grid_, gridValues_, enforcePdfNonnegative_)
        end
        
        function plotInterpolated(this)
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
            AbstractHyperhemisphericalDistribution.plotHemisphere;
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
                plotHighlightLeft = scatter3(this.grid(1,i),this.grid(2,i),this.grid(3,i),100,[0.8500, 0.3250, 0.0980],'filled');
                hhgdCurr = HemisphericalGridDistribution(this.grid, this.gridValues(:,i));
                set(gcf,'CurrentAxes',secondAx);
                plotPdfRight = hhgdCurr.plotInterpolated;
                drawnow;
            end
        end
        
        function plotInterpolatedFullSphere(this)
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
            plotHighlightLeft=[];
            plotPdfRight=[];
            warnStruct = warning('off', 'PDF:UseInterpolated');
            for i=1:size(this.grid,2)
                if ~isempty(plotPdfRight)
                    delete([plotPdfRight, plotHighlightLeft]);
                end
                set(gcf,'CurrentAxes',firstAx);
                plotHighlightLeft = scatter3(this.grid(1,i),this.grid(2,i),this.grid(3,i),100,[0.8500, 0.3250, 0.0980],'filled');
                hgdCurr = HemisphericalGridDistribution(this.grid, this.gridValues(:,i));
                hdgdCurr = hgdCurr.toFullSphere;
                set(gcf,'CurrentAxes',secondAx);
                plotPdfRight = hdgdCurr.plotInterpolated;
                drawnow;
            end
            warning(warnStruct);
        end
    end
    methods (Static)
        function s2conds2 = fromFunction(fun, noGridPoints, funDoesCartesianProduct, gridType)
            arguments
                fun function_handle
                noGridPoints (1,1) {mustBeInteger}
                % State if function does Cartesian product itself. Use this
                % if keeping one argument constant should be done by the
                % function, e.g., because it can be realized more
                % efficiently.
                funDoesCartesianProduct logical = false
                gridType char = 'eq_point_set_symm'
            end
            sdCondSd = SdHalfCondSdHalfGridDistribution.fromFunction(fun, noGridPoints, funDoesCartesianProduct, gridType,6);
            s2conds2 = S2HalfCondS2HalfGridDistribution(sdCondSd.grid, sdCondSd.gridValues, sdCondSd.enforcePdfNonnegative);
        end
    end
end