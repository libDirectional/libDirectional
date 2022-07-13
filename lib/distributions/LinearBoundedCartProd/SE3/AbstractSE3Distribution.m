classdef (Abstract) AbstractSE3Distribution < AbstractLinHemisphericalDistribution
 
    methods
        function h = plotMode(this)
            arguments
                this (1,1) AbstractSE3Distribution
            end
            % Visualized the mode in SE(3) space
            mode = this.mode();
            h = AbstractSE3Distribution.plotPoint(mode);
        end

        function h = plotState(this, orientationSamples, showMarginalized)
            arguments
                this (1,1) AbstractSE3Distribution
                orientationSamples (1,1) {mustBeInteger,mustBePositive} = 10
                showMarginalized = true
            end
            % Visualized the mode in SE(3) space
            samples = this.sample(orientationSamples);
            mode = this.mode();
            if showMarginalized
                h = this.marginalizePeriodic.plotState;
            else
                h = [];
            end
            for i=1:size(samples,2)
                if showMarginalized
                    linearPart = mode(5:end);
                else
                    linearPart = samples(5:end,i);
                end
                h = [h, AbstractSE3Distribution.plotPoint([samples(1:4,i);linearPart])]; %#ok<AGROW> 
            end
        end
    end
    methods (Static)
        function h = plotPoint(se3point)
            % Visualized just a point in the SE(3) domain (no uncertainties
            % are considered)
            arguments
                se3point (7,1) double % State in 
            end
            rotMat = quaternion(se3point(1:4)').rotmat('point');
            
            pos = se3point(5:end);
            h1 = farrow(pos(1),pos(2),pos(3),pos(1)+rotMat(1,1),pos(2)+rotMat(1,2),pos(3)+rotMat(1,3),'r');
            h2 = farrow(pos(1),pos(2),pos(3),pos(1)+rotMat(2,1),pos(2)+rotMat(2,2),pos(3)+rotMat(2,3),'g');
            h3 = farrow(pos(1),pos(2),pos(3),pos(1)+rotMat(3,1),pos(2)+rotMat(3,2),pos(3)+rotMat(3,3),'b');
            h = [h1, h2, h3];
            neededBoundaries = minmax([pos,pos+rotMat]);
            
            xl = xlim(); yl = ylim(); zl = zlim();
            
            optAx = round(pos/5)*5;
            if xl(1)<neededBoundaries(1,1) || xl(2)>neededBoundaries(1,2)
                xlim([optAx(1)-5,optAx(1)+5]);
            end
            if yl(1)<neededBoundaries(2,1) || yl(2)>neededBoundaries(2,2)
                ylim([optAx(2)-5,optAx(2)+5]);
            end
            if zl(1)<neededBoundaries(3,1) || zl(2)>neededBoundaries(3,2)
                zlim([optAx(3)-5,optAx(3)+5]);
            end
            [ele,~] = view;
            if ele==0
                view(45,45);
            end
        end
        function h = plotTrajectory(periodicStates,linStates,animate,delay)
            arguments
                periodicStates (4,:) double
                linStates (3,:) double
                animate (1,1) logical = false
                delay (1,1) double = 0.05
            end
            assert(size(periodicStates,2)==size(linStates,2));
            h = [];
            for i=1:size(periodicStates,2)
                h = [h, AbstractSE3Distribution.plotPoint([periodicStates(:,i);linStates(:,i)])]; %#ok<AGROW> 
                if animate
                    pause(delay);
                end
            end
        end
    end    
end

