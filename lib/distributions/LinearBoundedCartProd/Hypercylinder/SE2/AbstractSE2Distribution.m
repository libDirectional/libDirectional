classdef (Abstract) AbstractSE2Distribution < AbstractHypercylindricalDistribution
 
    methods        
        function h = plotSE2State(this, scalingFactor, circleColor, angleColor)
            arguments
                this (1,1) AbstractSE2Distribution
                scalingFactor (1,1) double = 1
                circleColor (3,1) double = [0    0.4470    0.7410]
                angleColor (3,1) double = [0.8500    0.3250    0.0980]
            end
            linearCovmat = this.linearCovariance;
            hybridMoment = this.hybridMoment;
            linearMean = hybridMoment(3:4);
            periocMean = atan2(hybridMoment(2),hybridMoment(1));
            periodicVar = 1-norm(hybridMoment(1:2));
            
            holdStatus = ishold();
            hold on
            xs = [linspace(0,2*pi,100),0];
            ps = scalingFactor*linearCovmat*[cos(xs);sin(xs)];
            h(1) = plot(ps(1,:)+linearMean(1),ps(2,:)+linearMean(2),'color',circleColor);hold on
            
            plotAngRange = 0.1*periodicVar*pi; % 0.1*variance (between 0 and 1) * pi (max -pi to +pi)
            xs = linspace(periocMean-plotAngRange, periocMean+plotAngRange,100);
            ps = scalingFactor*linearCovmat*[cos(xs);sin(xs)];
            scaledMeanVec = scalingFactor*linearCovmat*[cos(periocMean);sin(periocMean)];
            h(2) = quiver(linearMean(1),linearMean(2),scaledMeanVec(1),scaledMeanVec(2),'AutoScaleFactor',1,'color',angleColor);
            h(3) = plot(linearMean(1)+ps(1,:),linearMean(2)+ps(2,:),'--','color',angleColor);
            h(4) = plot([linearMean(1),linearMean(1)+ps(1,end)],[linearMean(2),linearMean(2)+ps(2,end)],'-','color',angleColor);
            h(5) = plot([linearMean(1),linearMean(1)+ps(1,1)],[linearMean(2),linearMean(2)+ps(2,1)],'-','color',angleColor);

            if ~holdStatus
                hold off
            end
        end
    end
    methods (Static)
        function dq = anglePosToDualQuaternion(x)
                % Convert a matrix with angle-pos along the first dimension to
                % dual quaternions
                arguments
                    x (3,:) double
                end
                rot = x(1,:);
                trans = x(2:3,:);
                dq_real = [cos(rot/2);sin(rot/2)];
                dq_real_tmp = [-dq_real(2,:);dq_real(1,:)];
                dq_dual = 0.5*[sum(dq_real.*trans,1);sum(dq_real_tmp.*trans,1)];
                dq = [dq_real;dq_dual];
        end
        function se2State = anglePosToSE2state(x)
            % Convert a single state into an object of a class that can
            % provide a variety of different representations for the SE2
            % state
            % Input: 
            %   rot (1xn vector): angles in radian
            %   trans (2xn matrix): translation vectors concatenated columnwise
            % Output:
            %   dq (4xn matrix): dual quaternions concatenated columnwise
            arguments
                x(3,1) double
            end
            se2State = SE2(x(1),x(2:end));
        end
        function [angleOrAnglePos, pos] = dualQuaternionToAnglePos(dq)
            arguments
                dq (4,:) double
            end
            % Input:
            %   dq (4xn matrix): dual quaternions concatenated columnwise
            % Output: 
            %   rot (1xn vector): angles in radian
            %   trans (2xn matrix): translation vectors concatenated columnwise
            angleOrAnglePos = mod(2*atan2(dq(2,:),dq(1,:)), 2*pi);
            q_1 = [dq(1,:);-dq(2,:)];
            q_2 = [dq(2,:);dq(1,:)];
            pos = 2*[sum(q_1.*dq(3:end,:),1);sum(q_2.*dq(3:end,:),1)];
            if nargout<=1
                angleOrAnglePos = [angleOrAnglePos;pos];
            end
        end
        function h = plotSE2trajectory(periodicStates,linStates,animate,delay)
            arguments
                periodicStates (1,:) double
                linStates (2,:) double
                animate (1,1) logical = false
                delay (1,1) double = 0.05
            end
            if ~animate
                h1 = quiver(linStates(1,:),linStates(2,:),cos(periodicStates),sin(periodicStates));
                h2 = scatter(linStates(1,:),linStates(2,:));
                h = [h1, h2];
                return
            end

            %First, set the axis right
            h = quiver(linStates(1,:),linStates(2,:),cos(periodicStates),sin(periodicStates),'AutoScale','off');
            % There is apparently no way to use the optimal scaling factor of matlab
            % and just erase some arrows. So we use manual scaling
            drawnow
            axis manual
            xlimEnd = xlim;
            ylimEnd = ylim;

            delete(h);

            for i = 1:numel(periodicStates)-1
                h1 = quiver(linStates(1,1:i),linStates(2,1:i),cos(periodicStates(1:i)),sin(periodicStates(1:i)),'AutoScale','off');
                h2 = scatter(linStates(1,1:i),linStates(2,1:i));
                xlim(xlimEnd);
                ylim(ylimEnd);
                pause(delay)
                delete([h1,h2]);
            end
            h1 = quiver(linStates(1,:),linStates(2,:),cos(periodicStates),sin(periodicStates),'AutoScale','off');
            h2 = scatter(linStates(1,:),linStates(2,:));
            h = [h1, h2];
        end
    end
    
end

