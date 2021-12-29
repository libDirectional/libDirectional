classdef SE2UKFMTest < matlab.unittest.TestCase    
    methods (Test)
        function testInitialization(testCase)
            testCase.verifyWarningFree(@()SE2UKFM(3));
        end
        function testSetState(testCase)
            rng default
            p = [1; 1.5];
            ang = 0.4;
            rotMat = [cos(ang),-sin(ang);sin(ang),cos(ang)];
            C = wishrnd(eye(3),3);
            
            iukf = SE2UKFM(3);
            iukf.setState(struct('p',p,'Rot',rotMat),C);
            testCase.verifyEqual(iukf.getPointEstimate(),[ang;p]);
            testCase.verifyEqual(iukf.C,C);
            iukf.setState([ang;p],C);
            testCase.verifyEqual(iukf.getPointEstimate(),[ang;p]);
            testCase.verifyEqual(iukf.C,C);
            
            stateAndCov = iukf.getEstimate();
            testCase.verifyEqual(stateAndCov.state.p,p);
            testCase.verifyEqual(stateAndCov.C,C);
        end
        
        function testPredictGoStraight(testCase)
            % With prediction update cycles, we should converge to a forced
            % (noise-free measurement)
            iukf = SE2UKFM(3,1e-3*ones(3,1));
            iukf.setState([0;0;0],eye(3));
            
            a = @localization_f;

            CBefore = iukf.C;
            iukf.predictNonlinear(a,eye(3),struct('gyro',0,'v',[1;0]));
            testCase.verifyEqual(iukf.getPointEstimate(),[0;1;0]);
            CAfter = iukf.C;
            testCase.verifyGreaterThanOrEqual(eig(CAfter-CBefore),0); % Verify uncertainty got bigger
            
            CBefore = iukf.C;
            iukf.predictNonlinear(a,eye(3),struct('gyro',0,'v',[1;0]));
            testCase.verifyEqual(iukf.getPointEstimate(),[0;2;0]);
            CAfter = iukf.C;
            testCase.verifyGreaterThanOrEqual(eig(CAfter-CBefore),0); % Verify uncertainty got bigger
        end
        
        function testPredictUpdateCyclePosOnly(testCase)
            % With prediction update cycles, we should converge to a forced
            % (noise-free measurement)
            C=[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu=[1;1;1]+pi/2;
            iukf = SE2UKFM(3);
            iukf.setState(mu,C);
            
            forcedPos = [10;20];
            for i=1:100
                iukf.predictIdentity(C);
                testCase.verifySize(iukf.getPointEstimate(), [3,1]);
                iukf.updatePositionMeasurement(0.5*C(2:3,2:3),forcedPos);
            end
            est = iukf.getPointEstimate;
            testCase.verifySize(est, [3,1]);
            testCase.verifyEqual(est(2:3),forcedPos,'AbsTol',0.1);            
        end
        
        function testPredictUpdateCycle(testCase)
            % With prediction update cycles, we should converge to a forced
            % (noise-free measurement)
            C=[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu=[1;1;1]+pi/2;
            iukf = SE2UKFM(3);
            iukf.setState(mu,C);
            
            forcedMeas = [1;10;20];
            for i=1:100
                iukf.predictIdentity(C);
                testCase.verifySize(iukf.getPointEstimate(), [3,1]);
                iukf.updateIdentity(0.5*C,forcedMeas);
            end
            est = iukf.getPointEstimate;
            testCase.verifySize(est, [3,1]);
            testCase.verifyEqual(est,forcedMeas,'AbsTol',0.1);            
        end
    end
end

