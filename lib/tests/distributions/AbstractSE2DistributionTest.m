classdef AbstractSE2DistributionTest< matlab.unittest.TestCase
   
    methods (Test)
        function testPlotting(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            se2wn = SE2WNDistribution([1;2;3],diag([1,5,5]));
            
            testCase.applyFixture(SuppressedWarningsFixture('MATLAB:hg:AutoSoftwareOpenGL')); % Supress warnings on jenkins server
            testCase.verifyWarningFree(@()se2wn.plotState());
        end
        
        function testPlotTrajectory(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('MATLAB:hg:AutoSoftwareOpenGL')); % Supress warnings on jenkins server
            testCase.verifyWarningFree(@()AbstractSE2Distribution.plotTrajectory(linspace(0,pi/4,10),[1:10;1:10],animate=true));
            testCase.verifyWarningFree(@()AbstractSE2Distribution.plotTrajectory(linspace(0,pi/4,10),[1:10;1:10],animate=false));
            testCase.verifyWarningFree(@()AbstractSE2Distribution.plotTrajectory(linspace(0,pi/4,10),[1:10;1:10],posColor=[1,0,1],angleColor=[0,1,1]));
            testCase.verifyWarningFree(@()AbstractSE2Distribution.plotTrajectory(linspace(0,pi/4,10),[1:10;1:10],fade=true));
        end
        
        function testAnglePosToDualQuaternion(testCase)
            x = 10*rand(3,10)-5;
            dqs1 = AbstractSE2Distribution.anglePosToDualQuaternion(x);
            dqs2 = NaN(4,size(x,2));
            for i=1:size(x,2)
                se2tmp = AbstractSE2Distribution.anglePosToSE2state(x(:,i));
                dqs2(:,i) = se2tmp.asDualQuaternion;
            end
            testCase.verifyEqual(dqs1,dqs2);
        end
        
        function testDualQuaternionToAnglePos(testCase)
            x = 10*rand(3,10)-5;
            dqs = AbstractSE2Distribution.anglePosToDualQuaternion(x);
            
            [xConverted(1,:),xConverted(2:3,:)] = AbstractSE2Distribution.dualQuaternionToAnglePos(dqs);
            testCase.verifyEqual(xConverted,[mod(x(1,:),2*pi);x(2:end,:)],'AbsTol',1e-15);
        end
    end
end