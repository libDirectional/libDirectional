classdef AbstractSE3DistributionTest < matlab.unittest.TestCase
    methods(Test)
        % Test conversions
        function testPlotMode(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('MATLAB:hg:AutoSoftwareOpenGL')); % Supress warnings on jenkins server
            for i=1:10
                cpd = SE3CartProdStackedDistribution(...
                    {HyperhemisphericalWatsonDistribution([1;1;1+i;1]/sqrt(3+(1+i)^2),1),GaussianDistribution([1+i/2;i;0],diag([3,2,3]))});
                testCase.verifyWarningFree(@()cpd.plotMode());
            end
        end
        function testPlotState(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('MATLAB:hg:AutoSoftwareOpenGL')); % Supress warnings on jenkins server
            cpd = SE3CartProdStackedDistribution({HyperhemisphericalWatsonDistribution([1;1;1;1]/sqrt(4),20),GaussianDistribution([1;0;0],blkdiag([1,0.5;0.5,1],1))});
            testCase.verifyWarningFree(@()cpd.plotState());
            testCase.applyFixture(SuppressedWarningsFixture('MATLAB:hg:AutoSoftwareOpenGL')); % Supress warnings on jenkins server
            cpd = SE3CartProdStackedDistribution({HyperhemisphericalWatsonDistribution([1;1;1;1]/sqrt(4),20),GaussianDistribution([10;10;10],blkdiag([1,0.5;0.5,1],1))});
            testCase.verifyWarningFree(@()cpd.plotState(10,false));
        end
        function testPlotTrajectory(testCase)
            close all
            offsets = 0:9;
            quats = [1;1;1;1]+[offsets;zeros(3,10)];
            quats = quats./vecnorm(quats);
            testCase.verifyWarningFree(@()AbstractSE3Distribution.plotTrajectory(quats,[1;2;3]+[offsets;zeros(1,10);offsets]));
            testCase.verifyWarningFree(@()AbstractSE3Distribution.plotTrajectory(quats,[1;2;3]+[offsets;zeros(1,10);offsets],true,0.05));
        end
    end
end
