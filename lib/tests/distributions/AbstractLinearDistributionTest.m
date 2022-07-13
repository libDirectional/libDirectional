classdef AbstractLinearDistributionTest< matlab.unittest.TestCase    
    % Test functionality implemented in the abstract class
    methods (Test)
        function testModeNumericalCustom1D(testCase)
            cd = CustomLinearDistribution(@(x)((x>-1)&(x<=0)).*(1+x)+((x>0)&(x<=1)).*(1-x),1);
            cd = cd.shift(0.5);
            testCase.verifyEqual(cd.modeNumerical(),0.5,'AbsTol',1e-5);
        end

        function testModeNumericalGaussian2D(testCase)
            mu = [5,1]';
            C = [2, 1; 1, 1];
            g = GaussianDistribution(mu,C);
            
            testCase.verifyEqual(g.modeNumerical(),mu,'AbsTol',1e-5);
        end
        
        function testModeNumericalGaussian2DMeanFarAway(testCase)
            mu = [5,10]';
            C = [2, 1; 1, 1];
            g = GaussianDistribution(mu,C);
            
            testCase.verifyWarning(@()g.modeNumerical(),'ModeNumerical:StoppedEarly');
            testCase.verifyEqual(g.modeNumerical([4;11]),mu,'AbsTol',1e-5);
        end
        
        function testModeNumericalGaussian3D(testCase)
            mu = [1,2,3]';
            C = [ 1.1 0.4 0; 0.4 0.9 0; 0 0 1];
            g = GaussianDistribution(mu,C);
            
            testCase.verifyEqual(g.modeNumerical(),mu,'AbsTol',1e-6);
        end
        
        function testMeanNumericalGaussian2D(testCase)
            mu = [5,10]';
            C = [2, 1; 1, 1];
            g = GaussianDistribution(mu,C);
            
            testCase.verifyEqual(g.meanNumerical(),mu,'AbsTol',1e-6);
        end
        
        function testMeanNumericalGaussian3D(testCase)
            mu = [1,2,3]';
            C = [ 1.1 0.4 0; 0.4 0.9 0; 0 0 1];
            g = GaussianDistribution(mu,C);
            
            testCase.verifyEqual(g.meanNumerical(),mu,'AbsTol',1e-8);
        end
        
        function testCovarianceNumericalGaussian2D(testCase)
            mu = [5,10]';
            C = [2, 1; 1, 1];
            g = GaussianDistribution(mu,C);
            
            testCase.verifyEqual(g.covarianceNumerical(),C,'AbsTol',1e-6);
        end
        
        function testCovarianceNumericalGaussian3D(testCase)
            mu = [1,2,3]';
            C = [ 1.1 0.4 0; 0.4 0.9 0; 0 0 1];
            g = GaussianDistribution(mu,C);
            
            testCase.verifyEqual(g.covarianceNumerical(),C,'AbsTol',1e-9);
        end

        function testPlotState(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('MATLAB:hg:AutoSoftwareOpenGL')); % Supress warnings on jenkins server
            gd = GaussianDistribution([1;2;3],blkdiag([1,0.5;0.5,1],1));
            testCase.verifyWarningFree(@()gd.plotState());
        end
    end
end