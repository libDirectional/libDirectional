classdef AbstractHypercylindricalDistributionTest< matlab.unittest.TestCase    
    % Test functionality implemented in the abstract class
    methods (Test)
        function testModeNumericalGaussian2D(testCase)
            mu = [5,1]';
            C = [2, 1; 1, 1];
            g = HypercylindricalWNDistribution(mu, C, 1);
            
            testCase.verifyEqual(g.modeNumerical(),mu,'AbsTol',1e-6);
        end
        
        function testModeNumericalGaussian2DMeanFarAway(testCase)
            mu = [5,10]';
            C = [2, 1; 1, 1];
            g = HypercylindricalWNDistribution(mu,C,1);
            
            testCase.verifyWarning(@()g.modeNumerical(),'ModeNumerical:StoppedEarly');
            testCase.verifyEqual(g.modeNumerical([4;11]),mu,'AbsTol',1e-6);
        end
        
        function testModeNumericalGaussian3D(testCase)
            mu = [1,2,3]';
            C = [ 1.1 0.4 0; 0.4 0.9 0; 0 0 1];
            g = HypercylindricalWNDistribution(mu,C,1);
            
            testCase.verifyEqual(g.modeNumerical(),mu,'AbsTol',1e-5);
        end
        
        function testHybridMomentNumericalGaussian2D(testCase)
            mu = [5,10]';
            C = [2, 1; 1, 1];
            g = HypercylindricalWNDistribution(mu,C,1);
            
            testCase.verifyEqual(g.hybridMomentNumerical(),g.hybridMoment(),'AbsTol',1e-7);
        end
        
        function testHybridMomentNumericalGaussian3D(testCase)
            mu = [1,2,3]';
            C = [ 1.1 0.4 0; 0.4 0.9 0; 0 0 1];
            g = HypercylindricalWNDistribution(mu,C,1);
            
            testCase.verifyEqual(g.hybridMomentNumerical(),g.hybridMoment(),'AbsTol',1e-7);
        end
        
        function testHybridMeanNumericalGaussian2D(testCase)
            mu = [5,10]';
            C = [2, 1; 1, 1];
            g = HypercylindricalWNDistribution(mu,C,1);
            
            testCase.verifyEqual(g.hybridMeanNumerical(),mu,'AbsTol',1e-7);
        end
        
        function testHybridMeanNumericalGaussian3D(testCase)
            mu = [1,2,3]';
            C = [ 1.1 0.4 0; 0.4 0.9 0; 0 0 1];
            g = HypercylindricalWNDistribution(mu,C,1);
            
            testCase.verifyEqual(g.hybridMeanNumerical(),mu,'AbsTol',1e-7);
        end
        
        function testLinearMeanNumerical(testCase)
            hwn = HypercylindricalWNDistribution([1,2],[2,0.3;0.3,1],1);
            testCase.verifyEqual(hwn.linearMeanNumerical(), 2,'AbsTol', 1e-10);
        end

        function testConditionOnPeriodic(testCase)
            hwn = HypercylindricalWNDistribution([1,2],[2,0.3;0.3,1],1);
            dist = hwn.conditionOnPeriodic(1.5);
            testCase.verifyEqual(diff(hwn.pdf([1.5*ones(1,11);-5:5])./dist.pdf(-5:5)), zeros(1,10),'AbsTol', 1e-10);
            dist = hwn.conditionOnPeriodic(1.5+2*pi);
            testCase.verifyEqual(diff(hwn.pdf([1.5*ones(1,11);-5:5])./dist.pdf(-5:5)), zeros(1,10),'AbsTol', 1e-10);
        end

        function testConditionOnLinear(testCase)
            hwn = HypercylindricalWNDistribution([1,2],[2,0.3;0.3,1],1);
            dist = hwn.conditionOnLinear(1.5);
            testCase.verifyEqual(diff(hwn.pdf([-5:5;1.5*ones(1,11)])./dist.pdf(-5:5)), zeros(1,10),'AbsTol', 1e-10);
            dist = hwn.conditionOnLinear(1.5+2*pi);
            testCase.verifyNotEqual(diff(hwn.pdf([-5:5;1.5*ones(1,11)])./dist.pdf(-5:5)), zeros(1,10));
        end
    end
end