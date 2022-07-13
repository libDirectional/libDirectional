classdef SE3CartProdStackedDistributionTest < matlab.unittest.TestCase
    
    methods(Test)
        % Test conversions
        function testConstructor(testCase)
            testCase.verifyWarningFree(@()...
                SE3CartProdStackedDistribution({HyperhemisphericalUniformDistribution(4),GaussianDistribution([1;2;3],diag([3,2,1]))}));
        end
        
        function testSampling(testCase)
            cpd = SE3CartProdStackedDistribution({HyperhemisphericalUniformDistribution(4),GaussianDistribution([1;2;0],diag([3,2,3]))});
            samples = cpd.sample(100);
            testCase.verifySize(samples, [7,100]);
        end
        
        function testPdf(testCase)
            cpd = SE3CartProdStackedDistribution({HyperhemisphericalUniformDistribution(4),GaussianDistribution([1;2;0],diag([3,2,3]))});
            testCase.verifySize(cpd.pdf(randn(7,100)), [1,100]);
            % Should be identical for identical inputs
            testCase.verifyEqual(diff(cpd.pdf(ones(7,100))), zeros(1,99));
        end
        
        function testMode(testCase)
            watson = HyperhemisphericalWatsonDistribution([2;1;3;1]/norm([2;1;3;1]),2);
            gaussian = GaussianDistribution([1;2;0],diag([3,2,3]));
            cpd = SE3CartProdStackedDistribution({watson,gaussian});
            testCase.verifyEqual(cpd.mode(), [watson.mode();gaussian.mode()]);
        end
    end
end
