classdef SE2CartProdStackedDistributionTest < matlab.unittest.TestCase
    
    methods(Test)
        % Test conversions
        function testConstructor(testCase)
            testCase.verifyWarningFree(@()...
                SE2CartProdStackedDistribution({CircularUniformDistribution,GaussianDistribution([1;2],diag([3,2]))}));
        end
       
        function testSampling(testCase)
            cpd = SE2CartProdStackedDistribution({CircularUniformDistribution,GaussianDistribution([1;2],diag([3,2]))});
            samples = cpd.sample(100);
            testCase.verifySize(samples, [3,100]);
        end
        
        function testPdf(testCase)
            cpd = SE2CartProdStackedDistribution({CircularUniformDistribution,GaussianDistribution([1;2],diag([3,2]))});
            testCase.verifySize(cpd.pdf(randn(3,100)), [1,100]);
            % Should be identical for identical inputs
            testCase.verifyEqual(diff(cpd.pdf(ones(3,100))), zeros(1,99));
        end
        
        function testHybridMean(testCase)
            cpd = SE2CartProdStackedDistribution({VMDistribution(3,2),GaussianDistribution([1;2],diag([3,2]))});
            testCase.verifyEqual(cpd.hybridMean(), [3;1;2]);
        end

    end
end
