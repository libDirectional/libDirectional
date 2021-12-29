classdef CartProdStackedDistributionTest < matlab.unittest.TestCase
    
    methods(Test)
        % Test conversions
        function testConstructor(testCase)
            testCase.verifyWarningFree(@()...
                CartProdStackedDistribution({CircularUniformDistribution,GaussianDistribution(0,1)}));
        end
       
        function testSampling(testCase)
            cpd = CartProdStackedDistribution({CircularUniformDistribution,GaussianDistribution(0,1)});
            samples = cpd.sample(100);
            testCase.verifySize(samples, [2,100]);
        end
        
        function testPdf(testCase)
            cpd = CartProdStackedDistribution({CircularUniformDistribution,GaussianDistribution(0,1)});
            testCase.verifySize(cpd.pdf(randn(2,100)), [1,100]);
            % Should be identical for identical inputs
            testCase.verifyEqual(diff(cpd.pdf(ones(2,100))), zeros(1,99));
        end
        
        function testHybridMean(testCase)
            cpd = CartProdStackedDistribution({VMDistribution(3,1),GaussianDistribution(2,1)});
            testCase.verifyEqual(cpd.hybridMean(), [3;2]);
        end
        
        function testShift(testCase)
            cpd = CartProdStackedDistribution({VMDistribution(3,1),GaussianDistribution(2,1)});
            cpdShifted = cpd.shift([2,4]);
            testCase.verifyEqual(cpdShifted.hybridMean(), [5;6]);
        end

        function testSetMode(testCase)
            cpd = CartProdStackedDistribution({VMDistribution(3,1),GaussianDistribution(2,1)});
            cpdShifted = cpd.setMode([2,4]);
            testCase.verifyEqual(cpdShifted.mode(), [2;4]);
        end
    end
end
