classdef AbstractDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testMultiplyOperator(testCase)
            dist = GaussianDistribution([0;1], [1,0.1; 0.1,1]);
            testCase.verifyEqual(dist * dist, dist.multiply(dist))
        end
        function testEqual(testCase)
            dist = GaussianDistribution([0;1], [1,0.1; 0.1,1]);
            testCase.verifyTrue(dist==dist);
        end
        function testUnequal(testCase)
            dist1 = GaussianDistribution([0;1], [1,0.1; 0.1,1]);
            dist2 = GaussianDistribution([0;1]+eps, [1,0.1; 0.1,1]);
            testCase.verifyFalse(dist1==dist2);
        end
    end
end