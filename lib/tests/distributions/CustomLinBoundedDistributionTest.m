classdef CustomLinBoundedDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testBasic(testCase)
            % Rectangle [0,1] x R
            gaussian = GaussianDistribution(1,1);
            cd = CustomLinBoundedDistribution(@(x)gaussian.pdf(x(2,:)),1,1);
            xs = linspace(-2,2,100);
            testCase.verifyEqual(cd.pdf([rand(1,100);xs]),gaussian.pdf(xs));
        end
    end
end