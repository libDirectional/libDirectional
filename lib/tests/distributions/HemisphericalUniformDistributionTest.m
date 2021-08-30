classdef HemisphericalUniformDistributionTest < matlab.unittest.TestCase
    methods (Test)                                
        function testPdf(testCase)
            sud = HemisphericalUniformDistribution();
            pointsInR3 = randn(3,100);
            pointsOnS2 = pointsInR3./vecnorm(pointsInR3);
            testCase.verifyEqual(sud.pdf(pointsOnS2),1/(2*pi)*ones(1,100));
        end
    end
end
