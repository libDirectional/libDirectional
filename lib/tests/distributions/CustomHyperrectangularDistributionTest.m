classdef CustomHyperrectangularDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testBasic(testCase)
            hud = HyperrectangularUniformDistribution([1,3;2,5]);
            cd = CustomHyperrectangularDistribution(@(x)hud.pdf(x), hud.bounds);
            [xMesh, yMesh] = meshgrid(linspace(1,3,50),linspace(2,5,50));
            testCase.verifyEqual(cd.pdf([xMesh(:)'; yMesh(:)']),1/6*ones(1,50^2));
        end
    end
end