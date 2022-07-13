classdef HypercylindricalStateSpaceSubdivisionDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testPloting(testCase)
            n = 100;
            muPeriodic = 4;
            muLinear = [1;2];
            rbd = HypercylindricalStateSpaceSubdivisionDistribution(FIGDistribution.fromDistribution(...
                VMDistribution(muPeriodic,1),n),...
                repmat(GaussianDistribution(muLinear, eye(2)),[n,1]));
            
            testCase.verifyWarningFree(@()rbd.plot());
        end
        
        function testFromDistribution(testCase)
            numberOfGridpoints = 1000;
            hcwn = HypercylindricalWNDistribution([1;1],[5,2;2,2],1);
            hcrbd = HypercylindricalStateSpaceSubdivisionDistribution.fromDistribution(hcwn,numberOfGridpoints);
            
            [xMesh, yMesh] = ndgrid(-4:0.1:4, -4:0.1:4);
            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)']),...
                hcwn.pdf([xMesh(:)';yMesh(:)']),'RelTol',0.002);
        end
    end
        
end
