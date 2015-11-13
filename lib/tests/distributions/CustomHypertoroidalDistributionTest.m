classdef CustomHypertoroidalDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testSimpleDistribution(testCase)
            twn=ToroidalWNDistribution([1;2],eye(2));
            chd=CustomHypertoroidalDistribution(@(xa)twn.pdf(xa),2);
            [xTest,yTest]=meshgrid(linspace(0,2*pi,100));
            testCase.verifyEqual(chd.pdf([xTest(:)';yTest(:)']),twn.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-10);
        end
        function testFromDistribution(testCase)
            twn=ToroidalWNDistribution([1;2],eye(2));
            chd=CustomHypertoroidalDistribution.fromDistribution(twn);
            [xTest,yTest]=meshgrid(linspace(0,2*pi,100));
            testCase.verifyEqual(chd.pdf([xTest(:)';yTest(:)']),twn.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-10);
        end
    end
end
