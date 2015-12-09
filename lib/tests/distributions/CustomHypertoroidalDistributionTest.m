classdef CustomHypertoroidalDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testSimpleDistribution2D(testCase)
            twn = ToroidalWNDistribution([1;2],eye(2));
            chd = CustomHypertoroidalDistribution(@(xa)twn.pdf(xa),2);
            ct1 = CustomToroidalDistribution(@(xa)twn.pdf(xa));
            ct2 = chd.toCustomToroidal();
            [xTest,yTest]=meshgrid(linspace(0,2*pi,50));
            testCase.verifyEqual(chd.pdf([xTest(:)';yTest(:)']),twn.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-10);
            testCase.verifyEqual(ct1.pdf([xTest(:)';yTest(:)']),twn.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-10);
            testCase.verifyEqual(ct2.pdf([xTest(:)';yTest(:)']),twn.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-10);
        end
        
        function testSimpleDistribution1D(testCase)
            wn = WNDistribution(2, 0.7);
            chd = CustomHypertoroidalDistribution(@(xa)wn.pdf(xa),1);
            ct1 = CustomCircularDistribution(@(xa)wn.pdf(xa));
            ct2 = chd.toCustomCircular();
            xTest = linspace(0,2*pi,50);
            testCase.verifyEqual(chd.pdf(xTest),wn.pdf(xTest),'AbsTol',1E-10);
            testCase.verifyEqual(ct1.pdf(xTest),wn.pdf(xTest),'AbsTol',1E-10);
            testCase.verifyEqual(ct2.pdf(xTest),wn.pdf(xTest),'AbsTol',1E-10);
        end        
        
        function testFromDistribution(testCase)
            twn=ToroidalWNDistribution([1;2],eye(2));
            chd=CustomHypertoroidalDistribution.fromDistribution(twn);
            [xTest,yTest]=meshgrid(linspace(0,2*pi,100));
            testCase.verifyEqual(chd.pdf([xTest(:)';yTest(:)']),twn.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-10);
        end
        
        function testShift(testCase)
            twn = ToroidalWNDistribution([3;5],eye(2));
            chd = CustomHypertoroidalDistribution.fromDistribution(twn);
            s = [1;2];
            chd2 = chd.shift(s);
            twn2 = twn.shift(s);
            [xTest,yTest]=meshgrid(linspace(0,2*pi,100));
            testCase.verifyEqual(chd2.pdf([xTest(:)';yTest(:)']),twn.pdf([xTest(:)' - s(1);yTest(:)' - s(2)]),'AbsTol',1E-10);
            testCase.verifyEqual(chd2.pdf([xTest(:)';yTest(:)']),twn2.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-10);
        end
    end
end
