classdef GSSCDistributionTest < matlab.unittest.TestCase
    
    methods(Test)
        function testPdfSoue(testCase)
            gsscdn2=GSSCDistribution(2,0.3,2,0.9);
            gsscdn2.pdf(linspace(-pi,3*pi,100));
            testCase.verifySize(gsscdn2.pdf(linspace(-pi,3*pi,100)),[1,100]);
            gsscdn3=GSSCDistribution(2,0.3,3,0.9);
            testCase.verifySize(gsscdn3.pdf(linspace(-pi,3*pi,100)),[1,100]);
        end
        function testIntegral(testCase)
            gsscdn2=GSSCDistribution(2,0.3,2,0.9);
            testCase.verifyEqual(gsscdn2.integral(),1,'AbsTol',1e-15);
            gsscdn3=GSSCDistribution(2,0.3,3,0.9);
            testCase.verifyEqual(gsscdn3.integral(),1,'AbsTol',1e-15);
        end
    end
end
