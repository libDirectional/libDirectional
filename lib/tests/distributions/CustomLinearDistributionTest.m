classdef CustomLinearDistributionTest < matlab.unittest.TestCase
    methods (Test)                                
        function testInitAndMean(testCase)
            g1 = GaussianDistribution([1;1],eye(2));
            g2 = GaussianDistribution([-3;-3],eye(2));
            gm = GaussianMixtureDistribution({g1,g2},[0.7,0.3]);
            
            cld = CustomLinearDistribution.fromDistribution(gm);
            
            CustomLinearDistributionTest.verifyPdfEqual(testCase, cld, gm, 1e-14);
        end
        
        function testIntegral(testCase)
            g1 = GaussianDistribution([1;1],eye(2));
            g2 = GaussianDistribution([-3;-3],eye(2));
            gm = GaussianMixtureDistribution({g1,g2},[0.7,0.3]);
            cld = CustomLinearDistribution.fromDistribution(gm);
            
            testCase.verifyEqual(cld.integral(),1,'AbsTol',1e-10);
        end
        
        function testNormalize(testCase)
            g1 = GaussianDistribution([1;1],eye(2));
            g2 = GaussianDistribution([-3;-3],eye(2));
            gm = GaussianMixtureDistribution({g1,g2},[0.7,0.3]);
            gm.w = gm.w/2;
            cld = CustomLinearDistribution.fromDistribution(gm);
            
            testCase.verifyEqual(cld.integral(),0.5,'AbsTol',1e-10);
        end
    end
    
    methods (Static)
        function verifyPdfEqual(testCase,dist1,dist2,tol)
            [x, y] = meshgrid(linspace(0, 2*pi, 10), linspace(0, 2*pi, 10));
            testCase.verifyEqual(dist1.pdf([x(:)'; y(:)']), dist2.pdf([x(:)'; y(:)']), 'AbsTol', tol);
        end
    end
end
