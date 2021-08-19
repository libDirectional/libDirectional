classdef GaussianMixtureDistributionTest < matlab.unittest.TestCase
    methods (Test)                                
        function testInitAndMean(testCase)
            g1 = GaussianDistribution([1;1],eye(2));
            gm = GaussianMixtureDistribution({g1,g1},[0.1,0.9]);
            testCase.verifyEqual(gm.mean,[1;1],'AbsTol',1e-14);
            
            g2 = GaussianDistribution([-3;-3],eye(2));
            gm = GaussianMixtureDistribution({g1,g2},[0.5,0.5]);
            testCase.verifyEqual(gm.mean,[-1;-1],'AbsTol',1e-14);
            
            g1 = GaussianDistribution([1;1;1],eye(3));
            gm = GaussianMixtureDistribution({g1,g1},[0.1,0.9]);
            testCase.verifyEqual(gm.mean,[1;1;1],'AbsTol',1e-14);
            
            g2 = GaussianDistribution([-3;-3;1],eye(3));
            gm = GaussianMixtureDistribution({g1,g2},[0.5,0.5]);
            testCase.verifyEqual(gm.mean,[-1;-1;1],'AbsTol',1e-14);
        end
        
        function testMeanAndCovNumerical(testCase)
            g1 = GaussianDistribution([-1;1],eye(2));
            g2 = GaussianDistribution([-3;-3],[2,-0.5;-0.5,3]);
            gm = GaussianMixtureDistribution({g1,g2},[0.4,0.6]);
            testCase.verifyEqual(gm.mean(),gm.meanNumerical(),'AbsTol',1e-8);
            testCase.verifyEqual(gm.covariance(),gm.covarianceNumerical(),'AbsTol',1e-8);
        end
        
        function testFormulaCov(testCase)
            g1 = GaussianDistribution([-1;1],eye(2));
            g2 = GaussianDistribution([-3;-3],[2,-0.5;-0.5,3]);
            gm = GaussianMixtureDistribution({g1,g2,g1},[0.3,0.3,0.4]);
            gaussArray = [gm.dists{:}];

            mu = [gaussArray.mu]*gm.w';
            zmSamples = [gaussArray.mu] - mu;
            
            covConvexCombination = sum(cat(3,gaussArray.C).*reshape(gm.w,[1,1,numel(gm.w)]),3);
            
            CsampleSet = NaN(gm.dim,gm.dim);
            for j=1:gm.dim
                for k=1:gm.dim
                    Ccurr = 0;
                    for l=1:numel(gm.w)
                        Ccurr = Ccurr + (zmSamples(j,l)*zmSamples(k,l)*gm.w(l));
                    end
                    CsampleSet(j,k)=Ccurr;
                end
            end
            covTotal = covConvexCombination+CsampleSet;
            testCase.verifyEqual(covTotal,gm.covariance,'AbsTol',1e-15);
        end
    end
   
end
