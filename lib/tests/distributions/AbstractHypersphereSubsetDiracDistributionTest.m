classdef AbstractHypersphereSubsetDiracDistributionTest < matlab.unittest.TestCase    
    methods(Test)
        % Test conversions
        function testMeanAxis(testCase)
            rng default
            q = [1,2,3,4]';
            q = q/norm(q);
            M = [quaternionMultiplication(q, [1 0 0 0]'), quaternionMultiplication(q, [0 1 0 0]'), quaternionMultiplication(q, [0 0 1 0]'), quaternionMultiplication(q, [0 0 0 1]')];
            Z = [-10 -2 -1 0]';
            
            bd = BinghamDistribution(Z,M);
            bdHemi = HyperhemisphericalBinghamDistribution(Z,M);
            wdFull = HypersphericalDiracDistribution.fromDistribution(bd,100001);
            wdHemi = HyperhemisphericalDiracDistribution.fromDistribution(bdHemi,100001);
            
            testCase.verifyEqual(wdFull.meanAxis(),bd.meanAxis(),'AbsTol',0.01);
            testCase.verifyEqual(wdHemi.meanAxis(),bd.meanAxis(),'AbsTol',0.01);
        end
    end
end
