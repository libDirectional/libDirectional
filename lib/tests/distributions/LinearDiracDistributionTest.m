classdef LinearDiracDistributionTest < matlab.unittest.TestCase
    
    methods(Test)
        % Test conversions
        function testConstructor(testCase)
            testCase.verifyWarningFree(@()LinearDiracDistribution(magic(2),[0.3,0.7]));
        end
        
        function testFromDistribution(testCase)
            rng default
            C = wishrnd(eye(3),3);
            hwn = GaussianDistribution([1;2;3],C);
            hwd = LinearDiracDistribution.fromDistribution(hwn,100000);
            testCase.verifyEqual(hwd.mean, hwn.mean, 'AbsTol',0.005);
            testCase.verifyEqual(hwd.covariance(), hwn.covariance(), 'Rel',0.01);
        end
        
        function testPlotting(testCase)
            ddist = LinearDiracDistribution(-2:2,ones(1,5)/5);
            testCase.verifyWarningFree(@()ddist.plot());
            
            ddist = LinearDiracDistribution(rand(2,5),ones(1,5)/5);
            testCase.verifyWarningFree(@()ddist.plot());
            
            ddist = LinearDiracDistribution(rand(3,5),ones(1,5)/5);
            testCase.verifyWarningFree(@()ddist.plot());
        end
        
        function testMeanAndCov(testCase)
            rng default
            gd = GaussianDistribution([1;2],[2,-0.3;-0.3,1]);
            ddist = LinearDiracDistribution(gd.sample(10000));
            testCase.verifyEqual(ddist.mean(),gd.mean(),'AbsTol',0.005);
            testCase.verifyEqual(ddist.covariance(),gd.covariance(),'AbsTol',0.05);
        end
        
        function testWeightedSamplesToMeanAndCov(testCase)
            % Test against libFiltering
            rng default
            C = wishrnd(eye(3),3);
            hwn = GaussianDistribution([1;2;3],C);
            
            samples = hwn.sample(100);
            weightsUnnorm = rand(1,100);
            weights = weightsUnnorm/sum(weightsUnnorm);
            
            [mu1, C1] = LinearDiracDistribution.weightedSamplesToMeanAndCov(samples, weights);
            [mu2, C2] = Utils.getMeanAndCov(samples, weights); 
            
            testCase.verifyEqual(mu1,mu2);
            testCase.verifyEqual(C1,C2);
        end
    end
end
