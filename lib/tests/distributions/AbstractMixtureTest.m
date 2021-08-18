classdef AbstractMixtureTest < matlab.unittest.TestCase
    methods (Test)
        function testSampleMetropolisHastingsBasicsOnlyT2(testCase)
            vmf = ToroidalWNDistribution([1;0],eye(2));
            mix = HypertoroidalMixture({vmf,vmf.shift([1;1])},[0.5,0.5]);
            n = 10;
            s = mix.sampleMetropolisHastings(n);
            testCase.verifySize(s, [mix.dim, n]);
            
            s2 = mix.sample(n);
            testCase.verifySize(s2, [mix.dim, n]);
        end
        function testSampleMetropolisHastingsBasicsOnlyS2(testCase)
            vmf1 = VMFDistribution([1 0]', 2);
            vmf2 = VMFDistribution([0 1]', 2);
            mix = HypersphericalMixture({vmf1,vmf2},[0.5,0.5]);
            n = 10;
            s = mix.sampleMetropolisHastings(n);
            testCase.verifySize(s, [mix.dim, n]);
            testCase.verifyEqual(sum(s.^2,1), ones(1,n), 'RelTol', 1E-10);
            
            s2 = mix.sample(n);
            testCase.verifySize(s2, [mix.dim, n]);
            testCase.verifyEqual(sum(s2.^2,1), ones(1,n), 'RelTol', 1E-10);
        end
        function testSampleMetropolisHastingsBasicsOnlyH2(testCase)
            vmf = VMFDistribution([1 0]', 2);
            mix = CustomHyperhemisphericalDistribution(@(x)vmf.pdf(x)+vmf.pdf(-x),2);
            n = 10;
            s = mix.sampleMetropolisHastings(n);
            testCase.verifySize(s, [mix.dim, n]);
            testCase.verifyEqual(sum(s.^2,1), ones(1,n), 'RelTol', 1E-10);
            
            s2 = mix.sample(n);
            testCase.verifySize(s2, [mix.dim, n]);
            testCase.verifyEqual(sum(s2.^2,1), ones(1,n), 'RelTol', 1E-10);
        end
    end    
   
end
