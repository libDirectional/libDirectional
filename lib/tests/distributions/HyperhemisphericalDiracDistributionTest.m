classdef HyperhemisphericalDiracDistributionTest < matlab.unittest.TestCase
    methods (Test)
        function testHyperhemisphericalDiracDistribution(testCase)
            d = [0.5,3,4,6,6;
                 2,2,5,3,0;
                 0.5,0.2,5.8,4.3,1.2];
            d = d./repmat(sqrt(sum(d.^2)), 3, 1);
            w = [0.1 0.1 0.1 0.1 0.6];
            hdd = HyperhemisphericalDiracDistribution(d,w);
            
            % test errors
            testCase.verifyError(@() hdd.pdf(1), 'PDF:UNDEFINED');
            testCase.verifyError(@() hdd.sampleMetropolisHastings(1), 'PDF:UNDEFINED');
            testCase.verifyError(@() hdd.entropyNumerical(), 'PDF:UNDEFINED');
            testCase.verifyError(@() hdd.integralNumerical(), 'PDF:UNDEFINED');
            
            % sanity check
            testCase.verifyClass(hdd, 'HyperhemisphericalDiracDistribution');
            testCase.verifyEqual(hdd.d, d);
            testCase.verifyEqual(hdd.w, w);
            testCase.verifyEqual(hdd.integral, 1, 'RelTol', 1E-10);
            
            % test sampling
            nSamples = 5;
            s = hdd.sample(nSamples);
            testCase.verifySize(s,[size(d,1),nSamples]);
            testCase.verifyEqual(s, mod(s,2*pi))

            % test apply function 
            same = hdd.applyFunction(@(x) x);
            testCase.verifyEqual(same.d, hdd.d, 'RelTol', 1E-10);
            testCase.verifyEqual(same.w, hdd.w, 'RelTol', 1E-10);
            mirrored = hdd.applyFunction(@(x) -x );
            testCase.verifyEqual(mirrored.d, -hdd.d, 'RelTol', 1E-10);
            testCase.verifyEqual(mirrored.w, hdd.w, 'RelTol', 1E-10);
            
            % test reweigh            
            f = @(x) 2*ones(1,size(x,2)); %does not change anything because of renormalization
            twdRew = hdd.reweigh(f);
            testCase.verifyClass(twdRew, 'HyperhemisphericalDiracDistribution');
            testCase.verifyEqual(twdRew.d, hdd.d);
            testCase.verifyEqual(twdRew.w, hdd.w);
            
            f = @(x) x(1,:);
            twdRew = hdd.reweigh(f);
            testCase.verifyClass(twdRew, 'HyperhemisphericalDiracDistribution');
            testCase.verifyEqual(twdRew.d, hdd.d);
            wNew = hdd.d(1,:).*hdd.w;
            testCase.verifyEqual(twdRew.w, wNew/sum(wNew));
            testCase.verifyEqual(twdRew.integral, 1, 'RelTol', 1E-10);
            
            % test entropy
            testCase.verifyWarning(@hdd.entropy, 'ENTROPY:DISCRETE')
        end

        function testMoment(testCase)
            q = [1,2,3,4]';
            q = q/norm(q);
            M = [quaternionMultiplication(q, [1 0 0 0]'), quaternionMultiplication(q, [0 1 0 0]'), quaternionMultiplication(q, [0 0 1 0]'), quaternionMultiplication(q, [0 0 0 1]')];
            Z = [-10 -2 -1 0]';
            B = HyperhemisphericalBinghamDistribution(Z,M);

            hdd = HyperhemisphericalDiracDistribution.fromDistribution(B, 100);
        end
    end
end