classdef HypersphericalDiracDistributionTest < matlab.unittest.TestCase
    properties
    end
    
    methods (Test)
        function testHypersphericalDiracDistribution(testCase)
            d = [0.5,3,4,6,6;
                 2,2,5,3,0;
                 0.5,0.2,5.8,4.3,1.2];
            d = d./repmat(sqrt(sum(d.^2)), 3, 1);
            w = [0.1 0.1 0.1 0.1 0.6];
            hdd = HypersphericalDiracDistribution(d,w);
            
            % test errors
            testCase.verifyError(@() hdd.pdf(1), 'PDF:UNDEFINED');
            testCase.verifyError(@() hdd.sampleMetropolisHastings(1), 'PDF:UNDEFINED');
            testCase.verifyError(@() hdd.entropyNumerical(), 'PDF:UNDEFINED');
            testCase.verifyError(@() hdd.integralNumerical(), 'PDF:UNDEFINED');
            
            % sanity check
            testCase.verifyClass(hdd, 'HypersphericalDiracDistribution');
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
            testCase.verifyClass(twdRew, 'HypersphericalDiracDistribution');
            testCase.verifyEqual(twdRew.d, hdd.d);
            testCase.verifyEqual(twdRew.w, hdd.w);
            
            f = @(x) x(1,:);
            twdRew = hdd.reweigh(f);
            testCase.verifyClass(twdRew, 'HypersphericalDiracDistribution');
            testCase.verifyEqual(twdRew.d, hdd.d);
            wNew = hdd.d(1,:).*hdd.w;
            testCase.verifyEqual(twdRew.w, wNew/sum(wNew));
            testCase.verifyEqual(twdRew.integral, 1, 'RelTol', 1E-10);
            
            % test entropy
            testCase.verifyWarning(@hdd.entropy, 'ENTROPY:DISCRETE')
        end
        
        function testToWD(testCase)
            rng default
            n = 20;
            d = 2*pi*rand(2,n);
            d = d./repmat(sqrt(sum(d.^2)), 2, 1);
            w = rand(1,n);
            hdd = HypersphericalDiracDistribution(d,w/sum(w));
            wd = hdd.toWD();
            testCase.verifyClass(wd, 'WDDistribution');
            testCase.verifyEqual(cos(wd.d), hdd.d(1,:), 'RelTol', 1E-10);
            testCase.verifyEqual(sin(wd.d), hdd.d(2,:), 'RelTol', 1E-10);
            testCase.verifyEqual(wd.w, hdd.w, 'RelTol', 1E-10);
        end
        
        function testFromDistribution(testCase)
            rng default
            vmf = VMFDistribution([1;1;1]/sqrt(3),1);
            wd = HypersphericalDiracDistribution.fromDistribution(vmf,50000);
            testCase.verifyEqual(wd.meanDirection, vmf.meanDirection, 'AbsTol',0.01);
        end
    end
end