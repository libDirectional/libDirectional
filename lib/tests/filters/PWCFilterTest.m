classdef PWCFilterTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testPWCFilter(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true; end          
            L = 2;
            filter = PWCFilter(L);
            wn = WNDistribution(1.3, 0.8);
            w = PWCDistribution.calculateParametersNumerically(@(x) wn.pdf(x),L);
            pwcd = PWCDistribution(w);
            
            % sanity check
            filter.setState(pwcd);
            pwcd1 = filter.getEstimate();
            testCase.verifyClass(pwcd1, 'PWCDistribution');
            testCase.verifyEqual(pwcd.w, pwcd1.w);
                       
            % check calculation of A
            a = @(x,w) mod(x+w,2*pi);
            if enableExpensive
                % todo check row/col
                A = PWCFilter.calculateSystemMatrixNumerically(L, a, wn);
                testCase.verifyEqual(sum(A,2), ones(L,1), 'RelTol', 1E-4);
                testCase.verifyGreaterThanOrEqual(A, zeros(L,L));
            end
            
            % check calculation of A with noise samples
            wd = wn.toDirac3();
            A2 = PWCFilter.calculateSystemMatrixNumericallySampledNoise(L, a, wd.d, wd.w);
            testCase.verifyEqual(sum(A2,2), ones(L,1), 'RelTol', 1E-4);
            testCase.verifyGreaterThanOrEqual(A2, zeros(L,L));
            
            % check calculation of H
            h = @(x,w) mod(x+w,2*pi);
            Lmeas = 3;
            if enableExpensive
                % todo check row/col
                H = PWCFilter.calculateMeasurementMatrixNumerically(L, Lmeas, h, wn);
                testCase.verifyEqual(sum(H,1), ones(1,L), 'RelTol', 1E-4);
                testCase.verifyGreaterThanOrEqual(H, zeros(Lmeas,L));
            end
            
            wd = wn.toDirac3();
            H2 = PWCFilter.calculateMeasurementMatrixNumericallySampledNoise(L, Lmeas, a, wd.d, wd.w);
            testCase.verifyEqual(sum(H2,1), ones(1,L), 'RelTol', 1E-4);
            testCase.verifyGreaterThanOrEqual(H2, zeros(Lmeas,L));
            
            likelihood = LikelihoodFactory.additiveNoiseLikelihood(@(x) x, wn);
            %H = PWCFilter.calculateMeasurementMatrixNumerically(n, likelihood);
            %testCase.verifyEqual(sum(B), ones(1,n), 'RelTol', 1E-10);
            
            % check prediction
            filter.setState(pwcd);
            filter.predict(A2);
            
            % check update
            filter.setState(pwcd);
            filter.update(H2,1);
            
            filter.setState(pwcd);
            filter.updateLikelihood(likelihood, 1);
        end
        
        function testUpdateLikelihood(testCase)    
            prior = VMDistribution(2,3);
            vmNoise = VMDistribution(0,1.7);
            h = @(x) x;
            likelihood = LikelihoodFactory.additiveNoiseLikelihood(h, vmNoise);
            z = 2.7;
            L = 1000;

            pwcFilter = PWCFilter(L);
            pwcFilter.setState(prior);
            pwcFilter.updateLikelihood(likelihood, z);
            pwc = pwcFilter.getEstimate();

            %VM filter is exact in this scenario
            vmFilter = VMFilter();
            vmFilter.setState(prior);
            vmFilter.updateIdentity(vmNoise, z);
            vm = vmFilter.getEstimate(); 
            
            testCase.verifyEqual(pwc.trigonometricMoment(1), vm.trigonometricMoment(1), 'RelTol', 1E-4);
            testCase.verifyEqual(pwc.trigonometricMoment(2), vm.trigonometricMoment(2), 'RelTol', 1E-4);
            testCase.verifyEqual(pwc.trigonometricMoment(3), vm.trigonometricMoment(3), 'RelTol', 1E-4);
        end
        
        function testUpdate(testCase)
            prior = VMDistribution(2,3);
            vmNoise = VMDistribution(0,1.7);
            h = @(x) x;
            likelihood = LikelihoodFactory.additiveNoiseLikelihood(h, vmNoise);
            z = 2.7;
            L = 100;
            Lmeas = 110;

            pwcFilter = PWCFilter(L);
            pwcFilter.setState(prior);
            H = PWCFilter.calculateMeasurementMatrixNumericallyFromLikelihood(L, Lmeas, likelihood);
            pwcFilter.update(H, z);
            pwc = pwcFilter.getEstimate();

            %VM filter is exact in this scenario
            vmFilter = VMFilter();
            vmFilter.setState(prior);
            vmFilter.updateIdentity(vmNoise, z);
            vm = vmFilter.getEstimate(); 
            
            testCase.verifyEqual(pwc.trigonometricMoment(1), vm.trigonometricMoment(1), 'RelTol', 1E-2);
            testCase.verifyEqual(pwc.trigonometricMoment(2), vm.trigonometricMoment(2), 'RelTol', 1E-2);
            testCase.verifyEqual(pwc.trigonometricMoment(3), vm.trigonometricMoment(3), 'RelTol', 1E-1);
            
            %compare to H matrix obtained using sampled noise
            %wdNoise = vmNoise.toDiracBT(100);
            %H2 = PWCFilter.calculateMeasurementMatrixNumericallySampledNoise(L, Lmeas, @(x,v) mod(x+v,2*pi), wdNoise.d, wdNoise.w);
            %testCase.verifyEqual(H,H2, 'AbsTol', 0.15) %H2 is not very accurate because of the limited number of samples and other artifacts
        end
        
        function testCalculateSystemMatrixNumerically(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true; end   
            
            if enableExpensive
                L = 2;
                a = @(x,w) mod(x+w,2*pi);
                noiseDistribution = WNDistribution(2,3);
                A = PWCFilter.calculateSystemMatrixNumerically(L, a, noiseDistribution);
                testCase.verifySize(A, [L,L]);
                testCase.verifyEqual(sum(A), ones(1,L), 'RelTol', 1E-5);
            end
        end
        
        function testCalculateMeasurementMatrixNumerically(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true; end   
            
            if enableExpensive
                L = 2;
                Lmeas = 3;
                h = @(x,w) mod(x+w,2*pi);
                noiseDistribution = WNDistribution(2,3);
                H = PWCFilter.calculateMeasurementMatrixNumerically(L, Lmeas, h, noiseDistribution);
                testCase.verifySize(H, [Lmeas,L]);
                testCase.verifyEqual(sum(H), ones(1,L), 'RelTol', 1E-5);
            end
        end
    end
end
