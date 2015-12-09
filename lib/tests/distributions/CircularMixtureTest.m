classdef CircularMixtureTest < matlab.unittest.TestCase
    properties
    end
    
    methods (Test)                                
        function testCircularMixture(testCase)
            mu1 = 3;
            sigma = 1.5;
            wn = WNDistribution(mu1,sigma);
            mu2 = 2;
            kappa = 0.2;
            vm = VMDistribution(mu2,kappa);
            w1 = 0.3;
            w2 = 1-w1;
            mixture = CircularMixture({wn, vm},[w1, w2]);
            testCase.verifyClass(mixture.dists,'cell');
            testCase.verifyNumElements(mixture.dists,2);
            testCase.verifyTrue(isa(mixture.dists{1},'AbstractCircularDistribution'));
            %% test pdf
            for x= 0:1:2*pi
                testCase.verifyEqual(mixture.pdf(x), w1*wn.pdf(x)+w2*vm.pdf(x), 'RelTol', 1E-10);
            end

            %% test integral
            testCase.verifyEqual(mixture.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(mixture.integralNumerical, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(mixture.integral(0,pi)+mixture.integral(pi,2*pi), 1, 'RelTol', 1E-10);
            testCase.verifyEqual(mixture.integral(-2*pi,4*pi), 3, 'RelTol', 1E-10);
            
            %% test angular moments
            testCase.verifyEqual(mixture.trigonometricMoment(0), mixture.trigonometricMomentNumerical(0), 'RelTol', 1E-5)
            testCase.verifyEqual(mixture.trigonometricMoment(1), mixture.trigonometricMomentNumerical(1), 'RelTol', 1E-5)
            testCase.verifyEqual(mixture.trigonometricMoment(2), mixture.trigonometricMomentNumerical(2), 'RelTol', 1E-5)
            testCase.verifyEqual(mixture.trigonometricMoment(3), mixture.trigonometricMomentNumerical(3), 'RelTol', 1E-5)
           
            %% test sampling
            n = 10;
            s = mixture.sample(10);
            testCase.verifySize(s,[1,n]);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
        end
    end    
   
end
