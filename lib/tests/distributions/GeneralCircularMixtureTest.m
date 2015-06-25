classdef GeneralCircularMixtureTest < matlab.unittest.TestCase
    properties
    end
    
    methods (Test)                                
        function testGeneralCircularMixture(testCase)
            mu1 = 3;
            sigma1 = 1.5;
            wn1 = WNDistribution(mu1,sigma1);
            mu2 = 2;
            sigma2 = 0.2;
            wn2 = WNDistribution(mu2,sigma2);
            w1 = 0.3;
            w2 = 1-w1;
            mixture = GeneralCircularMixture([wn1, wn2],[w1, w2]);
            
            %% test pdf
            for x= 0:1:2*pi
                testCase.verifyEqual(mixture.pdf(x), w1*wn1.pdf(x)+w2*wn2.pdf(x), 'RelTol', 1E-10);
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
            testCase.verifyEqual(size(s,1), 1);
            testCase.verifyEqual(size(s,2), n);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
        end
    end    
   
end
