classdef GeneralHypertoroidalMixtureTest < matlab.unittest.TestCase
    properties
    end
    
    methods (Test)                                
        function testGeneralHypertoroidalMixture(testCase)
            mu1 = [3;4];
            C = 2*[1,0.5;0.5,1];
            twn = HypertoroidalWNDistribution(mu1,C);
            mu2 = [1;2];
            kappa1=0.3;kappa2=1.5;
            lambda=0.5;
            tvm=ToroidalVMSineDistribution(mu2,[kappa1;kappa2],lambda);
            w1 = 0.3;
            w2 = 1-w1;
            mixture = GeneralHypertoroidalMixture({twn, tvm},[w1, w2]);
            
            %% test pdf
            testPoints=rand(2,100);
            testCase.verifyEqual(mixture.pdf(testPoints), w1*twn.pdf(testPoints)+w2*tvm.pdf(testPoints), 'RelTol', 1E-10);

            %% test integral
            testCase.verifyEqual(integral2(@(x,y)reshape(mixture.pdf([x(:)';y(:)']),size(x)),0,2*pi,0,2*pi), 1, 'RelTol', 1E-7);
            
            %% test angular moments
            testCase.verifyEqual(mixture.trigonometricMoment(0), mixture.trigonometricMomentNumerical(0), 'RelTol', 1E-5)
            testCase.verifyEqual(mixture.trigonometricMoment(1), mixture.trigonometricMomentNumerical(1), 'RelTol', 1E-5)
            testCase.verifyEqual(mixture.trigonometricMoment(2), mixture.trigonometricMomentNumerical(2), 'RelTol', 1E-5)
            testCase.verifyEqual(mixture.trigonometricMoment(3), mixture.trigonometricMomentNumerical(3), 'RelTol', 1E-5)
           
            %% test sampling
            mixture = GeneralHypertoroidalMixture({tvm,tvm,tvm},[0.3,0.3,0.4]);
            n = 10;
            s = mixture.sample(10);
            testCase.verifySize(s,[mixture.dim,n]);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
        end
    end    
   
end
