classdef HypertoroidalMixtureTest < matlab.unittest.TestCase
    properties
    end
    
    methods (Test)                                
        function test2D(testCase)
            mu1 = [3;4];
            C = 2*[1,0.5;0.5,1];
            twn = ToroidalWNDistribution(mu1,C);
            mu2 = [1;2];
            kappa1 = 0.3;
            kappa2 = 1.5;
            lambda = 0.5;
            tvm = ToroidalVMSineDistribution(mu2,[kappa1;kappa2],lambda);
            w1 = 0.3;
            w2 = 1-w1;
            mixture = HypertoroidalMixture({twn, tvm},[w1, w2]);
            testCase.verifyClass(mixture, 'HypertoroidalMixture');
            mixture2 = mixture.toToroidalMixture();
            testCase.verifyClass(mixture2, 'ToroidalMixture');
            
            % test pdf
            rng default
            testPoints = rand(2,100);
            testCase.verifyEqual(mixture.pdf(testPoints), w1*twn.pdf(testPoints)+w2*tvm.pdf(testPoints), 'RelTol', 1E-10);
            testCase.verifyEqual(mixture.pdf(testPoints), mixture2.pdf(testPoints), 'RelTol', 1E-10);

            % test integral
            testCase.verifyEqual(integral2(@(x,y)reshape(mixture.pdf([x(:)';y(:)']),size(x)),0,2*pi,0,2*pi), 1, 'RelTol', 1E-7);
            testCase.verifyEqual(mixture.integral(), 1, 'RelTol', 1E-7);
            
            % test angular moments
            testCase.verifyEqual(mixture.trigonometricMoment(0), mixture.trigonometricMomentNumerical(0), 'RelTol', 1E-5)
            testCase.verifyEqual(mixture.trigonometricMoment(1), mixture.trigonometricMomentNumerical(1), 'RelTol', 1E-5)
            testCase.verifyEqual(mixture.trigonometricMoment(2), mixture.trigonometricMomentNumerical(2), 'RelTol', 1E-5)
            testCase.verifyEqual(mixture.trigonometricMoment(3), mixture.trigonometricMomentNumerical(3), 'RelTol', 1E-5)
            
            % test sampling
            mixture = ToroidalMixture({tvm,tvm,tvm},[0.3,0.3,0.4]);
            n = 10;
            s = mixture.sample(10);
            testCase.verifySize(s,[mixture.dim,n]);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
        end
        
        function test1D(testCase)
            dist1 = CircularUniformDistribution();
            dist2 = WCDistribution(2,0.9);
            w1 = 0.2;
            w2 = 1-w1;
            mixture = HypertoroidalMixture({dist1, dist2},[w1, w2]);
            mixture2 = mixture.toCircularMixture();
            testCase.verifyClass(mixture2, 'CircularMixture');
           
            % test pdf
            rng default
            testPoints = rand(1,100);
            testCase.verifyEqual(mixture.pdf(testPoints), w1*dist1.pdf(testPoints)+w2*dist2.pdf(testPoints), 'RelTol', 1E-10);
            testCase.verifyEqual(mixture.pdf(testPoints), mixture2.pdf(testPoints), 'RelTol', 1E-10);

            % test integral
            testCase.verifyEqual(mixture.integral(), 1, 'RelTol', 1E-7);
            
            % test angular moments
            testCase.verifyEqual(mixture.trigonometricMoment(0), mixture.trigonometricMomentNumerical(0), 'RelTol', 1E-5)
            testCase.verifyEqual(mixture.trigonometricMoment(1), mixture.trigonometricMomentNumerical(1), 'RelTol', 1E-5)
            testCase.verifyEqual(mixture.trigonometricMoment(2), mixture.trigonometricMomentNumerical(2), 'RelTol', 1E-5)
            testCase.verifyEqual(mixture.trigonometricMoment(3), mixture.trigonometricMomentNumerical(3), 'RelTol', 1E-5)
            
            % test sampling
            n = 10;
            s = mixture.sample(10);
            testCase.verifySize(s,[mixture.dim,n]);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
        end
    end    
   
end
