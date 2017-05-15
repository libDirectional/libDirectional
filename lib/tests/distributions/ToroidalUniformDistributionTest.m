classdef ToroidalUniformDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testToroidalUniformDistribution(testCase)           
            tud = ToroidalUniformDistribution();
            x = repmat([1 2 3 4 5 6], 2, 1);

            %% test pdf
            testCase.verifyEqual(tud.pdf(x), 1/(2*pi)^2*ones(1,size(x,2)));
            
            %% test shift
            tudShifted = tud.shift([1;2]);
            testCase.verifyEqual(tudShifted.pdf(x), 1/(2*pi)^2*ones(1,size(x,2)));

            %% test trigonometric moments
            testCase.verifyEqual(tud.trigonometricMoment(0),tud.trigonometricMomentNumerical(0),'AbsTol', 1E-10);
            testCase.verifyEqual(tud.trigonometricMoment(0),ones(2,1),'RelTol', 1E-10);

            testCase.verifyEqual(tud.trigonometricMoment(1),tud.trigonometricMomentNumerical(1),'AbsTol', 1E-10);
            testCase.verifyEqual(tud.trigonometricMoment(1),zeros(2,1),'RelTol', 1E-10);

            testCase.verifyEqual(tud.trigonometricMoment(2),tud.trigonometricMomentNumerical(2),'AbsTol', 1E-10);
            testCase.verifyEqual(tud.trigonometricMoment(2),zeros(2,1),'RelTol', 1E-10);

            testCase.verifyEqual(tud.trigonometricMoment(3),tud.trigonometricMomentNumerical(3),'AbsTol', 1E-10);
            testCase.verifyEqual(tud.trigonometricMoment(3),zeros(2,1),'RelTol', 1E-10);

            %% test mean
            testCase.verifyWarning(@tud.circularMean,'MEAN:UNDEFINED');

            %% test entropy
            testCase.verifyEqual(tud.entropy(), tud.entropyNumerical(), 'RelTol', 1E-10);

            %% test sampling
            n = 10;
            s = tud.sample(n);
            testCase.verifyEqual(size(s,1), 2);
            testCase.verifyEqual(size(s,2), n);
            testCase.verifyGreaterThanOrEqual(s(:),zeros(size(s(:))));
            testCase.verifyLessThan(s(:),2*pi*ones(size(s(:))));
        end
    end
end