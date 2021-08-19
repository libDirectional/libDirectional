classdef HypertoroidalUniformDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testHypertoroidalUniformDistribution(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true; end
            if enableExpensive        
                dims = 1:3;
            else
                dims = 1:2;
            end
            
            for dim = dims
                hud = HypertoroidalUniformDistribution(dim);
                x = repmat([1 2 3 4 5 6], dim, 1);

                % test pdf
                testCase.verifyEqual(hud.pdf(x), 1/(2*pi)^dim*ones(1,size(x,2)));

                % test trigonometric moments
                testCase.verifyEqual(hud.trigonometricMoment(0),hud.trigonometricMomentNumerical(0),'AbsTol', 1E-10);
                testCase.verifyEqual(hud.trigonometricMoment(0),ones(dim,1),'RelTol', 1E-10);

                testCase.verifyEqual(hud.trigonometricMoment(1),hud.trigonometricMomentNumerical(1),'AbsTol', 1E-10);
                testCase.verifyEqual(hud.trigonometricMoment(1),zeros(dim,1),'RelTol', 1E-10);

                testCase.verifyEqual(hud.trigonometricMoment(2),hud.trigonometricMomentNumerical(2),'AbsTol', 1E-10);
                testCase.verifyEqual(hud.trigonometricMoment(2),zeros(dim,1),'RelTol', 1E-10);

                testCase.verifyEqual(hud.trigonometricMoment(3),hud.trigonometricMomentNumerical(3),'AbsTol', 1E-10);
                testCase.verifyEqual(hud.trigonometricMoment(3),zeros(dim,1),'RelTol', 1E-10);

                % test mean
                testCase.verifyWarning(@hud.circularMean,'MEAN:UNDEFINED');

                % test entropy
                testCase.verifyEqual(hud.entropy(), hud.entropyNumerical(), 'RelTol', 1E-10);

                % test sampling
                n = 10;
                s = hud.sample(n);
                testCase.verifyEqual(size(s,1), dim);
                testCase.verifyEqual(size(s,2), n);
                testCase.verifyGreaterThanOrEqual(s(:),zeros(size(s(:))));
                testCase.verifyLessThan(s(:),2*pi*ones(size(s(:))));
                
                % test integral 
                testCase.verifyEqual(hud.integral(), 1,'RelTol', 1E-10);
                testCase.verifyEqual(hud.integral(zeros(dim,1), 2*pi*ones(dim,1)), 1,'RelTol', 1E-10);
                testCase.verifyEqual(hud.integral(zeros(dim,1), 6*pi*ones(dim,1)), 3^dim,'RelTol', 1E-10);
                testCase.verifyEqual(hud.integral(2*pi*ones(dim,1), zeros(dim,1)), (-1)^dim,'RelTol', 1E-10);
                testCase.verifyEqual(hud.integral((1:dim)', ((1:dim).^2)'),hud.integralNumerical((1:dim)', ((1:dim).^2)'),'RelTol', 1E-10);
            end
        end
    end
end