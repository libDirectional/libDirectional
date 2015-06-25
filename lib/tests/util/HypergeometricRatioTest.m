classdef HypergeometricRatioTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testHypergeometricRatioAgainstInvUpsilon(testCase)
            %%
            concentrationMax = 1000;
            D = 3;
            x = [0 10 20 50 100 500];
            r = hypergeometricRatio(D, x);
            x_hat = hypergeometricRatioInverse(r, D, concentrationMax);
            
            %%
            testCase.verifyEqual(x, x_hat, 'AbsTol', 1e-4);
        end
    end
end
