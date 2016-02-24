classdef BesselratioInverseTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testBesselratio(testCase)
            % calculate bessselratio and invert (small kappa)
            for kappa = [0 0.001 0.01 0.1 1 2 3 10 20 50 100 200 300 400 500 600, 650];
                b0 = besselratio(0,kappa);
                b1 = besselratio(1,kappa);
                b2 = besselratio(2,kappa);
                % v = 0
                testCase.verifyEqual(besselratioInverse(0, b0, 'fisher'), kappa, 'RelTol', 1E-2);
                testCase.verifyEqual(besselratioInverse(0, b0, 'amosfsolve'), kappa, 'RelTol', 1E-10);
                testCase.verifyEqual(besselratioInverse(0, b0, 'matlabfsolve'), kappa, 'RelTol', 1E-10);
                testCase.verifyEqual(besselratioInverse(0, b0, 'sra'), kappa, 'RelTol', 1E-10);
                testCase.verifyEqual(besselratioInverse(0, b0, 'sraamos'), kappa, 'RelTol', 1E-10);
                testCase.verifyEqual(besselratioInverse(0, b0), kappa, 'RelTol', 1E-10);
                % v = 1
                testCase.verifyEqual(besselratioInverse(1, b1, 'amosfsolve'), kappa, 'RelTol', 1E-10);
                testCase.verifyEqual(besselratioInverse(1, b1, 'matlabfsolve'), kappa, 'RelTol', 1E-10);
                testCase.verifyEqual(besselratioInverse(1, b1, 'sra'), kappa, 'RelTol', 1E-10);
                % v = 2
                testCase.verifyEqual(besselratioInverse(2, b2, 'amosfsolve'), kappa, 'RelTol', 1E-10);
                testCase.verifyEqual(besselratioInverse(2, b2, 'matlabfsolve'), kappa, 'RelTol', 1E-10);
                testCase.verifyEqual(besselratioInverse(2, b2, 'sra'), kappa, 'RelTol', 1E-10);
            end
            
            % large kappa
            for kappa = [1E3 1E4 1E5];
                b0 = besselratio(0,kappa);
                b1 = besselratio(1,kappa);
                b2 = besselratio(2,kappa);
                % v = 0
                testCase.verifyEqual(besselratioInverse(0, b0, 'fisher'), kappa, 'RelTol', 1E-2);
                testCase.verifyEqual(besselratioInverse(0, b0, 'sraamos'), kappa, 'RelTol', 1E-10);
                testCase.verifyEqual(besselratioInverse(0, b0), kappa, 'RelTol', 1E-10);
                % v = 1
                testCase.verifyEqual(besselratioInverse(1, b1, 'sra'), kappa, 'RelTol', 1E-10);
                % v = 2
                testCase.verifyEqual(besselratioInverse(2, b2, 'sra'), kappa, 'RelTol', 1E-10);
            end
        end
    end
end
