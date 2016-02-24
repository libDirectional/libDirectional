classdef BesselratioTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testBesselratio(testCase)
            %% compare to quotient of besseli for kappa=0
            testCase.verifyEqual(besselratio(0,0), 0, 'RelTol', 1E-10);
            testCase.verifyEqual(besselratio(1,0), 0, 'RelTol', 1E-10);
            testCase.verifyEqual(besselratio(2,0), 0, 'RelTol', 1E-10);
                
            %% compare to quotient of besseli for small kappa
            for kappa = [0.001 0.01 0.1 1 2 3 10 20 50 100 200 300 400 500 600, 650];
                testCase.verifyEqual(besselratio(0,kappa), besseli(1,kappa)/besseli(0,kappa), 'RelTol', 1E-10);
                testCase.verifyEqual(besselratio(1,kappa), besseli(2,kappa)/besseli(1,kappa), 'RelTol', 1E-10);
                testCase.verifyEqual(besselratio(2,kappa), besseli(3,kappa)/besseli(2,kappa), 'RelTol', 1E-10);
                testCase.verifyEqual(besselratio(0,kappa)*besselratio(1,kappa), besseli(2,kappa)/besseli(0,kappa), 'RelTol', 1E-10);
            end
            
            %% check that besselratio is strictly increasing even for large kappa and compare to scaled MATLAB bessel functions
            for kappa = [1E3, 1E4, 1E5, 1E6, 1E7]
                testCase.verifyGreaterThan(besselratio(0,kappa), besselratio(0,0.9*kappa));
                testCase.verifyEqual(besselratio(0,kappa), besseli(1,kappa,1)/besseli(0,kappa,1), 'RelTol', 1E-10);
                testCase.verifyEqual(besselratio(1,kappa), besseli(2,kappa,1)/besseli(1,kappa,1), 'RelTol', 1E-10);
                testCase.verifyEqual(besselratio(2,kappa), besseli(3,kappa,1)/besseli(2,kappa,1), 'RelTol', 1E-10);
            end
            
            %% check lim kappa->infinity
            testCase.verifyEqual(besselratio(0,1E100), 1);
            testCase.verifyEqual(besselratio(1,1E100), 1);
            testCase.verifyEqual(besselratio(2,1E100), 1);
        end
    end
end
