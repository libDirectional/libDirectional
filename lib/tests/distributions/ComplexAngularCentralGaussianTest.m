classdef ComplexAngularCentralGaussianTest < matlab.unittest.TestCase
    properties
        N
        D
        
        cAG
        Z
    end
    
    methods(TestClassSetup)
        function createSampling(tc)
            % Define mean shape (different ideas)            
            C = diag([3 2 1]);
            C = 1/2*(C + C');
            tc.D = 3;
            
            % Define concentration parameter
            tc.N = 5000;
            
            tc.cAG = ComplexAngularCentralGaussian(C);
            
            % Sampling
            tc.Z = tc.cAG.sample(tc.N);        
        end
    end
    
    methods (Test)
        function sanityCheck(tc)
            % verification
            tc.verifyClass(tc.cAG, 'ComplexAngularCentralGaussian');
        end
        
        function testIntegral(tc)
            % test integral
            tc.verifyEqual(tc.cAG.integral(), 1, 'RelTol', 0.1);
        end
        
        function checkSizeOfSampleOutput(tc)
            % verification
            verifySize(tc, tc.Z, [tc.D, tc.N]);
        end
        
        function testFitting(tc)
            cAG_fit = ComplexAngularCentralGaussian.fit(tc.Z);
            
            % verification
            cAG_fit.C = cAG_fit.C/cAG_fit.C(1);
            tc.cAG.C = tc.cAG.C/tc.cAG.C(1);
            csDiv = ComplexBinghamDistribution.CauchySchwarzDivergence( ...
                cAG_fit.C, tc.cAG.C);
            epsilon = 1E-3;
            tc.verifyLessThan(csDiv, epsilon);
        end
    end
end
