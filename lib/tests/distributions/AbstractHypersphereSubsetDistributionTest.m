classdef AbstractHypersphereSubsetDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testMeanAxis(testCase)     
            q = [1,2,3,4]';
            q = q/norm(q);
            M = [quaternionMultiplication(q, [1 0 0 0]'), quaternionMultiplication(q, [0 1 0 0]'), quaternionMultiplication(q, [0 0 1 0]'), quaternionMultiplication(q, [0 0 0 1]')];
            Z = [-10 -2 -1 0]';
            bd = BinghamDistribution(Z,M);
            % Bingham is not properly normalized otherwise
            bd.F = bd.F*bd.integralNumerical();
            

            testCase.verifyEqual(bd.momentNumerical(), bd.moment(), 'AbsTol', 1e-3);
            testCase.verifyEqual(bd.meanAxis(), bd.mode(), 'AbsTol', 1e-10);
            testCase.verifyEqual(bd.meanAxisNumerical(), bd.mode(), 'AbsTol', 1e-7);
        end
    end
end