classdef QuaternionMultiplicationTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testQuaternionMultiplication(testCase)
            % test multiplication with identity quaternion
            q1 = [1;0;0;0];
            testCase.verifyEqual(q1, quaternionMultiplication(q1,q1), 'AbsTol', 1E-10);
            q2 = [1;1;1;1];
            testCase.verifyEqual(q2, quaternionMultiplication(q1,q2), 'AbsTol', 1E-10);

            % test axioms
            % i^2=j^2=k^2 = -1
            testCase.verifyEqual([-1;0;0;0], quaternionMultiplication([0;1;0;0],[0;1;0;0]), 'AbsTol', 1E-10);
            testCase.verifyEqual([-1;0;0;0], quaternionMultiplication([0;0;1;0],[0;0;1;0]), 'AbsTol', 1E-10);
            testCase.verifyEqual([-1;0;0;0], quaternionMultiplication([0;0;0;1],[0;0;0;1]), 'AbsTol', 1E-10);
            % i*j*k = -1 
            testCase.verifyEqual([-1;0;0;0], quaternionMultiplication([0;1;0;0], quaternionMultiplication([0;0;1;0],[0;0;0;1])), 'AbsTol', 1E-10);
            % i*j=k
            testCase.verifyEqual([0;0;0;1], quaternionMultiplication([0;1;0;0], [0;0;1;0]), 'AbsTol', 1E-10);
            % j*i=-k
            testCase.verifyEqual([0;0;0;-1], quaternionMultiplication([0;0;1;0], [0;1;0;0]), 'AbsTol', 1E-10);

            % test multiplication with inverse
            for a=0.5:0.5:2
                for b=-1:0.5:2
                    c = 7;
                    d = -12;
                    testCase.verifyEqual([1;0;0;0], quaternionMultiplication([a;b;c;d],[a;-b;-c;-d]/(a^2+b^2+c^2+d^2)), 'AbsTol', 1E-10);
                end
            end            
            
        end
    end
end
