classdef SE2Test < matlab.unittest.TestCase
   
    methods (Test)
        function testSE2(testCase)
            alpha = 0.7;
            translation = [2;3];
            s = SE2(alpha, translation);
            
            testCase.verifyEqual(s.alpha, alpha);
            testCase.verifyEqual(s.translation, translation);
            
            dq = s.asDualQuaternion();
            testCase.verifySize(dq, [4 1]);
            
            dqFull = s.asDualQuaternion('full');
            testCase.verifySize(dqFull, [8 1]);
            testCase.verifyEqual(dqFull, [dq(1), 0, 0, dq(2), 0, dq(3), dq(4), 0]'); %relation to compact representation
            
            dqMat = s.asDualQuaternionMatrix();
            testCase.verifySize(dqMat, [4 4]);
            testCase.verifyEqual(SE2.dualQuaternionToMatrix(dq), dqMat, 'RelTol', 1E-10);
            testCase.verifyEqual(SE2.matrixToDualQuaternion(dqMat), dq, 'RelTol', 1E-10);
            
            hMat = s.asHomogenousMatrix();
            testCase.verifySize(hMat, [3 3]);
            testCase.verifyEqual(translation, hMat(1:2, 3), 'RelTol', 1E-10);
            
            rotMat = s.rotationPartAsRotationMatrix();
            testCase.verifySize(rotMat, [2 2]);
            testCase.verifyEqual(det(rotMat), 1, 'RelTol', 1E-10);
            testCase.verifyEqual(rotMat*rotMat', eye(2,2), 'AbsTol', 1E-10);
            testCase.verifyEqual(rotMat, hMat(1:2, 1:2), 'RelTol', 1E-10);
        end
        
        function testDualQuaternionMultiply(testCase)
            dq1 = [1,0,0,0]';
            dq2 = [1,0,1,0]';
            
            testCase.verifyEqual(SE2.dualQuaternionMultiply(dq1, dq2), dq2);
            testCase.verifyEqual(SE2.dualQuaternionMultiply(dq2, dq1), dq2);
            
            s1 = SE2(0.3, [2,4]);
            s2 = SE2(0.5, [-2,3]);
            dq = SE2.dualQuaternionMultiply(s1.asDualQuaternion(), s2.asDualQuaternion());
            h = s1.asHomogenousMatrix * s2.asHomogenousMatrix;
            sH = SE2(acos(h(1,1)), h(1:2,3));
            testCase.verifyEqual(dq, sH.asDualQuaternion(), 'RelTol', 1E-10);
        end
       
        function testFromDualQuaternion(testCase)
            se2 = SE2(1,[1;1]);
            se2Converted = SE2.fromDualQuaternion(se2.asDualQuaternion());
            testCase.verifyEqual(se2.alpha,se2Converted.alpha);
            testCase.verifyEqual(se2.translation,se2Converted.translation);
        end
    end
end
