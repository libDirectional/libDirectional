classdef SE2BinghamFilterTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testSE2BinghamFilter(testCase)
            filter = SE2BinghamFilter();
            C = [1, -0.1, 3 4;
                -0.1 1 5 6; 
                3 5  -1 0; 
                4 6  0 -1];
            se = SE2BinghamDistribution(C);
            
            %% sanity check            
            filter.setState(se);
            se1 = filter.getEstimate();
            testCase.verifyClass(se1, 'SE2BinghamDistribution');
            testCase.verifyEqual(se1.C, se.C);

            %% predict identity
            filter.setState(se);
            seNoise = SE2BinghamDistribution(-2 * eye(4,4));
            filter.predictIdentity(seNoise);
            seIdentity = filter.getEstimate();
            testCase.verifyClass(seIdentity, 'SE2BinghamDistribution');
            
            %todo properly verify results
             
            %% update identity
            filter.setState(se);
            filter.updateIdentity(se, se.mode());
            seUpdateIdentity = filter.getEstimate();
            testCase.verifyClass(seUpdateIdentity, 'SE2BinghamDistribution');
            
            %todo properly verify results
            
            %% update identity with different measurement
            filter.setState(se);
            z = se.mode()+[0.1,0,3, 4]';
            filter.updateIdentity(seNoise, z/norm(z(1:2)));
            seUpdateIdentity2 = filter.getEstimate();
            testCase.verifyClass(seUpdateIdentity2, 'SE2BinghamDistribution');
            
            %todo properly verify results
        end
    end
end
