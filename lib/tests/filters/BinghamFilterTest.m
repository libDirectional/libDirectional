classdef BinghamFilterTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testBinghamFilter2d(testCase)
            filter = BinghamFilter();
            Z = [-5 0]';
            phi = 0.4;
            M = [cos(phi), -sin(phi); sin(phi), cos(phi)];
            B = BinghamDistribution(Z,M);
            
            %% sanity check
            filter.setState(B);
            B1 = filter.getEstimate();
            testCase.verifyClass(B1, 'BinghamDistribution');
            testCase.verifyEqual(B.M, B1.M);
            testCase.verifyEqual(B.Z, B1.Z);
            
            %% predict identity
            filter.setState(B);
            Bnoise = BinghamDistribution([-3 0]', [0 1; 1 0]);
            filter.predictIdentity(Bnoise);
            B2 = filter.getEstimate();
            testCase.verifyClass(B2, 'BinghamDistribution');
            testCase.verifyEqual(abs(B.M), abs(B2.M), 'RelTol', 1E-10); %each column of M is only determined up to sign
            testCase.verifyGreaterThanOrEqual(B2.Z, B.Z);
            
            %% update identity
            filter.setState(B);
            filter.updateIdentity(Bnoise, B.mode());
            B3 = filter.getEstimate();
            testCase.verifyClass(B3, 'BinghamDistribution');
            testCase.verifyEqual(abs(B.mode()), abs(B3.mode())); %antipodal symmetry
            testCase.verifyLessThanOrEqual(B3.Z, B.Z);
            
            %% update identity with different measurement
            filter.setState(B);
            z = B.mode()+[0.1,0]';
            filter.updateIdentity(Bnoise, z/norm(z));
            B4 = filter.getEstimate();
            testCase.verifyClass(B4, 'BinghamDistribution');
            testCase.verifyLessThanOrEqual(B4.Z, B.Z);
            
            %% predict nonlinear
            filter.setState(B);
            Bnoise = BinghamDistribution([-3 0]', [0 1; 1 0]);
            filter.predictNonlinear(@(x) x, Bnoise);
            B5 = filter.getEstimate();
            testCase.verifyClass(B5, 'BinghamDistribution');
            testCase.verifyEqual(abs(B.M), abs(B5.M), 'RelTol', 1E-10); %each column of M is only determined up to sign
            testCase.verifyGreaterThanOrEqual(B5.Z, B.Z);
        end
        
        function testBinghamFilter4d(testCase)
            filter = BinghamFilter();
            Z = [-5 -3 -2 0]';
            M = [0 1 0 0;
                 0 0 0 1;
                 0 0 1 0;
                 1 0 0 0];
            B = BinghamDistribution(Z,M);
            
            %% sanity check
            filter.setState(B);
            B1 = filter.getEstimate();
            testCase.verifyClass(B1, 'BinghamDistribution');
            testCase.verifyEqual(B.M, B1.M);
            testCase.verifyEqual(B.Z, B1.Z);
            
            %% predict identity
            filter.setState(B);
            Bnoise = BinghamDistribution([-2 -2 -2 0]', [0 0 0 1; 0 0 1 0; 0 1 0 0; 1 0 0 0 ]);
            filter.predictIdentity(Bnoise);
            BpredictIdentity = filter.getEstimate();
            testCase.verifyClass(BpredictIdentity, 'BinghamDistribution');
            testCase.verifyEqual(B.M, BpredictIdentity.M); 
            testCase.verifyGreaterThanOrEqual(BpredictIdentity.Z, B.Z);
            
            %% predict nonlinear
            filter.setState(B);
            filter.predictNonlinear(@(x) x, Bnoise);
            BpredictNonlinearIdentity = filter.getEstimate();
            testCase.verifyClass(BpredictNonlinearIdentity, 'BinghamDistribution');
            testCase.verifyEqual(abs(B.M), abs(BpredictNonlinearIdentity.M), 'AbsTol', 1E-10); %each column of M is only determined up to sign
            testCase.verifyGreaterThanOrEqual(BpredictNonlinearIdentity.Z, B.Z);
            testCase.verifyEqual(abs(BpredictIdentity.M), abs(BpredictNonlinearIdentity.M), 'AbsTol', 1E-10); %each column of M is only determined up to sign
            testCase.verifyEqual(BpredictIdentity.Z, BpredictNonlinearIdentity.Z, 'RelTol', 0.15);
            
            %% update identity
            filter.setState(B);
            filter.updateIdentity(Bnoise, B.mode());
            B3 = filter.getEstimate();
            testCase.verifyClass(B3, 'BinghamDistribution');
            testCase.verifyEqual(B.mode(), B3.mode());
            testCase.verifyLessThanOrEqual(B3.Z, B.Z);
            
            %% update identity with different measurement
            filter.setState(B);
            z = B.mode()+[0.1,0.1,0,0]';
            filter.updateIdentity(Bnoise, z/norm(z));
            B4 = filter.getEstimate();
            testCase.verifyClass(B4, 'BinghamDistribution');
            testCase.verifyLessThanOrEqual(B4.Z, B.Z);
        end
    end
end
