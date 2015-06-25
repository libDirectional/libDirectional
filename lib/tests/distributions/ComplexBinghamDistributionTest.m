classdef ComplexBinghamDistributionTest < matlab.unittest.TestCase
    properties
        N
        D
        
        cB
        Z
    end
    
    methods(TestClassSetup)
        function createSampling(tc)
            %% Define mean shape (different ideas)
            %%
            % Unitary matrix of right eigenvectors
            V = [1 exp(2i*pi*1/3) exp(2i*pi*2/3); 1 1i -1i; 0 0 1]';
            V = ComplexBinghamDistributionTest.gramSchmidtOrthogonalization(V);
            
            % Diagonal matrix of eigenvalues
            Diagonal = diag([50 2 1]);
            
            B = V*Diagonal*V';
            B = 1/2*(B + B');
            tc.D = 3;
            
            %% Define concentration parameter
            tc.N = 5000;
            
            tc.cB = ComplexBinghamDistribution(B);
            
            
            %% Sampling
            tc.Z = tc.cB.sample(tc.N);
            
        end
    end
    methods (Test)
        function sanityCheck(tc)
            %% verification
            tc.verifyClass(tc.cB, 'ComplexBinghamDistribution');
        end
        function testIntegral(tc)
            %% test integral
            tc.verifyEqual(tc.cB.integral(), 1, 'RelTol', 0.1);
            
        end
        function checkSizeOfSampleOutput(tc)
            %% verification
            verifySize(tc, tc.Z, [tc.D, tc.N]);
        end
        function testFitting(tc)
            %%
            
            cB_fit = ComplexBinghamDistribution.fit(tc.Z);
            
            %% verification
            csDiv = ComplexBinghamDistribution.CauchySchwarzDivergence(cB_fit, tc.cB);
            epsilon = 1E-3;
            tc.verifyLessThan(csDiv, epsilon);
        end
        function testCSDivergence(tc)
            
            I = eye(4);
            M_zero = zeros(4);
            
            %% verification
            csDiv = ComplexBinghamDistribution.CauchySchwarzDivergence(I, M_zero);
            tc.verifyEqual(csDiv, 0);
        end
    end
    methods (Static)
        function U = gramSchmidtOrthogonalization(A)
            % Calculate a unitary matrix using Gram Schmidt orthogonalization.
            %
            % :param A: Random real or complex matrix to orthogonalize,
            %   i.e. A = randn(3, 3);
            % :return U: Real orthonormal matrix or complex unitary matrix.
            
            Dim = size(A, 1);
            U = zeros(Dim);
            
            for d = 1:Dim
                v = A(:, d);
                for dd = 1:(d-1)
                    v = v - U(:, dd)' *  A(:, d) * U(:, dd);
                end
                U(:, d) = v / sqrt(v' * v);
            end
            
            % Test with: U'*U = eye(D)
        end
    end
end
