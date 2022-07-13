classdef HyperhemisphericalBinghamDistributionTest < matlab.unittest.TestCase
    methods (Test)
        function testHyperhemisphericalBinghamDistribution(testCase)
            for i=1:4
                switch i
                    case 1
                        M = eye(2,2);
                        Z = [-3 0]';
                    case 2
                        phi = 0.7;
                        M = [cos(phi),-sin(phi);sin(phi),cos(phi)];
                        Z = [-5 0]';
                    case 3
                        M = eye(4,4);
                        Z = [-10 -2 -1 0]';
                    case 4
                        q = [1,2,3,4]';
                        q = q/norm(q);
                        M = [quaternionMultiplication(q, [1 0 0 0]'), quaternionMultiplication(q, [0 1 0 0]'), quaternionMultiplication(q, [0 0 1 0]'), quaternionMultiplication(q, [0 0 0 1]')];
                        Z = [-10 -2 -1 0]';
                end

                B = HyperhemisphericalBinghamDistribution(Z,M);
                
                %% test pdf
                rng default
                testpoints = rand(B.dim, 20);
                testpoints = bsxfun(@rdivide,testpoints, sum(testpoints));
                for j=1:size(testpoints,2)
                    testCase.verifyEqual(B.pdf(testpoints(:,j)), 2/B.F*exp(testpoints(:,j)'*M*diag(Z)*M'*testpoints(:,j)), 'RelTol', 1E-10);
                end

                %% sanity check
                testCase.verifyClass(B, 'HyperhemisphericalBinghamDistribution');
                testCase.verifyEqual(B.M, M);
                testCase.verifyEqual(B.Z, Z);
                testCase.verifyEqual(B.dim, length(Z));

                %% test mode
                % respect antipodal symmetry
                modeNumerical = B.modeNumerical;
                difference = min(abs(modeNumerical-B.mode), abs(modeNumerical+B.mode)); 
                testCase.verifyEqual(difference, zeros(B.dim,1), 'AbsTol',1E-5)

                %% test integral
                if B.dim < 4
                    testCase.verifyEqual(B.integral(), 1, 'RelTol', 1E-1);
                end

                %% test multiplication
                Bmul = B.multiply(B);
                if isequal(M*M,M)
                    testCase.verifyEqual(Bmul.M, B.M, 'RelTol', 1E-10); %M is the same for both
                end
                renormconst = Bmul.pdf(B.mode()) / B.pdf(B.mode)^2;
                for j=1:size(testpoints,2)
                    testCase.verifyEqual(Bmul.pdf(testpoints(:,j)), renormconst*B.pdf(testpoints(:,j))^2, 'RelTol', 1E-10);
                end

                %% test composition
                Bcomp = B.compose(B);
                if length(M)==4
                    newMode = quaternionMultiplication(B.mode(), B.mode());
                else
                    newMode = complexMultiplication(B.mode(), B.mode());
                end
                if i~=4
                    %todo fix this for i=4
                    testCase.verifyEqual(abs(Bcomp.mode()),abs(newMode), 'RelTol', 1E-10);
                end
                testCase.verifyGreaterThanOrEqual(Bcomp.Z,B.Z); 
            end
        end
    end
end

