classdef SE2BinghamDistributionTest < matlab.unittest.TestCase
    
    properties
    end
    
    methods (Test)
        
        function testSE2BinghamDistribution (testCase)
            C1 = -diag([2 3]);
            C2 = [0.1 0.2; 0.01 0.3];
            C3 = -diag([1 2]);
            
            C = [C1 C2'; C2 C3];
            
            % test constructor
            seA = SE2BinghamDistribution(C);
            seB = SE2BinghamDistribution(C1,C2,C3);
            
            testCase.verifyEqual(seA.C, seB.C);
            testCase.verifyEqual(seA.C1, seB.C1);
            testCase.verifyEqual(seA.C2, seB.C2);
            testCase.verifyEqual(seA.C3, seB.C3);
            
            %% computeCovarianceMCMC
            rng default
            CovMatrix = seA.computeCovarianceMCMC();
            testCase.verifySize(CovMatrix, [4,4]); % size
            testCase.verifyGreaterThan(eig(CovMatrix), [0,0,0,0]'); % pos. definite
            
            %% mode
            m = seA.mode();
            testCase.verifySize(m, [4,1]); % size
            
            %% sampling
            n = 100;
            s = seA.sample(n);
            testCase.verifySize(s, [4 n]); % size
            testCase.verifyEqual(s(1,:).^2+s(2,:).^2, ones(1,n), 'RelTol', 1E-10); % normalization of first two entries
        end
        
        function testNormalizationConstant (testCase)
            % Tests whether the pdf integrates to the normalization constant.
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = false; end            
            
            if enableExpensive
                C1 = -diag([2 3]);
                C2 = [0.1 0.2; 0.01 0.3];
                C3 = -diag([1 2]);
                C = [C1 C2'; C2 C3];

                % Unnormalized PDF
                unPDF = @(al,x,y) ...
                    exp([sin(al); cos(al); x; y]' * C * [sin(al); cos(al); x; y]);
                pdfIntegrand =  @(al,x,y) arrayfun(unPDF, al, x, y); % As Integrand

                % Computation of reference value.
                refVal = integral3(pdfIntegrand, -pi,pi, -inf,inf,-inf,inf);

                % Computation of the actual value.
                s = SE2BinghamDistribution(C1, C2, C3);

                testCase.verifyEqual(refVal,s.NC, 'RelTol',1E-10);
            end
        end
        
        function testDeterministicSampling (testCase)
            % Tests deterministic sampling.
            % This test checks whether the deterministic sampling algorithm is
            % consistent with the parameter estimation algorithm.
                       
            % Set distribution parameters
            C1 = -diag([2 3]);
            C2 = [0.1 0.2; 0.01 0.3];
            C3 = -diag([1 2]);
            
            % Generate distribution
            s = SE2BinghamDistribution(C1, C2, C3);
            
            % Generate deterministic samples
            [samples, weights] = s.sampleDeterministic();
            testCase.verifyEqual(sum(weights), 1, 'RelTol', 1E-10); % sum of weights
            testCase.verifyEqual(weights, ones(size(weights))/length(weights), 'RelTol', 1E-10); %unform weights
            
            % Carry out fitting.
            sFitted = SE2BinghamDistribution.fit(samples);
            
            % Some interesting Debug matrices
            % Estimated Bingham Params: sFitted.C1 - sFitted.C2' * pinv(sFitted.C3) * sFitted.C2
            % Bingham Params: C1 - C2' * pinv(C3) * C2
            % Mean shift: -pinv(C3)*C2
            
            %C1 only matchs up to a constrant times the identity matrix
            C1normalized = C1-C1(2,2)*eye(2,2);
            C1fittedNormalized = sFitted.C1 - sFitted.C1(2,2)*eye(2,2);
            testCase.verifyEqual(C1normalized, C1fittedNormalized, 'AbsTol', 1E-10);
            
            testCase.verifyEqual(C2,sFitted.C2, 'AbsTol',1E-10);
            testCase.verifyEqual(C3,sFitted.C3, 'AbsTol',1E-10);            
        end        
    end
end

