classdef BinghamDistributionTest < matlab.unittest.TestCase
    % Contains tests for the Bingham Distribution
    
    properties
    end
    
    methods (Test)
        function testBinghamDistribution(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true; end
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

                B = BinghamDistribution(Z,M);
                
                %% test pdf
                rng default
                testpoints = rand(B.dim, 20);
                testpoints = bsxfun(@rdivide,testpoints, sum(testpoints));
                for j=1:size(testpoints,2)
                    testCase.verifyEqual(B.pdf(testpoints(:,j)), 1/B.F*exp(testpoints(:,j)'*M*diag(Z)*M'*testpoints(:,j)), 'RelTol', 1E-10);
                end

                %% sanity check
                testCase.verifyClass(B, 'BinghamDistribution');
                testCase.verifyEqual(B.M, M);
                testCase.verifyEqual(B.Z, Z);
                testCase.verifyEqual(B.dim, length(Z));

                %% test mode
                % respect antipodal symmetry
                modeNumerical = B.modeNumerical;
                difference = min(abs(modeNumerical-B.mode), abs(modeNumerical+B.mode)); 
                testCase.verifyEqual(difference, zeros(B.dim,1), 'AbsTol',1E-5)

                %% test integral
                if B.dim < 4 || enableExpensive
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

                %% test deterministic sampling
                % This tests seem to weak which is mostly due to use of
                % approximations.
                lambda = 0.5;
                [samples, weights] = B.sampleDeterministic(lambda);
                Bfitted = BinghamDistribution.fit(samples,weights);
                testCase.verifyEqual(Bfitted.Z, B.Z, 'RelTol', 5*1E-2);
                testCase.verifyEqual(abs(Bfitted.M), abs(B.M), 'AbsTol', 1E-1);
                % check that covariance is retained
                sampleCov = samples*diag(weights)*samples';
                testCase.verifyEqual(sampleCov, B.moment(), 'AbsTol', 0.01); % this is kind of inaccurate, maybe because of the saddlepoint approx?

                %% test stochastic sampling
                rng default
                n = 1000;
                samples = B.sampleGlover(n);
                testCase.verifySize(samples, [B.dim, n]);
                testCase.verifyEqual(sum(samples.*samples), ones(1,n), 'RelTol', 1E-10);                
                Bfitted = BinghamDistribution.fit(samples);
                testCase.verifyEqual(Bfitted.Z, B.Z, 'RelTol', 0.2); % this can be quite imprecise
                testCase.verifyEqual(abs(Bfitted.M), abs(B.M), 'AbsTol', 0.2); % this can be quite imprecise
                % the following test fails with very low probability for
                % Glover's sampling because rejected samples are repeated
                testCase.verifyEqual(size(unique(samples', 'rows'),1), n); % make sure all samples are unique
                
                samples = B.sampleKent(n);
                testCase.verifySize(samples, [B.dim, n]);
                testCase.verifyEqual(sum(samples.*samples), ones(1,n), 'RelTol', 1E-10);                
                Bfitted = BinghamDistribution.fit(samples);
                testCase.verifyEqual(Bfitted.Z, B.Z, 'RelTol', 0.2); % this can be quite imprecise
                testCase.verifyEqual(abs(Bfitted.M), abs(B.M), 'AbsTol', 0.2); % this can be quite imprecise
                testCase.verifyEqual(size(unique(samples', 'rows'),1), n); % make sure all samples are unique

                %% test scatter/covariance matrix
                S = B.moment();
                % compare to Glover's method for calculating the scatter
                % matrix
                S2 = zeros(B.dim, B.dim);
                for j=1:B.dim
                    sigma = B.dF(j)/B.F;
                    v = B.M(:,j);
                    S2 = S2 + sigma*(v*v');
                end
                testCase.verifyEqual(S,S2,'RelTol', 1E-2);
                testCase.verifyEqual(trace(S),1,'RelTol', 1E-1); % this is quite imprecise
                BfromS = BinghamDistribution.fitToMoment(S);
                testCase.verifyEqual(BfromS.Z, B.Z, 'RelTol', 1E-1);
                testCase.verifyEqual(abs(BfromS.M), abs(B.M), 'AbsTol', 1E-10);
                
                %% test sampleWeighted
                n = 10;
                [s,w] = B.sampleWeighted(n);
                testCase.verifyEqual(size(s,1),B.dim);
                testCase.verifyEqual(size(s,2),n);
                testCase.verifyEqual(size(w,1),1);
                testCase.verifyEqual(size(w,2),n);
                testCase.verifyEqual(sum(s.^2,1), ones(1,n), 'RelTol', 1E-10);
                testCase.verifyEqual(sum(w),1, 'RelTol', 1E-10);
                
                %% test bounds for d=2
                if B.dim == 2
                    for p = 0.1:0.1:0.9
                        alpha = B.bounds(p);
                        f = @(phi) B.pdf([cos(phi); sin(phi)]);
                        mB = B.mode();
                        mBangle = mod(atan2(mB(2), mB(1)),2*pi);
                        testCase.verifyEqual(integral(f, mBangle-alpha, mBangle+alpha), p, 'RelTol', 1E-3)
                    end
                end
                
                %% test covariance
                C = B.gaussianCovariance();
                testCase.verifySize(C, [B.dim B.dim]);
                testCase.verifyEqual(C,C');
                testCase.verifyGreaterThanOrEqual(eig(C), zeros(B.dim,1));
            end
        end
        
        function testVmRelation(testCase)
            phi = 1.9;
            M = [cos(phi),-sin(phi);sin(phi),cos(phi)]';
            Z = [-8 0]';
            B = BinghamDistribution(Z,M);
            
            % convert to VM and compare pdf
            m=B.mode;
            vm = VMDistribution(2*atan2(m(2),m(1)), -Z(1)/2);
            x = 0:0.2:2*pi;
            testCase.verifyEqual(vm.pdf(x), B.pdf([cos(x/2);sin(x/2)]),'RelTol',1E-10);
            
            % compare deterministic samples
            wd = vm.toDirac3();
            samples = B.sampleDeterministic('uniform');
            samples1D = mod(atan2(samples(2,:), samples(1,:))*2,2*pi);
            testCase.verifyEqual(sort(wd.d), sort(samples1D), 'RelTol',1E-10);
        end
        
        function testNormConst(testCase)
            % test normalization constant and derivatives
            
            %% 4d, compare to mhg, only use small values
            modes4d = {'default','saddlepoint','mhg', 'wood', 'glover'};
            compF = @(Z) BinghamDistribution.computeUnitSphereSurface(length(Z))* mhg(100, 2, 0.5, length(Z)/2, Z);
            for i=1:length(modes4d)
                mode = modes4d{i};
                Z = [-4 -3 -1 0]';
                testCase.verifyEqual(BinghamDistribution.computeF(Z,mode), compF(Z), 'RelTol', 1E-1);
                Z = [-10 -3 -1 0]';
                testCase.verifyEqual(BinghamDistribution.computeF(Z,mode), compF(Z), 'RelTol', 1E-1);
                Z = [-11 -11 -11 0]';
                testCase.verifyEqual(BinghamDistribution.computeF(Z,mode), compF(Z), 'RelTol', 1E-1);
                if ~strcmp(mode, 'glover') %there is a bug when entries are zero
                    Z = [-21 -20 0 0]';
                    testCase.verifyEqual(BinghamDistribution.computeF(Z,mode), compF(Z), 'RelTol', 1E-1);
                    Z = [-23 0 0 0]';
                    testCase.verifyEqual(BinghamDistribution.computeF(Z,mode), compF(Z), 'RelTol', 1E-1);
                end
            end
            
            %% 3d, compare to mhg, only use small values
            modes3d = {'default','saddlepoint','mhg','glover'};
            compF = @(Z) BinghamDistribution.computeUnitSphereSurface(length(Z))* mhg(100, 2, 0.5, length(Z)/2, Z);
            for i=1:length(modes3d)
                mode = modes3d{i};
                Z = [-3 -2 0]';
                testCase.verifyEqual(BinghamDistribution.computeF(Z,mode), compF(Z), 'RelTol', 1E-1);
                Z = [-10 -2 0]';
                testCase.verifyEqual(BinghamDistribution.computeF(Z,mode), compF(Z), 'RelTol', 1E-1);
                Z = [-11 -11 0]';
                testCase.verifyEqual(BinghamDistribution.computeF(Z,mode), compF(Z), 'RelTol', 1E-1);
                Z = [-21 -20 0]';
                testCase.verifyEqual(BinghamDistribution.computeF(Z,mode), compF(Z), 'RelTol', 1E-1);
                if ~strcmp(mode, 'glover') %there is a bug when entries are zero
                    Z = [-23 0 0]';
                    testCase.verifyEqual(BinghamDistribution.computeF(Z,mode), compF(Z), 'RelTol', 1E-1);
                end
            end
            
            %% 2d
            modes2d = {'default','bessel','hypergeom','saddlepoint','mhg','glover'};
            for i=1:length(modes2d)
                mode = modes2d{i};
                Z = [-1 0]';
                testCase.verifyEqual(BinghamDistribution.computeF(Z,mode), compF(Z), 'RelTol', 1E-1);
                Z = [-10 0]';
                testCase.verifyEqual(BinghamDistribution.computeF(Z,mode), compF(Z), 'RelTol', 1E-1);
                Z = [-20 0]';
                testCase.verifyEqual(BinghamDistribution.computeF(Z,mode), compF(Z), 'RelTol', 1E-1);
            end
            
            %% derivatives
            epsilon = 1E-5;
            compDz = @(dim,totalDim) epsilon*[zeros(dim-1,1); 1; zeros(totalDim-dim,1)];
            compDiffDim = @(Z,dim,totalDim) (compF(Z+compDz(dim,totalDim))-compF(Z-compDz(dim,totalDim)))/(2*epsilon);
            compDiff4 = @(Z) [compDiffDim(Z,1,4),compDiffDim(Z,2,4),compDiffDim(Z,3,4),compDiffDim(Z,4,4)];
            compDiff2 = @(Z) [compDiffDim(Z,1,2),compDiffDim(Z,2,2)];
            
            %% derivatives 4d
            modes4dDiff = {'default','saddlepoint','finitedifferences-default', 'finitedifferences-saddlepoint', 'finitedifferences-wood', 'finitedifferences-mhg'};
            % finitedifferences-glover does not work at this time
            for i=1:length(modes4dDiff)
                mode = modes4dDiff{i};
                Z = [-4 -3 -1 0]';
                testCase.verifyEqual(BinghamDistribution.computeDF(Z,mode), compDiff4(Z), 'RelTol', 1E-1);
                Z = [-10 -3 -1 0]';
                testCase.verifyEqual(BinghamDistribution.computeDF(Z,mode), compDiff4(Z), 'RelTol', 1E-1);
                Z = [-11 -11 -11 0]';
                testCase.verifyEqual(BinghamDistribution.computeDF(Z,mode), compDiff4(Z), 'RelTol', 1E-1);
                % todo: problems with 0 entries in Z?!
                %Z = [-21 -20 0 0]';
                %testCase.verifyEqual(BinghamDistribution.computeDF(Z,mode), compDiff4(Z), 'RelTol', 1E-1);
                %Z = [-23 0 0 0]';
                %testCase.verifyEqual(BinghamDistribution.computeDF(Z,mode), compDiff4(Z), 'RelTol', 1E-1);
            end
            
            %% derivatives 2d
            modes2dDiff = {'default','bessel','hypergeom','saddlepoint','finitedifferences-default','finitedifferences-bessel','finitedifferences-hypergeom','finitedifferences-mhg'};
            % finitedifferences-saddlepoint works but is very inaccurate
            % finitedifferences-glover does not work at this time
            for i=1:length(modes2dDiff)
                mode = modes2dDiff{i};
                Z = [-1 0]';
                testCase.verifyEqual(BinghamDistribution.computeDF(Z,mode), compDiff2(Z), 'RelTol', 1E-1);
                Z = [-10 0]';
                testCase.verifyEqual(BinghamDistribution.computeDF(Z,mode), compDiff2(Z), 'RelTol', 1E-1);
                Z = [-20 0]';
                testCase.verifyEqual(BinghamDistribution.computeDF(Z,mode), compDiff2(Z), 'RelTol', 1E-1);
            end
        end

        function testParameterEstimation(testCase)
            % Obtain random orthogonal matrices for M and and random
            % vectors for Z, perform deterministic sampling and then
            % parameter estimation.
            for d = 2:5 % dimension
                rng default;
                for i = 1:20
                    B = BinghamDistributionTest.getRandomBingham(d);
                    
                    [s,w] = B.sampleDeterministic();
                    Bhemi = BinghamDistribution.fit(s,w); % samples on one hemisphere only
                    Bmirrored = BinghamDistribution.fit([s -s], [w w]/2); % version with mirrored samples
                    fitOptions.algorithm = 'fsolve';
                    Bgaussnewton = BinghamDistribution.fit(s,w, fitOptions); % samples on one hemisphere only
                    Bscatter = BinghamDistribution.fitToMoment(B.moment());
                    
                    % All values here are to weak. The saddlepoint
                    % approximation seems to be particularly poor in low
                    % dimensions.
                    testCase.verifyEqual(B.Z, Bhemi.Z, 'RelTol', 2*1E-1);
                    testCase.verifyEqual(abs(B.M), abs(Bhemi.M), 'RelTol', 1E-10); %too weak test, all columns should be equal up to sign
                    testCase.verifyEqual(B.Z, Bmirrored.Z, 'RelTol', 2*1E-1);
                    testCase.verifyEqual(abs(B.M), abs(Bmirrored.M), 'RelTol', 1E-10); %too weak test, all columns should be equal up to sign
                    testCase.verifyEqual(B.Z, Bgaussnewton.Z, 'RelTol', 0.15);
                    testCase.verifyEqual(abs(B.M), abs(Bgaussnewton.M), 'RelTol', 1E-10); %too weak test, all columns should be equal up to sign
                    testCase.verifyEqual(B.Z, Bscatter.Z, 'RelTol', 2*1E-1);
                    testCase.verifyEqual(abs(B.M), abs(Bscatter.M), 'RelTol', 1E-10); %too weak test, all columns should be equal up to sign
                    
                end
            end
        end
        
        function testFitToMoment(testCase)
            options.algorithm = 'fminunc';
            options.fMethod = 'wood';
            options.dFmethod = 'finitedifferences-wood';
            
            S = diag([0.2622, 0.2606, 0.2417, 0.2354]);
            B = BinghamDistribution.fitToMoment(S, options);
            S2 = B.moment();
            testCase.verifyEqual(S, S2, 'RelTol', 1E-3);
            testCase.verifyEqual(trace(S2), 1, 'RelTol', 1E-10);
            
            % case where trace==1 does not hold for the original matrix
            S3 = diag([0.2622, 0.3140, 0.2417, 0.2354]);
            B = BinghamDistribution.fitToMoment(S3, options);
            S4 = B.moment();
            testCase.verifyEqual(S3/sum(diag(S3)), S4, 'RelTol', 2E-4);
            testCase.verifyEqual(trace(S4), 1, 'RelTol', 1E-10);
        end
        
        function testComposition(testCase)
            function S = composeGlover(S1, S2) %implementation from Glover's libbingham, MATLAB version
                a11 = S1(1,1);
                a12 = S1(1,2);
                a13 = S1(1,3);
                a14 = S1(1,4);
                a22 = S1(2,2);
                a23 = S1(2,3);
                a24 = S1(2,4);
                a33 = S1(3,3);
                a34 = S1(3,4);
                a44 = S1(4,4);

                b11 = S2(1,1);
                b12 = S2(1,2);
                b13 = S2(1,3);
                b14 = S2(1,4);
                b22 = S2(2,2);
                b23 = S2(2,3);
                b24 = S2(2,4);
                b33 = S2(3,3);
                b34 = S2(3,4);
                b44 = S2(4,4);

                S = zeros(4,4);
                S(1,1) = a11*b11 - 2*a12*b12 - 2*a13*b13 - 2*a14*b14 + a22*b22 + 2*a23*b23 + 2*a24*b24 + a33*b33 + 2*a34*b34 + a44*b44;
                S(1,2) = a11*b12 + a12*b11 + a13*b14 - a14*b13 - a12*b22 - a22*b12 - a13*b23 - a23*b13 - a14*b24 - a24*b14 - a23*b24 + a24*b23 - a33*b34 + a34*b33 - a34*b44 + a44*b34;
                S(1,3) = a11*b13 + a13*b11 - a12*b14 + a14*b12 - a12*b23 - a23*b12 - a13*b33 + a22*b24 - a24*b22 - a33*b13 - a14*b34 - a34*b14 + a23*b34 - a34*b23 + a24*b44 - a44*b24;
                S(1,4) = a11*b14 + a12*b13 - a13*b12 + a14*b11 - a12*b24 - a24*b12 - a22*b23 + a23*b22 - a13*b34 - a34*b13 - a23*b33 + a33*b23 - a14*b44 - a24*b34 + a34*b24 - a44*b14;
                S(2,2) = 2*a12*b12 + a11*b22 + a22*b11 + 2*a13*b24 - 2*a14*b23 + 2*a23*b14 - 2*a24*b13 - 2*a34*b34 + a33*b44 + a44*b33;
                S(2,3) = a12*b13 + a13*b12 + a11*b23 + a23*b11 - a12*b24 + a14*b22 - a22*b14 + a24*b12 + a13*b34 - a14*b33 + a33*b14 - a34*b13 + a24*b34 + a34*b24 - a23*b44 - a44*b23;
                S(2,4) = a12*b14 + a14*b12 + a11*b24 + a12*b23 - a13*b22 + a22*b13 - a23*b12 + a24*b11 - a14*b34 + a34*b14 + a13*b44 + a23*b34 - a24*b33 - a33*b24 + a34*b23 - a44*b13;
                S(3,3) = 2*a13*b13 + 2*a14*b23 - 2*a23*b14 + a11*b33 + a33*b11 - 2*a12*b34 + 2*a34*b12 - 2*a24*b24 + a22*b44 + a44*b22;
                S(3,4) = a13*b14 + a14*b13 - a13*b23 + a23*b13 + a14*b24 - a24*b14 + a11*b34 + a12*b33 - a33*b12 + a34*b11 + a23*b24 + a24*b23 - a12*b44 - a22*b34 - a34*b22 + a44*b12;
                S(4,4) = 2*a14*b14 - 2*a13*b24 + 2*a24*b13 + 2*a12*b34 - 2*a23*b23 - 2*a34*b12 + a11*b44 + a22*b33 + a33*b22 + a44*b11;
                S(2,1) = S(1,2);
                S(3,1) = S(1,3);
                S(4,1) = S(1,4);
                S(3,2) = S(2,3);
                S(4,2) = S(2,4);
                S(4,3) = S(3,4);
            end
            
            d = 4;
            rng default
            for i = 1:30
                % generate random Bingham distributions
                B1 = BinghamDistributionTest.getRandomBingham(d);
                B2 = BinghamDistributionTest.getRandomBingham(d);

                % Compare scatter Matrix after composition with Glover's
                % scatter matrix (this also tests the MLE performed in
                % between).
                S = B1.compose(B2).moment();
                Sglover = composeGlover(B1.moment(), B2.moment());
                testCase.verifyEqual(S,Sglover, 'AbsTol', 0.15);
            end
        end
        
        function testQuantization(testCase)
            rng default
            B = BinghamDistributionTest.getRandomBingham(2);
            n = 11;
            [s,w] = B.sampleOptimalQuantization(n);
            testCase.verifyEqual(sum(w), 1, 'RelTol', 1E-10);
            testCase.verifySize(s, [1 2*n]);
            testCase.verifySize(w, [1 2*n]);
            wd = WDDistribution(s,w);
            testCase.verifyEqual(wd.trigonometricMoment(1), 0, 'AbsTol', 1E-10);
            sAug = [cos(s); sin(s)];
            C = sAug * diag(w) * sAug';
            testCase.verifyEqual(C, B.moment(), 'RelTol', 1E-3);
        end
    end
    
    methods (Static)
        function B = getRandomBingham(dim)
            C = rand(dim,dim);
            [M,~] = qr(C); % obtain orthogonal matrix
            Z = [sort(-exp(5*rand(dim-1,1))); 0];
            B = BinghamDistribution(Z,M);
        end
    end
    
end

