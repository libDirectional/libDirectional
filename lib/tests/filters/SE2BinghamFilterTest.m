classdef SE2BinghamFilterTest < matlab.unittest.TestCase
   
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
        
        function testUpdateNonlinear(testCase)
            function p = likelihoodBf(z,x)
                gaussianMeasNoiseCov = 0.5*eye(2);
                [~,pos] = AbstractSE2Distribution.dualQuaternionToAnglePos(x);
                p = mvnpdf(pos',z',gaussianMeasNoiseCov)';
            end
            sampleForFitting = 10000;
            initialPriorPeriodic = VMDistribution(0,1);
            initialPriorLinear = GaussianDistribution([0;0],eye(2));
            initialPrior = SE2CartProdStackedDistribution({initialPriorPeriodic;initialPriorLinear});
            
            se2bf = SE2BinghamFilter();
            se2bf.setState(SE2BinghamDistribution.fit(initialPrior.sample(sampleForFitting)));
            
            stateAndMeasForced = [1;1];
            for i=1:100
                se2bf.updateProgressive(@likelihoodBf,stateAndMeasForced);
            end
            est = se2bf.getPointEstimate(true);
            testCase.verifyEqual(est(2:3),stateAndMeasForced,'AbsTol',1e-4);
        end
        
        function testPredictNonlinear(testCase)
            rng default
            muPriorLinear = [5;2];
            initialPriorLinear = GaussianDistribution(muPriorLinear, 0.1*eye(2));
            
            vmSysNoise = VMDistribution(0, 10);
            gaussianSysNoise = GaussianDistribution([0; 0], 0.1*eye(2));
            sysNoise = SE2CartProdStackedDistribution( ...
                {vmSysNoise, gaussianSysNoise});

            stepSize = 2.5;
            genNextStateWithoutNoise = ...
                        @(x)[x(1, :); x(2:3, :) + stepSize * [cos(x(1)); sin(x(1))]];

            noSamplesForFitting = 10000;
            sysNoiseForFilter = SE2BinghamDistribution.fit(sysNoise.sample(noSamplesForFitting));
            
            se2bf = SE2BinghamFilter();
            for currPriorAngle = linspace(-pi/2,5/2*pi,10)
                initialPriorPeriodic = VMDistribution(currPriorAngle, 1);
                initialPrior = SE2CartProdStackedDistribution( ...
                    {initialPriorPeriodic; initialPriorLinear});
                
                priorForFilter = SE2BinghamDistribution.fit(initialPrior.sample(noSamplesForFitting));               
                se2bf.setState(priorForFilter);
                testCase.verifyEqual(se2bf.getPointEstimate,[mod(currPriorAngle,2*pi); muPriorLinear], 'AbsTol', 0.1);
                
                se2bf.predictNonlinear(genNextStateWithoutNoise, sysNoiseForFilter);
                testCase.verifyEqual(se2bf.getPointEstimate,...
                    genNextStateWithoutNoise([mod(currPriorAngle,2*pi); muPriorLinear]), 'AbsTol', 0.1);
            end
        end
    end
end
