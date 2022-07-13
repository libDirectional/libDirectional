classdef StateSpaceSubdivisionGaussianDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testMultiplyS1xR1identicalprecise(testCase)
            % Other eval also in HypercylincdricalStateSpaceSubdivisionGaussian
            n = 100;
            rbd1 = StateSpaceSubdivisionGaussianDistribution(FIGDistribution.fromDistribution(...
                CircularUniformDistribution,n),...
                repmat(GaussianDistribution(0, 1),[n,1]));
            
            rbd2 = StateSpaceSubdivisionGaussianDistribution(FIGDistribution.fromDistribution(...
                CircularUniformDistribution,n),...
                repmat(GaussianDistribution(2, 1),[n,1]));
            
            rbdUp = rbd1.multiply(rbd2);
            % Verify uncertainty in linear domain got lower
            for i=1:100
                testCase.verifyLessThan(det(rbdUp.linearDistributions(1).C),det(rbd1.linearDistributions(1).C));
                testCase.verifyEqual(rbdUp.linearDistributions(i).mean, 1);
            end
        end
        
        function testMultiplyS2xR3rough(testCase)
            n = 100;
            rbd1 = StateSpaceSubdivisionGaussianDistribution(HemisphericalGridDistribution.fromDistribution(...
                HyperhemisphericalUniformDistribution(3),n),...
                repmat(GaussianDistribution([0;0;0], 1000*eye(3)),[n,1]));
            rbd2 = StateSpaceSubdivisionGaussianDistribution(HemisphericalGridDistribution.fromDistribution(...
                HyperhemisphericalUniformDistribution(3),n),...
                repmat(GaussianDistribution([2;2;2], 1000*eye(3)),[n,1]));
            
            rbdUp = rbd1.multiply(rbd2);
            % Verify uncertainty in linear domain got lower
            for i=1:100
                testCase.verifyLessThan(det(rbdUp.linearDistributions(1).C),det(rbd1.linearDistributions(1).C));
                testCase.verifyEqual(rbdUp.linearDistributions(i).mean, [1;1;1]);
            end
        end
        
        function testHybridMean(testCase)
            n = 100;
            muPeriodic = 4;
            muLinear = [1;2;3];
            rbd = StateSpaceSubdivisionGaussianDistribution(FIGDistribution.fromDistribution(...
                VMDistribution(muPeriodic,1),n),...
                repmat(GaussianDistribution(muLinear, 1000*eye(3)),[n,1]));
            testCase.verifyEqual(rbd.hybridMean(),[muPeriodic;muLinear],'AbsTol',1e-15);
        end
        
        function testLinearMean(testCase)
            n = 100;
            muPeriodic = 4;
            muLinear = [1;2;3];
            rbd = StateSpaceSubdivisionGaussianDistribution(FIGDistribution.fromDistribution(...
                VMDistribution(muPeriodic,1),n),...
                repmat(GaussianDistribution(muLinear, 1000*eye(3)),[n,1]));
            testCase.verifyEqual(rbd.linearMean(),muLinear,'AbsTol',1e-15);
        end
        
        function testModeWarningUniform(testCase)
            n = 100;
            muLinear = [1;2;3];
            rbd = StateSpaceSubdivisionGaussianDistribution(FIGDistribution.fromDistribution(...
                CircularUniformDistribution(),n),...
                repmat(GaussianDistribution(muLinear, 1000*eye(3)),[n,1]));
            testCase.verifyWarning(@()rbd.mode(),'Mode:PotentiallyMultimodal');
        end
        
        function testMode(testCase)
            n = 100;
            muPeriodic = 4;
            muLinear = [1;2;3];
            rbd = StateSpaceSubdivisionGaussianDistribution(FIGDistribution.fromDistribution(...
                VMDistribution(muPeriodic,10),n),...
                repmat(GaussianDistribution(muLinear, eye(3)),[n,1]));
            
            testCase.verifyWarningFree(@()rbd.mode(),'Mode:PotentiallyMultimodal');
            % Since we land somewhere on the grid, it is expected we that
            % we are not indefinitely precise.
            testCase.verifyEqual(rbd.mode(),[muPeriodic;muLinear],'AbsTol',pi/n);
        end
    end
        
end