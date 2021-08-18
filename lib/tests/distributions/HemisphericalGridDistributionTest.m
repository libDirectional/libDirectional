classdef HemisphericalGridDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testWarningAsymm(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vmf = VMFDistribution(1/sqrt(2)*[-1;0;-1],2);
            mixture = HypersphericalMixture(...
                {VMFDistribution(1/sqrt(2)*[-1;0;-1],2),VMFDistribution(1/sqrt(2)*[1;0;-1],2)},[0.5,0.5]);
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            testCase.verifyWarning(@()HemisphericalGridDistribution.fromDistribution(vmf, 1012), 'fromDistribution:asymmetricOnHypersphere');
            testCase.verifyWarning(@()HemisphericalGridDistribution.fromDistribution(mixture, 1012), 'fromDistribution:asymmetricOnHypersphere');
        end
        
        function testApproxVMMixtureS2(testCase)
            dist = HypersphericalMixture(...
                {VMFDistribution(1/sqrt(2)*[-1;0;1],2),VMFDistribution(1/sqrt(2)*[1;0;-1],2)},[0.5,0.5]);
            
            hhgd = HemisphericalGridDistribution.fromDistribution(dist, 1012);
            testCase.verifyEqual(hhgd.gridValues', 2*dist.pdf(hhgd.getGrid()));
        end
        
        function testApproxBinghamS2(testCase)
            M = eye(3);
            Z = [-2 -1 0]';
            dist = BinghamDistribution(Z,M);
            % Improve normalization constant for Bingham distribution
            dist.F = dist.F*dist.integralNumerical;
            
            hhgd = HemisphericalGridDistribution.fromDistribution(dist, 1012);
            testCase.verifyEqual(hhgd.gridValues', 2*dist.pdf(hhgd.getGrid()));
        end
        
        function testApproxBinghamS3(testCase)
            M = eye(4);
            Z = [-2 -1 -0.5 0]';
            dist = BinghamDistribution(Z,M);
            % Improve normalization constant for Bingham distribution
            dist.F = dist.F*dist.integralNumerical;
            
            hhgd = HemisphericalGridDistribution.fromDistribution(dist, 1012);
            testCase.verifyEqual(hhgd.gridValues', 2*dist.pdf(hhgd.getGrid()));
        end
        
        % Test operations for prediction and filter steps
        function testMultiplyVMFMixtureS2(testCase)
            % Test 3d case by validating against SphericalGridDistribution
            for kappa1 = 0.1:0.3:4
                for kappa2 = 0.1:0.3:4
                    dist1 = HypersphericalMixture({VMFDistribution(1/sqrt(2)*[-1;0;1], kappa1),...
                        VMFDistribution(-1/sqrt(2)*[-1;0;1], kappa1)},[0.5,0.5]);
                    dist2 = HypersphericalMixture({VMFDistribution([0;-1;0], kappa2),...
                        VMFDistribution([0;1;0], kappa2)},[0.5,0.5]);
                    hhgd1 = HemisphericalGridDistribution.fromDistribution(dist1, 1000, 'eq_point_set');
                    hhgd2 = HemisphericalGridDistribution.fromDistribution(dist2, 1000, 'eq_point_set');
                    hhgdFiltered = hhgd1.multiply(hhgd2);
                    
                    hgd1 = HypersphericalGridDistribution.fromDistribution(dist1, 2000, 'eq_point_set_symm');
                    hgd2 = HypersphericalGridDistribution.fromDistribution(dist2, 2000, 'eq_point_set_symm');
                    hgdFiltered = hgd1.multiply(hgd2);
                    
                    currentGrid = hgdFiltered.getGrid();
                    testCase.verifyEqual(hhgdFiltered.getGrid(), currentGrid(:,1:1000));
                    testCase.verifyEqual(0.5*hhgdFiltered.gridValues, hgdFiltered.gridValues(1:1000), 'AbsTol', 1e-11);
                end
            end
        end
        function testMultiplyError(testCase)
            dist1 = HypersphericalMixture({VMFDistribution(1/sqrt(2)*[-1;0;1], 1),...
                VMFDistribution(-1/sqrt(2)*[-1;0;1], 1)},[0.5,0.5]);
            f1 = HemisphericalGridDistribution.fromDistribution(dist1, 84, 'eq_point_set');
            f2 = f1;
            f2.gridValues = f2.gridValues(1:end-1);
            currGrid = f2.getGrid();
            f2.grid = currGrid(:,1:end-1);
            testCase.verifyError(@()f1.multiply(f2),'Multiply:IncompatibleGrid');
        end
        function testToFullSphere(testCase)
            dist = HypersphericalMixture({VMFDistribution(1/sqrt(2)*[-1;0;1], 1),...
                VMFDistribution(-1/sqrt(2)*[-1;0;1], 1)},[0.5,0.5]);
            hgd = HypersphericalGridDistribution.fromDistribution(dist, 84, 'eq_point_set_symm');
            hhgd = HemisphericalGridDistribution.fromDistribution(dist, 42, 'eq_point_set_symm');
            
            hhgd2hgd = hhgd.toFullSphere;
            testCase.verifyEqual(hhgd2hgd.getGrid(), hgd.getGrid());
            testCase.verifyEqual(hhgd2hgd.gridValues, hgd.gridValues);
        end
    end
end
