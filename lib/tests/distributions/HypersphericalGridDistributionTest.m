classdef HypersphericalGridDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testApproxVMMixtureS2(testCase)
            dist = HypersphericalMixture(...
                {VMFDistribution(1/sqrt(2)*[-1;0;1],2),VMFDistribution([0;-1;0],2)},[0.5,0.5]);
            
            hgd = HypersphericalGridDistribution.fromDistribution(dist, 1012);
            sgd = SphericalGridDistribution.fromDistribution(dist, 1012);
            testCase.verifyEqual(hgd.getGrid(),sgd.getGrid());
            testCase.verifyEqual(hgd.gridValues,sgd.gridValues);
        end
        
        function testApproxVMMixtureSd(testCase)
            rng(1) % Make it reproducible
            for i=2:6
                % Multiplication has same mean as mixture of two
                % equally-weighted ones
                vmf1 = VMFDistribution(HypersphericalUniformDistribution(i).sample(1),2);
                vmf2 = VMFDistribution(HypersphericalUniformDistribution(i).sample(1),2);
                vmfMult = vmf1.multiply(vmf2);
                dist = HypersphericalMixture({vmf1,vmf2},[0.5,0.5]);
                hgd = HypersphericalGridDistribution.fromDistribution(dist, 1012);
                testCase.verifyEqual(hgd.meanDirection,vmfMult.meanDirection,'AbsTol',5e-2);
            end
        end
        
        function testApproxBinghamS2(testCase)
            M = eye(3);
            Z = [-2 -1 0]';
            dist = BinghamDistribution(Z,M);
            % Improve normalization constant for Bingham distribution
            dist.F=dist.F*dist.integralNumerical;
            
            
            hgd = HypersphericalGridDistribution.fromDistribution(dist, 1012);
            sgd = SphericalGridDistribution.fromDistribution(dist, 1012);
            
            % First verify that SphericalGridDistribution approxiamtes it
            % correctly, then verify HypersphericalGridDistribution has
            % equal grid and grid values
            verifyPdfEqual(testCase,sgd,dist,1e-6);
            testCase.verifyEqual(hgd.getGrid(),sgd.getGrid());
            testCase.verifyEqual(hgd.gridValues,sgd.gridValues);
        end
        
        function testMeanDirectionSd(testCase)
            for i=2:5
                rng(1) % Make it reproducible
                mu = HypersphericalUniformDistribution(i).sample(1);
                vmf = VMFDistribution(mu, 2);
                hgd = HypersphericalGridDistribution.fromDistribution(vmf, 1012);
                testCase.verifyEqual(hgd.meanDirection, mu, 'AbsTol', 2e-2);
            end
        end
        % Test operations for prediction and filter steps
        function testMultiplyVMFS2(testCase)
            % Test 3d case by validating against SphericalGridDistribution
            for kappa1 = 0.1:0.3:4
                for kappa2 = 0.1:0.3:4
                    dist1 = VMFDistribution(1/sqrt(2)*[-1;0;1], kappa1);
                    dist2 = VMFDistribution([0;-1;0], kappa2);
                    hgd1 = HypersphericalGridDistribution.fromDistribution(dist1, 1000, 'eq_point_set');
                    hgd2 = HypersphericalGridDistribution.fromDistribution(dist2, 1000, 'eq_point_set');
                    hgdFiltered = hgd1.multiply(hgd2);
                    
                    sgd1 = SphericalGridDistribution.fromDistribution(dist1, 1000, 'eq_point_set');
                    sgd2 = SphericalGridDistribution.fromDistribution(dist2, 1000, 'eq_point_set');
                    sgdFiltered = sgd1.multiply(sgd2);
                    
                    testCase.verifyEqual(hgdFiltered.gridValues, sgdFiltered.gridValues, 'AbsTol', 1e-6);
                end
            end
        end
        function testMultiplyVMFSd(testCase)
            % Test 3d case by validating against SphericalGridDistribution
            kappa1 = 2;
            kappa2 = 1;
            for i=2:6
                rng(1) % Make it reproducible
                % Multiplication has same mean as mixture of two
                % equally-weighted ones
                vmf1 = VMFDistribution(HypersphericalUniformDistribution(i).sample(1),kappa1);
                vmf2 = VMFDistribution(HypersphericalUniformDistribution(i).sample(1),kappa2);
                vmfMult = vmf1.multiply(vmf2);
                
                hgd1 = HypersphericalGridDistribution.fromDistribution(vmf1, 1000, 'eq_point_set');
                hgd2 = HypersphericalGridDistribution.fromDistribution(vmf2, 1000, 'eq_point_set');
                hgdMult = hgd1.multiply(hgd2);
                testCase.verifyEqual(hgdMult.meanDirection,vmfMult.meanDirection,'AbsTol',5e-2);
            end
        end
        function testMultiplyError(testCase)
            dist1 = VMFDistribution(1/sqrt(2)*[-1;0;1], 1);
            f1 = HypersphericalGridDistribution.fromDistribution(dist1, 84, 'eq_point_set');
            f2 = f1;
            f2.gridValues = f2.gridValues(1:end-1);
            gridFull = f2.getGrid();
            f2.grid = gridFull(:,1:end-1);
            testCase.verifyError(@()f1.multiply(f2),'Multiply:IncompatibleGrid');
        end
        function testSymmetrizeVMFMixtureS2(testCase)
            dist = HypersphericalMixture(...
                {VMFDistribution([0;1;0],2),VMFDistribution([0;-1;0],2)},[0.5,0.5]);
            f = HypersphericalGridDistribution.fromDistribution(dist, 50, 'eq_point_set_symm');
            testCase.verifyEqual(f.gridValues(size(f.getGrid(),2)/2+1:end),f.gridValues(1:size(f.getGrid(),2)/2));
            fAsymm = f;
            fAsymm.gridValues(26:27) = fAsymm.gridValues([27,26]);
            testCase.verifyNotEqual(fAsymm.gridValues(size(fAsymm.getGrid(),2)/2+1:end),fAsymm.gridValues(1:size(fAsymm.getGrid(),2)/2));
            fSymm = fAsymm.symmetrize;
            testCase.verifyEqual(fSymm.gridValues(size(fSymm.getGrid(),2)/2+1:end),fSymm.gridValues(1:size(fSymm.getGrid(),2)/2));
            testCase.verifyNotEqual(fSymm.gridValues,fAsymm.gridValues);
            % Because the average is and distributed to both points, the new
            % density is not equal to original one.
            testCase.verifyNotEqual(fSymm.gridValues,f.gridValues);
        end
        function testSymmetrizeWatsonS3(testCase)
            dist = WatsonDistribution(1/sqrt(2)*[1;1;0],1);
            f = HypersphericalGridDistribution.fromDistribution(dist, 50, 'eq_point_set_symm');
            testCase.verifyEqual(f.gridValues(size(f.getGrid(),2)/2+1:end),f.gridValues(1:size(f.getGrid(),2)/2));
            fAsymm = f;
            fAsymm.gridValues([26,46]) = fAsymm.gridValues([46,26]);
            testCase.verifyNotEqual(fAsymm.gridValues(size(fAsymm.getGrid(),2)/2+1:end),fAsymm.gridValues(1:size(fAsymm.getGrid(),2)/2));
            fSymm = fAsymm.symmetrize;
            testCase.verifyEqual(fSymm.gridValues(size(fSymm.getGrid(),2)/2+1:end),fSymm.gridValues(1:size(fSymm.getGrid(),2)/2));
            testCase.verifyNotEqual(fSymm.gridValues,fAsymm.gridValues);
            % Because the average is and distributed to both points, the new
            % density is not equal to original one.
            testCase.verifyNotEqual(fSymm.gridValues,f.gridValues);
        end
        function testSymmetrizeError(testCase)
            dist = VMFDistribution(1/sqrt(2)*[-1;0;1], 1);
            f = HypersphericalGridDistribution.fromDistribution(dist, 84, 'eq_point_set');
            testCase.verifyError(@()f.symmetrize,'Symmetrize:AsymmetricGrid');
        end
    end
    methods
        function verifyPdfEqual(testCase,dist1,dist2,tol)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            fixture1 = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            fixture2 = testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(dist1.pdf([x; y; z]), dist2.pdf([x; y; z]), 'AbsTol', tol);
            fixture1.teardown;
            fixture2.teardown;
        end
    end
end
