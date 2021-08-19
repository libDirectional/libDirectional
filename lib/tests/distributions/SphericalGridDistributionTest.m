classdef SphericalGridDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testApproxVMMixture(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            dist = HypersphericalMixture(...
                {VMFDistribution(1/sqrt(2)*[-1;0;1],2),VMFDistribution([0;-1;0],2)},[0.5,0.5]);
            
            sgd = SphericalGridDistribution.fromDistribution(dist, 1012, 'eq_point_set');
            verifyPdfEqual(testCase,sgd,dist,1e-11);
            sgd = SphericalGridDistribution.fromDistribution(dist, 1012, 'eq_point_set', true);
            verifyPdfEqual(testCase,sgd,dist,1e-11);
            sgd = SphericalGridDistribution.fromDistribution(dist, 1012, 'eq_point_set', false);
            verifyPdfEqual(testCase,sgd,dist,1e-11);
            
            % For spherical harmonics grid, normalization function does not
            % properly normalize (just keeps it in the reasonable range).
            testCase.applyFixture(SuppressedWarningsFixture('Transformation:notEq_Point_set'));
            testCase.applyFixture(SuppressedWarningsFixture('SphericalGridDistribution:CannotNormalizeShGrid'));
            sgd = SphericalGridDistribution.fromDistribution(dist, 1012, 'sh_grid');
            verifyPdfEqual(testCase,sgd,dist,1e-11);
            sgd = SphericalGridDistribution.fromDistribution(dist, 1012, 'sh_grid', true);
            verifyPdfEqual(testCase,sgd,dist,1e-11);
            sgd = SphericalGridDistribution.fromDistribution(dist, 1012, 'sh_grid', false);
            verifyPdfEqual(testCase,sgd,dist,1e-11);
        end
        
        function testApproxBingham(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            M = eye(3);
            Z = [-2 -1 0]';
            dist = BinghamDistribution(Z,M);
            % Improve normalization constant for Bingham distribution
            dist.F=dist.F*dist.integralNumerical;
            
            testCase.applyFixture(SuppressedWarningsFixture('Transformation:notEq_Point_set'));
            for grid_pattern = {'eq_point_set', 'sh_grid'}
                if strcmp([grid_pattern{:}],'sh_grid')
                    testCase.applyFixture(SuppressedWarningsFixture('Transformation:notEq_Point_set'));
                    testCase.applyFixture(SuppressedWarningsFixture('SphericalGridDistribution:CannotNormalizeShGrid'));
                end
                
                sgd = SphericalGridDistribution.fromDistribution(dist, 1012, [grid_pattern{:}]);
                verifyPdfEqual(testCase,sgd,dist,1e-6);

                sgd = SphericalGridDistribution.fromDistribution(dist, 1012, [grid_pattern{:}], true);
                verifyPdfEqual(testCase,sgd,dist,1e-6);

                sgd = SphericalGridDistribution.fromDistribution(dist, 1012, [grid_pattern{:}], false);
                verifyPdfEqual(testCase,sgd,dist,1e-6);
            end
        end
        
        function testErrorHigherDim(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            dist = VMFDistribution(1/sqrt(2)*[-1;0;1;0],2);
            testCase.verifyError(@()SphericalGridDistribution.fromDistribution(dist, 1012, 'eq_point_set'),...
                'fromDistribution:WrongDimension');
        end
        
        function testNormalization(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            dist = HypersphericalMixture(...
                {VMFDistribution(1/sqrt(2)*[-1;0;1],2),VMFDistribution([0;-1;0],2)},[0.5,0.5]);
            sgd = SphericalGridDistribution.fromDistribution(dist, 1012, 'eq_point_set');
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(sgd.integral,1,'AbsTol',1e-2);
        end
        function testPlottingWarningFree(testCase)
            dist = HypersphericalMixture(...
                {VMFDistribution(1/sqrt(2)*[-1;0;1],2),VMFDistribution([0;-1;0],2)},[0.5,0.5]);
            sgd = SphericalGridDistribution.fromDistribution(dist, 1012, 'eq_point_set');
            testCase.verifyWarningFree(@()sgd.plot);
            testCase.verifyWarningFree(@()sgd.plotInterpolated);
        end
        function testMeanDirection(testCase)
            mu = 1/sqrt(2)*[-1;0;1];
            dist = VMFDistribution(mu,2);
            sgd = SphericalGridDistribution.fromDistribution(dist, 1000000);
            testCase.verifyEqual(sgd.meanDirection, mu, 'AbsTol', 1e-6);
        end
        % Test operations for prediction and filter steps
        function testMultiplyVMF(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            for enforceNonnegative=[false,true]
                for kappa1 = 0.1:0.3:4
                    for kappa2 = 0.1:0.3:4
                        dist1 = VMFDistribution(1/sqrt(2)*[-1;0;1], kappa1);
                        dist2 = VMFDistribution([0;-1;0], kappa2);
                        f1 = SphericalGridDistribution.fromDistribution(dist1, 1000, 'eq_point_set', enforceNonnegative);
                        f2 = SphericalGridDistribution.fromDistribution(dist2, 1000, 'eq_point_set', enforceNonnegative);
                        fFiltered = f1.multiply(f2);
                        distResult = dist1.multiply(dist2);
                        % The tolerance for the spherical harmonics may be
                        % different here, prevent warnings
                        fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
                        verifyPdfEqual(testCase,fFiltered,distResult,1e-5);
                        fixture.teardown();
                    end
                end
            end
        end
        function testMultiplyError(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            dist1 = VMFDistribution(1/sqrt(2)*[-1;0;1], 1);
            dist2 = VMFDistribution([0;-1;0], 2);
            f1 = SphericalGridDistribution.fromDistribution(dist1, 84, 'eq_point_set', true);
            testCase.applyFixture(SuppressedWarningsFixture('Transformation:notEq_Point_set'));
            f2 = SphericalGridDistribution.fromDistribution(dist2, 84, 'sh_grid', false);
            testCase.verifyError(@()f1.multiply(f2),'Multiply:IncompatibleGrid');
        end
    end
    methods
        function verifyPdfEqual(testCase,dist1,dist2,tol)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            fixture = testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            if strcmp(dist1.gridType,'sh_grid')
                fixture(2) = testCase.applyFixture(SuppressedWarningsFixture('SphericalGridDistribution:CannotNormalizeShGrid'));
                fixture(3) = testCase.applyFixture(SuppressedWarningsFixture('PDF:NeedNormalizationShGrid'));
            end
            testCase.verifyEqual(dist1.pdf([x; y; z]), dist2.pdf([x; y; z]), 'AbsTol', tol);
            for f=fixture
                f.teardown;
            end
        end
    end
end
