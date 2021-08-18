classdef SphericalHarmonicsFilterTest < matlab.unittest.TestCase
    
    methods(Test)
        function testUpdateIdentity(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            for transformation = {'identity', 'sqrt'}
                shdFilter = SphericalHarmonicsFilter(30, [transformation{:}]);
                vmfFilter = VMFFilter();
                
                vmf1 = VMFDistribution([0; 1; 0], 1);
                vmf2 = VMFDistribution([0; 0; 1], 0.1);
                shd1 = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(vmf1, 30);
                shd2 = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(vmf2, 30);
                
                fixture = testCase.applyFixture(SuppressedWarningsFixture('SphericalHarmonicsFilter:rotationRequired'));
                vmfFilter.setState(vmf1);
                vmfFilter.updateIdentity(vmf2, [1; 0; 0]);
                
                shdFilter.setState(shd1);
                shdFilter.updateIdentity(shd2, [1; 0; 0]);
                fixture.teardown;
                
                testCase.verifyEqual(vmfFilter.getEstimateMean, shdFilter.getEstimateMean, 'AbsTol', 1E-10);
            end
        end
        function testUpdateUsingLikelihood(testCase)
            rng(1)
            posTrue = -1 / sqrt(3) * ones(3, 1);
            
            %% generate measurements according to truncated gaussian along x and y axis
            measX = NaN(1, 5);
            sigmaX = .3;
            measY = NaN(1, 5);
            sigmaY = .3;
            measZ = NaN(1, 5);
            sigmaZ = .3;
            
            for i = 1:numel(measX)
                measX(i) = normrnd(posTrue(1), sigmaX);
                if measX(i) > 1 || measX(i) < -1
                    i = i - 1; %#ok<FXSET>
                end
            end
            for i = 1:numel(measY)
                measY(i) = normrnd(posTrue(2), sigmaY);
                if measY(i) > 1 || measY(i) < -1
                    i = i - 1; %#ok<FXSET>
                end
            end
            for i = 1:numel(measZ)
                measZ(i) = normrnd(posTrue(3), sigmaZ);
                if measZ(i) > 1 || measZ(i) < -1
                    i = i - 1; %#ok<FXSET>
                end
            end
            for transformation = {'identity', 'sqrt'}
                shFilter = SphericalHarmonicsFilter(11, [transformation{:}]); % Initialize with uniform
                
                %% Update first using measurements of the x-axis, then of the y-axis and then of the z-axis
                for i = 1:numel(measX)
                    shFilter.updateNonlinear(@(z, x)normpdf(z(1), x(1, :), sigmaX), [measX(i); 0; 0])
                end
                for i = 1:numel(measY)
                    shFilter.updateNonlinear(@(z, x)normpdf(z(2), x(2, :), sigmaY), [0; measY(i); 0])
                end
                for i = 1:numel(measZ)
                    shFilter.updateNonlinear(@(z, x)normpdf(z(3), x(3, :), sigmaZ), [0; 0; measZ(i)])
                end
                testCase.verifyEqual(shFilter.getEstimateMean, posTrue, 'AbsTol', .3);
            end
        end
        function testUpdateUsingLikelihoodMultiple(testCase)
            % Test if approximately equal to update sequentially or at
            % once.
            for transformation = {'identity', 'sqrt'}
                shFilter1 = SphericalHarmonicsFilter(10, [transformation{:}]);
                shFilter2 = SphericalHarmonicsFilter(10, [transformation{:}]);
                
                sigmaX = 0.3;
                sigmaY = 0.3;
                sigmaZ = 0.3;
                shFilter1.updateNonlinear(@(z, x)normpdf(z(1), x(1, :), sigmaX), [-1 / sqrt(3); 0; 0]);
                shFilter1.updateNonlinear(@(z, x)normpdf(z(2), x(2, :), sigmaX), [0; -1 / sqrt(3); 0]);
                shFilter1.updateNonlinear(@(z, x)normpdf(z(3), x(3, :), sigmaX), [0; 0; -1 / sqrt(3)]);
                
                shFilter2.updateNonlinear({@(z, x)normpdf(z(1), x(1, :), sigmaX), ...
                    @(z, x)normpdf(z(2), x(2, :), sigmaY), @(z, x)normpdf(z(3), x(3, :), sigmaZ)}, ...
                    {[-1 / sqrt(3); 0; 0], [0; -1 / sqrt(3); 0], [0; 0; -1 / sqrt(3)]});
                testCase.verifyEqual(shFilter2.getEstimateMean, shFilter1.getEstimateMean, 'AbsTol', 1E-5);
            end
        end
        function testPredictionSqrtVsID(testCase)
            degree = 21;
            densityInit = VMFDistribution([1; 1; 0]/sqrt(2), 2);
            sysNoise = VMFDistribution([0; 0; 1], 1);
            
            shdInitId = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(densityInit, degree, 'identity');
            shdInitSqrt = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(densityInit, degree, 'sqrt');
            shdNoiseId = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(sysNoise, degree, 'identity');
            shdNoiseSqrt = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(sysNoise, degree, 'sqrt');
            
            shFilter1 = SphericalHarmonicsFilter(degree, 'identity');
            shFilter2 = SphericalHarmonicsFilter(degree, 'sqrt');
            
            shFilter1.setState(shdInitId);
            shFilter2.setState(shdInitSqrt);
            
            shFilter1.predictIdentity(shdNoiseId);
            shFilter2.predictIdentity(shdNoiseSqrt);
            
            testCase.verifyEqual(shFilter1.state.totalVariationDistanceNumerical(shFilter2.state), 0, 'AbsTol', 1.5E-15);
        end
    end
end
