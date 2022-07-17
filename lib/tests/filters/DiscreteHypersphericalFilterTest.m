classdef DiscreteHypersphericalFilterTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testDiscreteFilter(testCase)
            for d=2:4
                nParticles = 300;
                filter = DiscreteHypersphericalFilter(nParticles, d);
                vmf = VMFDistribution([1;zeros(d-1,1)], 2);
                vmfFilter = VMFFilter();

                % test VMF approximation
                filter.setState(vmf);
                vmfFilter.setState(vmf);
                %testCase.verifyEqual(filter.getEstimateMean(), vmfFilter.getEstimateMean(), 'AbsTol', 2E-2);
                testCase.verifyLessThan(acos(filter.getEstimateMean()' * vmfFilter.getEstimateMean()), 0.01);

                % test update
                z = ones(d,1)/sqrt(d);
                vmf2 = VMFDistribution([zeros(d-1,1); 1], 0.8);
                vmfFilter.updateIdentity(vmf2, z);
                filter.updateNonlinear(@(z,x) VMFDistribution(z, vmf2.kappa).pdf(x), z);
                %testCase.verifyEqual(filter.getEstimateMean(), vmfFilter.getEstimateMean(), 'AbsTol', 2E-2);
                testCase.verifyLessThan(acos(filter.getEstimateMean()' * vmfFilter.getEstimateMean()), 0.01);

                % test prediction with arbitrary noise
                noiseSamples = vmf2.sampleDeterministic();
                noiseWeights = ones(1, size(noiseSamples,2));
                f = @(x, w) (x + w)/norm(x+w);
                filter.predictNonlinearNonAdditive(f, noiseSamples, noiseWeights);
                if d==3
                    filter.predictNonlinearNonAdditive(f, noiseSamples, noiseWeights, 'delaunay');
                end
                filter.predictNonlinearNonAdditive(f, noiseSamples, noiseWeights, 'knn');

                % test prediction with vmf noise
                filter.setState(vmf);
                vmfFilter.setState(vmf);
                filter.predictNonlinear(@(x) x, vmf2);
                vmfFilter.predictIdentity(vmf2);
                testCase.verifyLessThan(acos(filter.getEstimateMean()' * vmfFilter.getEstimateMean()), 0.01);
            end
        end
        
        function testDiscretizationMethods(testCase)
            nParticles = 50;
            
            for dim = 2:5
                switch dim
                    case 2
                        methods = {'equalarea', 'random', 'optimize'};
                    case 3
                        methods = {'equalarea', 'random', 'icosahedronsubdivision', 'cubesubdivision', 'gridsphere', 'spiral', 'stratified', 'reisz'}; 
                    case 4
                        methods = {'equalarea', 'random', 'schaefer', 'tesseractsubdivision'};
                    case 5
                        methods = {'equalarea', 'random'};
                    otherwise
                        error('should never happen')
                end
                
                for i = 1:length(methods)
                    filter = DiscreteHypersphericalFilter(nParticles, dim, methods{i});
                    testCase.verifyEqual(sum(filter.d.^2), ones(1,size(filter.d,2)), 'RelTol', 1E-10);
                    testCase.verifyGreaterThanOrEqual(size(filter.d,2), nParticles);
                end
            end
            
            %test platonic solids
            for n = [4 6 8 12 20]
                filter = DiscreteHypersphericalFilter(n, 3, 'platonic');
                testCase.verifyEqual(sum(filter.d.^2), ones(1,size(filter.d,2)), 'RelTol', 1E-10);
                testCase.verifyEqual(mean(filter.d,2), zeros(3,1), 'AbsTol', 1E-10);
                testCase.verifyEqual(size(filter.d,2), n);
            end
            
            %test regular 4 polytopes
            for n = [24 120 600]
                filter = DiscreteHypersphericalFilter(n, 4, 'regular4polytope');
                testCase.verifyEqual(sum(filter.d.^2), ones(1,size(filter.d,2)), 'RelTol', 1E-10);
                testCase.verifyEqual(mean(filter.d,2), zeros(4,1), 'AbsTol', 1E-10);
                testCase.verifyEqual(size(filter.d,2), n);
            end
        end
    end
end
