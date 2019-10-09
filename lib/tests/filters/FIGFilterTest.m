classdef FIGFilterTest < matlab.unittest.TestCase
    
    methods(Test)
        function testInitialization(testCase)
            filter = FIGFilter(101);
            testCase.verifyClass(filter.gd, 'FIGDistribution');
            testCase.verifyLength(filter.gd.gridValues, 101);
            testCase.verifyEqual(filter.gd.transformation, 'identity');
        end
        
        function testSetState(testCase)
            filter = FIGFilter(101);
            gd = FIGDistribution.fromDistribution(VMDistribution(2, 3), 101, 'identity');
            filter.setState(gd);
            testCase.verifyEqual(filter.gd, gd);
            testCase.verifyWarning(@()filter.setState(VMDistribution(2, 3)), 'setState:nonGrid');
            gdtmp = gd;
            gdtmp.gridValues = 1:3;
            testCase.verifyWarning(@()filter.setState(gdtmp), 'setState:noOfGridValuesDiffer');
            gdtmp = gd;
            gdtmp.transformation = 'sqrt';
            testCase.verifyWarning(@()filter.setState(gdtmp), 'setState:transDiffer');
            testCase.verifyWarning(@()filter.setState(VMDistribution(2, 3)), 'setState:nonGrid');
        end
        
        function testPrediction(testCase)
            % Validity of operation is tested in Fourier class, only test
            % filter interface
            for transformation = {'identity', 'sqrt'}
                filter = FIGFilter(101, [transformation{:}]);
                vm = VMDistribution(3, 2);
                warnstruct = warning('off', 'Predict:automaticConversion');
                filter.predictIdentity(vm);
                warning(warnstruct);
                fd1 = FIGDistribution.fromDistribution(CircularUniformDistribution(), 101, [transformation{:}]);
                fd2 = fd1.convolve(FIGDistribution.fromDistribution(vm, 101, [transformation{:}]));
                testCase.verifyLength(filter.gd.gridValues, 101);
                testCase.verifyEqual(filter.gd.gridValues, fd2.gridValues, 'AbsTol', 1E-8);
            end
        end
        
        function testUpdateIdentity(testCase)
            % As above, only test filter interface
            for transformation = {'identity', 'sqrt'}
                filter = FIGFilter(101, [transformation{:}]);
                vmMultiply = VMDistribution(3, 2);
                vmFilter = VMDistribution(0, 2);
                filter.updateIdentity(vmFilter, 3);
                gd1 = FIGDistribution.fromDistribution(CircularUniformDistribution(), 101, [transformation{:}]);
                gd2 = gd1.multiply(FIGDistribution.fromDistribution(vmMultiply, 101, [transformation{:}]));
                testCase.verifyLength(filter.gd.gridValues, 101);
                testCase.verifyEqual(filter.gd.gridValues, gd2.gridValues, 'AbsTol', 1E-8);
            end
        end
        
        function testUpdateNonlinear(testCase)
            z = 3;
            likelihood = @(z, x)1 ./ (abs(x-z) + .1);
            for transformation = {'identity', 'sqrt'}
                filter = FIGFilter(23, [transformation{:}]);
                hfd1 = FIGDistribution.fromDistribution(VMDistribution(3, 2), 23, [transformation{:}]);
                filter.setState(hfd1);
                warning('off', 'Normalization:notNormalized')
                filter.updateNonlinear(likelihood, z);
                warning('on', 'Normalization:notNormalized')
                testCase.verifyEqual(filter.getEstimateMean, z, 'AbsTol', 0.2); % The mean should stay approximately equal, but not entirely
            end
        end
    end
    
end
