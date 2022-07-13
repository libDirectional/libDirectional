classdef FIGFilterTest < matlab.unittest.TestCase
    
    methods(Test)
        function testInitialization(testCase)
            filter = FIGFilter(101);
            testCase.verifyClass(filter.gd, 'FIGDistribution');
            testCase.verifyLength(filter.gd.gridValues, 101);
            testCase.verifyEqual(filter.gd.enforcePdfNonnegative, true);
        end
        
        function testSetState(testCase)
            filter = FIGFilter(101);
            gd = FIGDistribution.fromDistribution(VMDistribution(2, 3), 101);
            filter.setState(gd);
            testCase.verifyEqual(filter.gd, gd);
            testCase.verifyWarning(@()filter.setState(VMDistribution(2, 3)), 'setState:nonGrid');
            gdtmp = gd;
            gdtmp.gridValues = 1:3;
            testCase.verifyWarning(@()filter.setState(gdtmp), 'setState:noOfGridValuesDiffer');
            testCase.verifyWarning(@()filter.setState(VMDistribution(2, 3)), 'setState:nonGrid');
        end
        
        function testPrediction(testCase)
            % Validity of operation is tested in Fourier class, only test
            % filter interface
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Predict:automaticConversion'));
            for currEnforcement = [false, true]
                filter = FIGFilter(101, currEnforcement);
                vm = VMDistribution(3, 2);
                filter.predictIdentity(vm);
                fd1 = FIGDistribution.fromDistribution(CircularUniformDistribution(), 101, currEnforcement);
                fd2 = fd1.convolve(FIGDistribution.fromDistribution(vm, 101, currEnforcement));
                testCase.verifyLength(filter.gd.gridValues, 101);
                testCase.verifyEqual(filter.gd.gridValues, fd2.gridValues, 'AbsTol', 1E-8);
            end
        end
        
        function testUpdateIdentity(testCase)
            % As above, only test filter interface
            for currEnforcement = [false, true]
                filter = FIGFilter(101, currEnforcement);
                vmMultiply = VMDistribution(3, 2);
                vmFilter = VMDistribution(0, 2);
                filter.updateIdentity(vmFilter, 3);
                gd1 = FIGDistribution.fromDistribution(CircularUniformDistribution(), 101, currEnforcement);
                gd2 = gd1.multiply(FIGDistribution.fromDistribution(vmMultiply, 101, currEnforcement));
                testCase.verifyLength(filter.gd.gridValues, 101);
                testCase.verifyEqual(filter.gd.gridValues, gd2.gridValues, 'AbsTol', 1E-8);
            end
        end
        
        function testUpdateNonlinear(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            z = 3;
            likelihood = @(z, x)1 ./ (abs(x-z) + .1);
            for currEnforcement = [false, true]
                filter = FIGFilter(23, currEnforcement);
                hfd1 = FIGDistribution.fromDistribution(VMDistribution(3, 2), 23, currEnforcement);
                filter.setState(hfd1);
                filter.updateNonlinear(likelihood, z);
                % The mean should stay approximately equal, but not entirely
                testCase.verifyEqual(filter.getEstimateMean, z, 'AbsTol', 0.2);
                % Verify that it is normalized by validating no
                % normalization is perform when converting to Fourier-based
                % representation
                testCase.verifyWarningFree(@()filter.getEstimate.integral);
            end
        end
    end
    
end
