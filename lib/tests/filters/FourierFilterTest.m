classdef FourierFilterTest < matlab.unittest.TestCase
   
    methods (Test)           
        function testInitialization(testCase)
            filter=FourierFilter(101);
            testCase.verifyClass(filter.fd, 'FourierDistribution');
            testCase.verifyLength(filter.fd.a,51);
            testCase.verifyLength(filter.fd.b,50);
            testCase.verifyEqual(filter.fd.transformation,'sqrt');
            testCase.verifyEqual(filter.fd.a(1),sqrt(2/pi),'AbsTol',1E-10);
        end
        
        function testSetState(testCase)
            filter=FourierFilter(101);
            fd=FourierDistribution.fromDistribution(VMDistribution(2,3),101);
            filter.setState(fd);
            testCase.verifyEqual(filter.fd,fd);
            testCase.verifyWarning(@()filter.setState(VMDistribution(2,3)),'setState:nonFourier');
            fdtmp=fd;
            fdtmp.a=fd.a(1:end-1);fdtmp.b=fd.b(1:end-1);
            testCase.verifyWarning(@()filter.setState(fdtmp),'setState:noOfCoeffsDiffer');
            fdtmp=fd;
            fdtmp.transformation='identity';
            testCase.verifyWarning(@()filter.setState(fdtmp),'setState:transDiffer');
            testCase.verifyWarning(@()filter.setState(VMDistribution(2,3)),'setState:nonFourier');
        end
        
        function testPrediction(testCase)
            % Validity of operation is tested in Fourier class, only test
            % filter interface
            for transformation={'identity','sqrt'}
                filter=FourierFilter(101,[transformation{:}]);
                vm=VMDistribution(3,2);
                warnstruct=warning('off','Predict:automaticConversion');
                filter.predictIdentity(vm);
                warning(warnstruct);
                fd1=FourierDistribution.fromDistribution(CircularUniformDistribution(),101,[transformation{:}]);
                fd2=fd1.convolve(FourierDistribution.fromDistribution(vm,101,[transformation{:}]));
                testCase.verifyLength(filter.fd.a,51);
                testCase.verifyLength(filter.fd.b,50);
                testCase.verifyEqual(filter.fd.a,fd2.a,'AbsTol',1E-8);
                testCase.verifyEqual(filter.fd.b,fd2.b,'AbsTol',1E-8);
            end
        end
        
        function testUpdate(testCase)
            % As above, only test filter interface
            for transformation={'identity','sqrt'}
                filter=FourierFilter(101,[transformation{:}]);
                vmMultiply=VMDistribution(3,2);
                vmFilter=VMDistribution(0,2);
                warnstruct=warning('off','Update:automaticConversion');
                filter.updateIdentity(vmFilter,3);
                warning(warnstruct);
                fd1=FourierDistribution.fromDistribution(CircularUniformDistribution(),101,[transformation{:}]);
                fd2=fd1.multiply(FourierDistribution.fromDistribution(vmMultiply,101,[transformation{:}]));
                fd3=fd2.truncate(101);
                testCase.verifyLength(filter.fd.a,51);
                testCase.verifyLength(filter.fd.b,50);
                testCase.verifyEqual(filter.fd.a,fd3.a,'AbsTol',1E-8);
                testCase.verifyEqual(filter.fd.b,fd3.b,'AbsTol',1E-8);
            end
        end
        
        function testPredictNonlinearForLinear(testCase)
            densityInit=VMDistribution(3,5);
            fNoiseDist=VMDistribution(0.5,10);
            noCoeffs=31;
            for transformation={'identity','sqrt'}
                fourierFilterLin=FourierFilter(noCoeffs,[transformation{:}]);
                fourierFilterLin.setState(FourierDistribution.fromDistribution(densityInit,noCoeffs,[transformation{:}]));
                fourierFilterNl=FourierFilter(noCoeffs,[transformation{:}]);
                fourierFilterNl.setState(FourierDistribution.fromDistribution(densityInit,noCoeffs,[transformation{:}]));

                warnstruct=warning('off','Predict:automaticConversion');
                fourierFilterLin.predictIdentity(fNoiseDist)
                fourierFilterNl.predictNonlinear(@(x)x,fNoiseDist)
                testCase.verifyEqual(fourierFilterLin.fd.kldNumerical(fourierFilterNl.fd),0,'AbsTol',1E-8);

                fNoiseDistShifted=fNoiseDist.shift(1);
                fourierFilterLin.predictIdentity(fNoiseDistShifted)
                fourierFilterNl.predictNonlinear(@(x)x+1,fNoiseDist)
                testCase.verifyEqual(fourierFilterLin.fd.kldNumerical(fourierFilterNl.fd),0,'AbsTol',1E-8);
                warning(warnstruct);
            end
        end
        
        function testGetEstimateMean(testCase)
            vm=VMDistribution(3,0.1);
            for transformation={'identity','sqrt'}
                filter=FourierFilter(31,[transformation{:}]);
                filter.setState(FourierDistribution.fromDistribution(vm,2*length(filter.fd.a)-1,[transformation{:}]));
                testCase.verifyEqual(filter.getEstimateMean,vm.mu,'AbsTol',1E-8);
            end
        end
    end
end

