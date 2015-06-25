classdef FourierSqrtFilterTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)           
        function testInitialization(testCase)
            filter=FourierSqrtFilter(101);
            testCase.verifyClass(filter.fd, 'FourierDistribution');
            testCase.verifyEqual(length(filter.fd.a),51);
            testCase.verifyEqual(length(filter.fd.b),50);
            testCase.verifyEqual(filter.fd.a(1),sqrt(2/pi),'AbsTol',1E-10);
        end
        
        function testSetState(testCase)
            filter=FourierSqrtFilter(101);
            filter.setState(FourierDistribution.fromDistribution(VMDistribution(2,3),101));
            testCase.verifyEqual(filter.fd,FourierDistribution.fromDistribution(VMDistribution(2,3),101));
        end
        
        function testPrediction(testCase)
            % Validity of operation is tested in Fourier class, only test
            % filter wrapper
            filter=FourierSqrtFilter(101);
            vm=VMDistribution(3,2);
            filter.predictIdentity(vm);
            fd1=FourierDistribution.fromDistribution(CUDistribution(),101,'sqrt');
            fd2=fd1.convolve(FourierDistribution.fromDistribution(vm,101,'sqrt'));
            testCase.verifyEqual(length(filter.fd.a),51);
            testCase.verifyEqual(length(filter.fd.b),50);
            testCase.verifyEqual(filter.fd.a,fd2.a,'AbsTol',1E-8);
            testCase.verifyEqual(filter.fd.b,fd2.b,'AbsTol',1E-8);
        end
        
        function testUpdate(testCase)
            % As above, only test filter wrapper
            filter=FourierSqrtFilter(101);
            vmMultiply=VMDistribution(3,2);
            vmFilter=VMDistribution(0,2);
            filter.updateIdentity(vmFilter,3);
            fd1=FourierDistribution.fromDistribution(CUDistribution(),101,'sqrt');
            fd2=fd1.multiply(FourierDistribution.fromDistribution(vmMultiply,101,'sqrt'));
            fd3=fd2.truncate(101);
            testCase.verifyEqual(length(filter.fd.a),51);
            testCase.verifyEqual(length(filter.fd.b),50);
            testCase.verifyEqual(filter.fd.a,fd3.a,'AbsTol',1E-8);
            testCase.verifyEqual(filter.fd.b,fd3.b,'AbsTol',1E-8);
        end
    end
end

