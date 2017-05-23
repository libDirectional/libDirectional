classdef HypertoroidalFourierDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testConstructor(testCase)
            hfd=HypertoroidalFourierDistribution([0,0,0;0,1/sqrt(2*pi)^2,0;0,0,0]);
            testCase.verifySize(hfd.C,[3,3]);
            testCase.verifyEqual(hfd.transformation,'sqrt');
            testCase.verifyWarning(@()HypertoroidalFourierDistribution(1/sqrt(2*pi)),'fourierCoefficients:singleCoefficient');
            testCase.verifyError(@()HypertoroidalFourierDistribution(WNDistribution(0,1)),'fourierCoefficients:invalidCoefficientMatrix');
        end
        
        function testNormalization2D(testCase)
            unnormalizedCoeffs2D=fftshift(fftn(rand(3,7)+0.5));
            unnormalizedCoeffs2D(2,4)=1;
            warningSettings=warning('off','Normalization:notNormalized');
            hfdId=HypertoroidalFourierDistribution(unnormalizedCoeffs2D,'identity');
            hfdSqrt=HypertoroidalFourierDistribution(unnormalizedCoeffs2D,'sqrt');
            warning(warningSettings);
            testCase.verifyEqual(integral2(@(x,y)reshape(hfdId.pdf([x(:)';y(:)']),size(x)),0,2*pi,0,2*pi),1,'RelTol',1E-4);
            testCase.verifyEqual(integral2(@(x,y)reshape(hfdSqrt.pdf([x(:)';y(:)']),size(x)),0,2*pi,0,2*pi),1,'RelTol',1E-4);
            % Test warnings
            testCase.verifyWarning(@()HypertoroidalFourierDistribution([0,0,0;0,-1,0;0,0,0],'identity'),'Normalization:negative');
            testCase.verifyError(@()HypertoroidalFourierDistribution([0,0,0;0,1E-201,0;0,0,0],'identity'),'Normalization:almostZero');
        end
        
        function testNormalization3D(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = false; end
            if enableExpensive
                unnormalizedCoeffs3D=fftshift(fftn(rand(3,11,7)+0.5));
                unnormalizedCoeffs3D(2,6,4)=1;
                warningSettings=warning('off','Normalization:notNormalized');
                hfdId=HypertoroidalFourierDistribution(unnormalizedCoeffs3D,'identity');
                hfdSqrt=HypertoroidalFourierDistribution(unnormalizedCoeffs3D,'sqrt');
                warning(warningSettings);
                testCase.verifyEqual(integral3(@(x,y,z)reshape(hfdId.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,0,2*pi,0,2*pi),1,'RelTol',1E-4);
                testCase.verifyEqual(integral3(@(x,y,z)reshape(hfdSqrt.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,0,2*pi,0,2*pi),1,'RelTol',1E-4);
            end
        end
        
        function testTruncation1D(testCase)
            hfd1=HypertoroidalFourierDistribution.fromDistribution(WNDistribution(1,1),101);
            hfd2=hfd1.truncate(51);
            testCase.verifySize(hfd2.C,[51,1]);
            xvals=linspace(0,2*pi,100);
            testCase.verifyEqual(hfd1.pdf(xvals),hfd2.pdf(xvals),'AbsTol',1E-8);
        end
        
        function testTruncation3D(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = false; end
            
            coeffs=[15,15,15];
            mu=[0;0;0];
            C=[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwnd=HypertoroidalWNDistribution(mu,C);
            hfdId=HypertoroidalFourierDistribution.fromDistribution(hwnd,coeffs,'identity');
            hfdSqrt=HypertoroidalFourierDistribution.fromDistribution(hwnd,coeffs,'sqrt');
            hfdIdTruncWithoutNormalization=hfdId;
            hfdIdTruncWithoutNormalization.C=hfdId.C(7:9,7:9,7:9);
            hfdSqrtTruncWithoutNormalization=hfdSqrt;
            hfdSqrtTruncWithoutNormalization.C=hfdSqrt.C(7:9,7:9,7:9);
            % Test warnings
            testCase.verifyWarning(@()hfdId.truncate([21,21,21]),'Truncate:TooFewCoefficients');
            testCase.verifyWarning(@()hfdSqrt.truncate([21,21,21]),'Truncate:TooFewCoefficients');
            
            if enableExpensive 
                % No problems exist when using identity transformation
                testCase.verifyEqual(integral3(@(x,y,z)reshape(hfdIdTruncWithoutNormalization.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,0,2*pi,0,2*pi),1,'RelTol',1E-4);
                % Show that in this example, naive truncation does not work for
                % sqrt case
                testCase.verifyTrue(abs(integral3(@(x,y,z)reshape(hfdSqrtTruncWithoutNormalization.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,0,2*pi,0,2*pi)-1)>0.01);
                % Show that .truncate handles it correctly
                hfdSqrtTrunc=hfdSqrt.truncate([3,3,3]);
                testCase.verifyEqual(integral3(@(x,y,z)reshape(hfdSqrtTrunc.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,0,2*pi,0,2*pi),1,'RelTol',1E-4);
            end
        end
        
        function testFromFunction2D(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = false; end
            if enableExpensive
                [xTest,yTest]=meshgrid(-pi:0.1:3*pi);
            else
                [xTest,yTest]=meshgrid(0:1:2*pi);
            end
            coeffs=[13,15];
            mu=[1;0];
            for kappa1=[0.2,1]
                for kappa2=[0.5,2]
                    for lambda=[0,1]
                        tvm=ToroidalVMSineDistribution(mu,[kappa1;kappa2],lambda);
                        hfdId=HypertoroidalFourierDistribution.fromFunction(...
                            @(x,y)reshape(tvm.pdf([x(:)';y(:)']),size(x)),coeffs,'identity');
                        hfdSqrt=HypertoroidalFourierDistribution.fromFunction(...
                            @(x,y)reshape(tvm.pdf([x(:)';y(:)']),size(x)),coeffs,'sqrt');
                        warnstruct=warning('off','Normalization:cannotTest');
                        hfdLog=HypertoroidalFourierDistribution.fromFunction(...
                            @(x,y)reshape(tvm.pdf([x(:)';y(:)']),size(x)),coeffs,'log');
                        warning(warnstruct);
                        testCase.verifyClass(hfdId,'HypertoroidalFourierDistribution');
                        testCase.verifyClass(hfdSqrt,'HypertoroidalFourierDistribution');
                        testCase.verifyClass(hfdLog,'HypertoroidalFourierDistribution');
                        testCase.verifySize(hfdId.C,coeffs);
                        testCase.verifySize(hfdSqrt.C,coeffs);
                        testCase.verifySize(hfdLog.C,coeffs);
                        testCase.verifyEqual(tvm.pdf([xTest(:)';yTest(:)']),hfdId.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-5);
                        testCase.verifyEqual(tvm.pdf([xTest(:)';yTest(:)']),hfdSqrt.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-6);
                        warnstruct=warning('off','pdf:mayNotBeNormalized');
                        testCase.verifyEqual(tvm.pdf([xTest(:)';yTest(:)']),hfdLog.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-6);
                        warning(warnstruct);
                    end
                end
            end
            testCase.verifyError(@()HypertoroidalFourierDistribution.fromFunction(...
                @(x,y)reshape(tvm.pdf([x(:)';y(:)']),size(x)),coeffs,'abc'),'fromFunctionValues:unrecognizedTranformation')
        end
        
        function testFromFunction3D(testCase)
            rng default
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = false; end
            if enableExpensive
                testPoints=rand(3,1000);
            else
                testPoints=rand(3,30);
            end
            C=[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu=2*pi*rand(3,1);
            hwnd=HypertoroidalWNDistribution(mu,C);
            coeffs=[25,27,23];
            hfdId=HypertoroidalFourierDistribution.fromFunction(...
                @(x,y,z)reshape(hwnd.pdf([x(:)';y(:)';z(:)']),size(x)),coeffs,'identity');
            hfdSqrt=HypertoroidalFourierDistribution.fromFunction(...
                @(x,y,z)reshape(hwnd.pdf([x(:)';y(:)';z(:)']),size(x)),coeffs,'sqrt');
            testCase.verifyClass(hfdId,'HypertoroidalFourierDistribution')
            testCase.verifySize(hfdId.C,coeffs);
            testCase.verifyEqual(hfdId.pdf(testPoints),hwnd.pdf(testPoints),'AbsTol',1E-6);
            testCase.verifyEqual(hfdSqrt.pdf(testPoints),hwnd.pdf(testPoints),'AbsTol',1E-5);
        end
        
        function testFromFunction4D(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = false; end
            if enableExpensive
                C=[0.7,0.4,0.2,-0.5;0.4,0.6,0.1,0;0.2,0.1,1,-0.3;-0.5,0,-0.3,0.9]*2;
                mu=2*pi*rand(4,1);
                hwnd=HypertoroidalWNDistribution(mu,C);
                coeffs=[19,17,15,21];
                hfdId=HypertoroidalFourierDistribution.fromFunction(...
                    @(x,y,z,w)reshape(hwnd.pdf([x(:)';y(:)';z(:)';w(:)']),size(x)),coeffs,'identity');
                hfdSqrt=HypertoroidalFourierDistribution.fromFunction(...
                    @(x,y,z,w)reshape(hwnd.pdf([x(:)';y(:)';z(:)';w(:)']),size(x)),coeffs,'sqrt');
                testCase.verifyClass(hfdId,'HypertoroidalFourierDistribution')
                testCase.verifyClass(hfdSqrt,'HypertoroidalFourierDistribution')
                testCase.verifySize(hfdId.C,coeffs);
                testCase.verifySize(hfdSqrt.C,coeffs);
                testPoints=rand(4,100);
                testCase.verifyEqual(hfdId.pdf(testPoints),hwnd.pdf(testPoints,4),'AbsTol',1E-6);
                testCase.verifyEqual(hfdSqrt.pdf(testPoints),hwnd.pdf(testPoints,4),'AbsTol',1E-5);
            end
        end
        
        function testFromDistribution2D(testCase)
            % Test that from Distribution and fromFunction result in equal
            % approximations
            kappa1=0.3;kappa2=1.5;
            lambda=0.5;
            coeffs=[5,7];
            tvm=ToroidalVMSineDistribution([1;2],[kappa1;kappa2],lambda);
            hfd1id=HypertoroidalFourierDistribution.fromFunction(@(x,y)reshape(tvm.pdf([x(:)';y(:)']),size(x)),coeffs,'identity');
            hfd1sqrt=HypertoroidalFourierDistribution.fromFunction(@(x,y)reshape(tvm.pdf([x(:)';y(:)']),size(x)),coeffs,'sqrt');
            hfd2id=HypertoroidalFourierDistribution.fromDistribution(tvm,coeffs,'identity');
            hfd2sqrt=HypertoroidalFourierDistribution.fromDistribution(tvm,coeffs,'sqrt');
            testCase.verifyClass(hfd2id,'HypertoroidalFourierDistribution');
            testCase.verifyClass(hfd2sqrt,'HypertoroidalFourierDistribution');
            testCase.verifySize(hfd2id.C,coeffs);
            testCase.verifySize(hfd2sqrt.C,coeffs);
            % Verify approximation by validating coefficients
            testCase.verifyEqual(hfd2id.C,hfd1id.C,'AbsTol',1E-10);
            testCase.verifyEqual(hfd2sqrt.C,hfd1sqrt.C,'AbsTol',1E-10);
        end
        
        function testFromDistribution3D(testCase)
            % Test that error is thrown if dimensionality is wrong
            testCase.verifyError(...
                @()HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution([1;1;1],eye(3)),[15,15]),...
                'fromDistribution:invalidObject');
            % Test that it yields same results as fromFunction
            C=[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu=[3;5;2];
            coeffs=[21,21,21];
            hwnd=HypertoroidalWNDistribution(mu,C);
            hfd1id=HypertoroidalFourierDistribution.fromFunction(@(x,y,z)reshape(hwnd.pdf([x(:)';y(:)';z(:)']),size(x)),coeffs,'identity');
            hfd1sqrt=HypertoroidalFourierDistribution.fromFunction(@(x,y,z)reshape(hwnd.pdf([x(:)';y(:)';z(:)']),size(x)),coeffs,'sqrt');
            hfd2id=HypertoroidalFourierDistribution.fromDistribution(hwnd,coeffs,'identity');
            hfd2sqrt=HypertoroidalFourierDistribution.fromDistribution(hwnd,coeffs,'sqrt');
            testCase.verifyClass(hfd2id,'HypertoroidalFourierDistribution');
            testCase.verifyClass(hfd2sqrt,'HypertoroidalFourierDistribution');
            testCase.verifySize(hfd2id.C,coeffs);
            testCase.verifySize(hfd2sqrt.C,coeffs);
            testCase.verifyEqual(hfd2id.C,hfd1id.C,'AbsTol',1E-10);
            testCase.verifyEqual(hfd2sqrt.C,hfd1sqrt.C,'AbsTol',1E-10);
        end
        
        function testFromDistributionHWN(testCase)
            % Verify closed form solution against FFT version
            % 2D case
            mu=[1;2];
            C=2*[1,0.5;0.5,1.3];
            hwn=HypertoroidalWNDistribution(mu,C);
            hfd1id=HypertoroidalFourierDistribution.fromDistribution(hwn,[35,35],'identity');
            hfd2id=HypertoroidalFourierDistribution.fromFunction(@(x,y)reshape(hwn.pdf([x(:)';y(:)']),size(x)),[35,35],'identity');
            testCase.verifyEqual(hfd1id.C,hfd2id.C,'AbsTol',1E-7);
            % 3D case
            mu=[1;2;4];
            C=[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn=HypertoroidalWNDistribution(mu,C);
            hfd1id=HypertoroidalFourierDistribution.fromDistribution(hwn,[19,19,19],'identity');
            hfd2id=HypertoroidalFourierDistribution.fromFunction(@(x,y,z)reshape(hwn.pdf([x(:)';y(:)';z(:)']),size(x)),[19,19,19],'identity');
            testCase.verifyEqual(hfd1id.C,hfd2id.C,'AbsTol',1E-7);
        end
        
        function testMultiply2DWithoutTruncation(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = false; end
            if enableExpensive
                [xTest,yTest]=meshgrid(-pi:0.1:3*pi);
            else
                [xTest,yTest]=meshgrid(0:1:2*pi);
            end
            warning('off', 'pdf:mayNotBeNormalized');
            warning('off', 'Normalization:cannotTest');
            
            tvm1=ToroidalVMSineDistribution([1;3],[0.3;0.5],0.5);
            tvm2=ToroidalVMSineDistribution([1;4],[0.8;1.5],0.2);
            
            hfd1id=HypertoroidalFourierDistribution.fromDistribution(tvm1,[17,15],'identity');
            hfd2id=HypertoroidalFourierDistribution.fromDistribution(tvm2,[15,17],'identity');
            hfd1sqrt=HypertoroidalFourierDistribution.fromDistribution(tvm1,[17,15],'sqrt');
            hfd2sqrt=HypertoroidalFourierDistribution.fromDistribution(tvm2,[15,17],'sqrt');
            hfd1log=HypertoroidalFourierDistribution.fromDistribution(tvm1,[17,15],'log');
            hfd2log=HypertoroidalFourierDistribution.fromDistribution(tvm2,[15,17],'log');
            hfdMultid=hfd1id.multiply(hfd2id,[17,15]+[15,17]-1);
            hfdMultsqrt=hfd1sqrt.multiply(hfd2sqrt,[17,15]+[15,17]-1);
            
            warnstruct=[warning('off', 'Truncate:TooFewCoefficients'),warning('off', 'Multiply:NotNormalizing')];
            hfdMultlog=hfd1log.multiply(hfd2log,[17,17]);
            warning(warnstruct);
            
            testCase.verifyClass(hfdMultid,'HypertoroidalFourierDistribution');
            testCase.verifyClass(hfdMultsqrt,'HypertoroidalFourierDistribution');
            testCase.verifyClass(hfdMultlog,'HypertoroidalFourierDistribution');
            % Verify that they integrate to 1 (except log)
            testCase.verifyEqual(...
                integral2(@(x,y)reshape(hfdMultid.pdf([x(:)';y(:)']),size(x)),0,2*pi,0,2*pi),1,'AbsTol',1E-6);
            testCase.verifyEqual(...
                integral2(@(x,y)reshape(hfdMultsqrt.pdf([x(:)';y(:)']),size(x)),0,2*pi,0,2*pi),1,'AbsTol',1E-6);
            % Normalize multiplication using integral to obtain approximate
            % ground truth
            normConst=integral2(@(x,y)reshape(tvm1.pdf([x(:)';y(:)']).*tvm2.pdf([x(:)';y(:)']),size(x)),0,2*pi,0,2*pi);
            valTrueApprox=reshape(tvm1.pdf([xTest(:)';yTest(:)']).*tvm2.pdf([xTest(:)';yTest(:)']),size(xTest))/normConst;
            valFourierId=reshape(hfdMultid.pdf([xTest(:)';yTest(:)']),size(xTest));
            valFourierSqrt=reshape(hfdMultsqrt.pdf([xTest(:)';yTest(:)']),size(xTest));
            valFourierLog=reshape(hfdMultlog.pdf([xTest(:)';yTest(:)']),size(xTest));
            % Verify correct function values
            testCase.verifyEqual(valFourierId,valTrueApprox,'AbsTol',1E-6);
            testCase.verifyEqual(valFourierSqrt,valTrueApprox,'AbsTol',1E-7);
            normConstLog=1/integral2(@(x,y)reshape(hfdMultlog.pdf([x(:)';y(:)']),size(x)),0,2*pi,0,2*pi);
            testCase.verifyEqual(valFourierLog*normConstLog,valTrueApprox,'AbsTol',1E-6);
            
            warning('on', 'pdf:mayNotBeNormalized');
            warning('on', 'Normalization:cannotTest');
        end
        
        function testMultiply3DWithoutTruncation(testCase)
            C=[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu=[1;2;3];
            hwnd1=HypertoroidalWNDistribution(mu,C);
            C=[1.6,0.8,1.3;0.8,0.7,0.6;1.3,0.6,1.1];
            hwnd2=HypertoroidalWNDistribution(mu,C);
            coeffs=19*ones(1,3);
            
            hfd1id=HypertoroidalFourierDistribution.fromFunction(...
                @(x,y,z)reshape(hwnd1.pdf([x(:)';y(:)';z(:)']),size(x)),coeffs,'identity');
            hfd2id=HypertoroidalFourierDistribution.fromFunction(...
                @(x,y,z)reshape(hwnd2.pdf([x(:)';y(:)';z(:)']),size(x)),coeffs,'identity');
            hfdResultId=hfd1id.multiply(hfd2id,2*coeffs-1);
            
            [x,y,z]=meshgrid(linspace(-4*pi,pi,5));
            testPoints=[x(:)';y(:)';z(:)'];
            % Calculate approximate normalization factor because integral3
            % is expensive.
            approxNormFactorId=hfdResultId.pdf([1;1;1])/(hfd1id.pdf([1;1;1])*hfd2id.pdf([1;1;1]));
            testCase.verifyEqual(approxNormFactorId*(hfd1id.pdf(testPoints).*hfd2id.pdf(testPoints)),hfdResultId.pdf(testPoints),'AbsTol',1E-10);
            
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = false; end
            if enableExpensive            
                hfd1sqrt=HypertoroidalFourierDistribution.fromFunction(...
                    @(x,y,z)reshape(hwnd1.pdf([x(:)';y(:)';z(:)']),size(x)),coeffs,'sqrt');
                hfd2sqrt=HypertoroidalFourierDistribution.fromFunction(...
                    @(x,y,z)reshape(hwnd2.pdf([x(:)';y(:)';z(:)']),size(x)),coeffs,'sqrt');                
                hfdResultSqrt=hfd1sqrt.multiply(hfd2sqrt,2*coeffs-1);
                approxNormFactorSqrt=hfdResultSqrt.pdf([1;2;1])/(hfd1sqrt.pdf([1;2;1])*hfd2sqrt.pdf([1;2;1]));
                testCase.verifyEqual(approxNormFactorSqrt*(hfd1sqrt.pdf(testPoints).*hfd2sqrt.pdf(testPoints)),hfdResultSqrt.pdf(testPoints),'AbsTol',1E-10);
            end
        end
        
        function testConvolve2D(testCase)
            tvm1=ToroidalVMSineDistribution([1;2],[0.3;0.5],0.5);
            tvm2=ToroidalVMSineDistribution([1;4],[0.8;1.5],0.2);
            hfd1id=HypertoroidalFourierDistribution.fromDistribution(tvm1,[17,15],'identity');
            hfdtid=HypertoroidalFourierDistribution.fromDistribution(tvm2,[15,17],'identity');
            hfd1sqrt=HypertoroidalFourierDistribution.fromDistribution(tvm1,[17,15],'sqrt');
            hfdtsqrt=HypertoroidalFourierDistribution.fromDistribution(tvm2,[15,17],'sqrt');
            warnstruct=warning('off','Normalization:cannotTest');
            hfd1log=HypertoroidalFourierDistribution.fromDistribution(tvm1,[17,15],'log');
            hfdtlog=HypertoroidalFourierDistribution.fromDistribution(tvm2,[15,17],'log');
            testCase.verifyError(@()hfd1log.convolve(hfdtlog,[15,13]),'transformation:unrecognizedTransformation');
            warning(warnstruct);
            
            hfd2id=hfd1id.convolve(hfdtid,[15,13]);
            hfd2sqrt=hfd1sqrt.convolve(hfdtsqrt,[15,13]);
            testCase.verifyClass(hfd2id,'HypertoroidalFourierDistribution')
            testCase.verifyClass(hfd2sqrt,'HypertoroidalFourierDistribution')
            testCase.verifySize(hfd2id.C,[15,13])
            testCase.verifySize(hfd2sqrt.C,[15,13])
            hfd2id=hfd1id.convolve(hfdtid);
            hfd2sqrt=hfd1sqrt.convolve(hfdtsqrt);
            testCase.verifySize(hfd2id.C,[17,15]);
            testCase.verifySize(hfd2sqrt.C,[17,15]);
            % Verify that they integrate to 1
            testCase.verifyEqual(...
                integral2(@(x,y)reshape(hfd2id.pdf([x(:)';y(:)']),size(x)),0,2*pi,0,2*pi),1,'AbsTol',1E-6);
            testCase.verifyEqual(...
                integral2(@(x,y)reshape(hfd2sqrt.pdf([x(:)';y(:)']),size(x)),0,2*pi,0,2*pi),1,'AbsTol',1E-6);
            [xTest,yTest]=meshgrid(0:2*pi/100:2*pi-2*pi/100);
            valFourierId=reshape(hfd2id.pdf([xTest(:)';yTest(:)']),size(xTest));
            valFourierSqrt=reshape(hfd2sqrt.pdf([xTest(:)';yTest(:)']),size(xTest));
            % Calculate cyclic discrete convolution, the result is not 
            % normalized, so use values of Fourier for normalization
            % Cyclic convolution would be easier with FFT & IFFT but since
            % is approach is used in FourierDistribution, a different
            % approach is used for validation
            tvm1vals=reshape(tvm1.pdf([xTest(:)';yTest(:)']),size(xTest));
            tvm2valspadded=padarray(reshape(tvm2.pdf([xTest(:)';yTest(:)']),size(xTest)),size(tvm1vals),'circular');
            convolutionTmp=conv2(tvm1vals,tvm2valspadded);
            valConvUnnorm=convolutionTmp(size(tvm1vals,1)+1:2*size(tvm1vals,1),size(tvm1vals,2)+1:2*size(tvm1vals,2));
            testCase.verifyEqual(valConvUnnorm/sum(sum(valConvUnnorm))*sum(sum(valFourierId)),valFourierId,'AbsTol',1E-6);
            testCase.verifyEqual(valConvUnnorm/sum(sum(valConvUnnorm))*sum(sum(valFourierSqrt)),valFourierSqrt,'AbsTol',1E-6);
        end
        function testShift(testCase)
            Ccell={2*[1,0.5;0.5,1],[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1],...
                [0.7,0.4,0.2,-0.5;0.4,0.6,0.1,0;0.2,0.1,1,-0.3;-0.5,0,-0.3,0.9]*2};
            offsets={4*pi*rand(2,1)-pi,4*pi*rand(3,1)-pi,4*pi*rand(4,1)-pi};
            coeffsPerDim=13;
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = false; end
            if enableExpensive
                maxDim = 4;
            else
                maxDim = 3;
            end
            for dim=2:maxDim
                hfdId=HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(zeros(dim,1),Ccell{dim-1}),coeffsPerDim*ones(1,dim),'identity');
                hfdSqrt=HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(zeros(dim,1),Ccell{dim-1}),coeffsPerDim*ones(1,dim),'sqrt');
                hfdIdShiftWN=HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(offsets{dim-1},Ccell{dim-1}),coeffsPerDim*ones(1,dim),'identity');
                hfdSqrtShiftWN=HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(offsets{dim-1},Ccell{dim-1}),coeffsPerDim*ones(1,dim),'sqrt');
                hfdIdShiftFD=hfdId.shift(offsets{dim-1});
                hfdSqrtShiftFD=hfdSqrt.shift(offsets{dim-1});
                testCase.verifyEqual(hfdIdShiftFD.C,hfdIdShiftWN.C,'AbsTol',abs(max(hfdIdShiftWN.C(:))/1000));
                testCase.verifyEqual(hfdSqrtShiftFD.C,hfdSqrtShiftWN.C,'AbsTol',abs(max(hfdSqrtShiftWN.C(:))/1000));
            end
        end
        function testIntegral2D(testCase)
            % Test against implementation in toroidal (test case that this
            % works correctly exists in ToroidalFourierDistributionTest)
            kappa1=0.3;kappa2=1.5;
            lambda=0.5;
            coeffs=[5,7];
            tvm=ToroidalVMSineDistribution([1;2],[kappa1;kappa2],lambda);
            hfdId=HypertoroidalFourierDistribution.fromFunction(@(x,y)reshape(tvm.pdf([x(:)';y(:)']),size(x)),coeffs,'identity');
            hfdSqrt=HypertoroidalFourierDistribution.fromFunction(@(x,y)reshape(tvm.pdf([x(:)';y(:)']),size(x)),coeffs,'sqrt');

            tfdId=ToroidalFourierDistribution(hfdId.C,hfdId.transformation);
            tfdSqrt=ToroidalFourierDistribution(hfdSqrt.C,hfdSqrt.transformation);
            testCase.verifyEqual(hfdId.integral([0;0],[pi;pi]),tfdId.integral([0;0],[pi;pi]),'AbsTol',1E-4);
            testCase.verifyEqual(hfdSqrt.integral([0;0],[pi;pi]),tfdSqrt.integral([0;0],[pi;pi]),'AbsTol',1E-4);
            testCase.verifyEqual(hfdId.integral([0;0],[pi;2*pi]),tfdId.integral([0;0],[pi;2*pi]),'AbsTol',1E-4);
            testCase.verifyEqual(hfdSqrt.integral([0;0],[pi;2*pi]),tfdSqrt.integral([0;0],[pi;2*pi]),'AbsTol',1E-4);
            testCase.verifyEqual(hfdId.integral([0;-1],[3*pi;5*pi]),tfdId.integral([0;-1],[3*pi;5*pi]),'AbsTol',1E-4);
            testCase.verifyEqual(hfdSqrt.integral([0;-1],[3*pi;5*pi]),tfdSqrt.integral([0;-1],[3*pi;5*pi]),'AbsTol',1E-4);
        end
        function testIntegral3D(testCase)
            % Test against integral
            C=[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu=[1;2;5];
            coeffs=[5,9,7];
            warningSettings=warning('off','Normalization:notNormalized'); % Very coarse approximation, lack of normalization is to be expected
            hfdId=HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(mu,C),coeffs,'identity');
            hfdSqrt=HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(mu,C),coeffs,'sqrt');
            warning(warningSettings);
            
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = false; end
            if enableExpensive
                testCase.verifyEqual(hfdId.integral([0;0;0], [pi;pi;pi]),integral3(@(x,y,z)reshape(hfdId.pdf([x(:)';y(:)';z(:)']),size(x)),0,pi,0,pi,0,pi),'AbsTol',1E-4);
                testCase.verifyEqual(hfdSqrt.integral([0;0;0], [pi;pi;pi]),integral3(@(x,y,z)reshape(hfdSqrt.pdf([x(:)';y(:)';z(:)']),size(x)),0,pi,0,pi,0,pi),'AbsTol',1E-4);
            end
        end
        function testCovariance2dimd2D(testCase)
            % Test against integral
            C=[0.7,0.4;0.4,0.6];
            mu=[1;2];
            coeffs=[15,17];
            twn=ToroidalWNDistribution(mu,C);
            hfdId=HypertoroidalFourierDistribution.fromDistribution(twn,coeffs,'identity');
            hfdSqrt=HypertoroidalFourierDistribution.fromDistribution(twn,coeffs,'sqrt');
            
            testCase.verifyEqual(hfdId.covariance2dimD,twn.covariance4D,'AbsTol',1E-4);
            testCase.verifyEqual(hfdSqrt.covariance2dimD,twn.covariance4D,'AbsTol',1E-4);
        end
        function testCovariance2dimd3D(testCase)
            C=[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu=[1;2;5];
            coeffs=[15,19,17];
            hfdId=HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(mu,C),coeffs,'identity');
            hfdSqrt=HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(mu,C),coeffs,'sqrt');
            
            % Test by taking only two correlated dimensions and comparing
            % to that result for each pair of dimensions
            Cov6DId=hfdId.covariance2dimD;
            Cov6DSqrt=hfdSqrt.covariance2dimD;
            twn=ToroidalWNDistribution(mu(1:2),C(1:2,1:2));
            testCase.verifyEqual(Cov6DId(1:4,1:4),twn.covariance4D,'AbsTol',1E-4);
            testCase.verifyEqual(Cov6DSqrt(1:4,1:4),twn.covariance4D,'AbsTol',1E-4);
            twn=ToroidalWNDistribution(mu([1,3]),C([1,3],[1,3]));
            testCase.verifyEqual(Cov6DId([1,2,5,6],[1,2,5,6]),twn.covariance4D,'AbsTol',1E-4);
            testCase.verifyEqual(Cov6DSqrt([1,2,5,6],[1,2,5,6]),twn.covariance4D,'AbsTol',1E-4);
            twn=ToroidalWNDistribution(mu(2:3),C(2:3,2:3));
            testCase.verifyEqual(Cov6DId(3:6,3:6),twn.covariance4D,'AbsTol',1E-4);
            testCase.verifyEqual(Cov6DSqrt(3:6,3:6),twn.covariance4D,'AbsTol',1E-4);
        end
        function testMarginals2D(testCase)
            tvm=ToroidalVMSineDistribution([1;2],[2;3],3);
            hfd=HypertoroidalFourierDistribution.fromDistribution(tvm,[53,53],'identity');
            for d=1:2
                vm=tvm.marginalizeTo1D(d);
                fd1=hfd.marginalizeTo1D(d);
                fd2=hfd.marginalizeOut(1+(d==1));
                testCase.verifyEqual(vm.hellingerDistanceNumerical(fd1),0,'AbsTol',1E-3);
                testCase.verifyEqual(vm.hellingerDistanceNumerical(fd2),0,'AbsTol',1E-3);
            end
        end
        function testPlotting(testCase)
            dist=ToroidalVMSineDistribution([1;2],[2;3],3);
            tfd=ToroidalFourierDistribution.fromDistribution(dist,[51,201],'identity');
            figure(987)
            h=tfd.plot;
            noPoints=102;
            [alpha,beta]=meshgrid(linspace(0,2*pi,noPoints));
            fvals=reshape(tfd.pdf([alpha(:)'; beta(:)']),size(alpha));
            testCase.verifyEqual(get(h,'ZData'),fvals,'AbsTol',1E-15);
            close(987)
        end
        function testConv1D(testCase)
            for transformation={'identity','sqrt'}
                fd1=FourierDistribution.fromDistribution(VMDistribution(0,1),13,[transformation{:}]);
                fd2=FourierDistribution.fromDistribution(VMDistribution(2,1),13,[transformation{:}]);
                hfd1=HypertoroidalFourierDistribution(fd1.c',[transformation{:}]);
                hfd2=HypertoroidalFourierDistribution(fd2.c',[transformation{:}]);
                fdConv=fd1.convolve(fd2);
                hfdConv=hfd1.convolve(hfd2);
                testCase.verifyEqual(hfdConv.C,fdConv.c','AbsTol',1E-8);
            end
        end
        function testMult1D(testCase)
            for transformation={'identity','sqrt'}
                fd1=FourierDistribution.fromDistribution(VMDistribution(0,1),13,[transformation{:}]);
                fd2=FourierDistribution.fromDistribution(VMDistribution(2,1),13,[transformation{:}]);
                hfd1=HypertoroidalFourierDistribution(fd1.c',[transformation{:}]);
                hfd2=HypertoroidalFourierDistribution(fd2.c',[transformation{:}]);
                fdMult=fd1.multiply(fd2,2*(numel(fd1.a)+numel(fd1.b))-1);
                hfdMult=hfd1.multiply(hfd2,2*size(hfd1.C,1)-1);
                testCase.verifyEqual(hfdMult.C,fdMult.c','AbsTol',1E-10);
            end
        end
    end
end
