classdef SphericalHarmonicDistributionRealTest< matlab.unittest.TestCase    
    methods (Test)
        function testNormalization(testCase)
            testCase.verifyWarning(@()SphericalHarmonicDistributionReal(1),'Normalization:notNormalized');
            testCase.verifyError(@()SphericalHarmonicDistributionReal(0),'Normalization:almostZero');
            warningSettings=warning('off','Normalization:notNormalized');
            unnormalizedCoeffs=rand(3,5);
            shd=SphericalHarmonicDistributionReal(unnormalizedCoeffs);
            warning(warningSettings);
            testCase.verifyEqual(shd.integral,1,'AbsTol',1E-6);
            % Enforce unnormalized coefficients and compare ratio
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            valsNormalized=shd.pdf([x;y;z]);
            shd.coeffMat=unnormalizedCoeffs;
            valsUnnormalized=shd.pdf([x;y;z]);
            testCase.verifyEqual(diff(valsNormalized./valsUnnormalized),zeros(1,size(x,2)-1),'AbsTol',1E-6);
        end
        function testTruncation(testCase)
            warningSettings=warning('off','Normalization:notNormalized');
            rshd=SphericalHarmonicDistributionReal(rand(4,7));
            warning(warningSettings);
            testCase.verifyWarning(@()rshd.truncate(4),'Truncate:TooFewCoefficients');
            warningSettings=[warningSettings,warning('off','Truncate:TooFewCoefficients')];
            rshd2=rshd.truncate(4);
            testCase.verifySize(rshd2.coeffMat,[5,9]);
            testCase.verifyTrue(all(isnan(rshd2.coeffMat(5,:))|rshd2.coeffMat(5,:)==0));
            rshd3=rshd.truncate(5);
            warning(warningSettings);
            testCase.verifySize(rshd3.coeffMat,[6,11]);
            testCase.verifyTrue(all(all(isnan(rshd3.coeffMat(5:6,:))|rshd3.coeffMat(5:6,:)==0)));
            rshd4=rshd2.truncate(3);
            testCase.verifySize(rshd4.coeffMat,[4,7]);
            rshd5=rshd3.truncate(3);
            testCase.verifySize(rshd5.coeffMat,[4,7]);
            
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(rshd2.pdf([x;y;z]),rshd.pdf([x;y;z]),'AbsTol',1E-6);
            testCase.verifyEqual(rshd3.pdf([x;y;z]),rshd.pdf([x;y;z]),'AbsTol',1E-6);
            testCase.verifyEqual(rshd4.pdf([x;y;z]),rshd.pdf([x;y;z]),'AbsTol',1E-6);
            testCase.verifyEqual(rshd5.pdf([x;y;z]),rshd.pdf([x;y;z]),'AbsTol',1E-6);
        end
        % Test some basis functions by comparing them with the Table 4 in
        % Chapter 4 of "Group Theoretical Techniques in Quantum Chemistry" 
        % by C. D. H. Chrisholm (some factors adjusted)
        function testl0m0(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi)); % Initialize with 1/sqrt(4*pi) and overwrite to prevent normalization
            shd.coeffMat=[1,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN;0,0,0,0,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),ones(1,10)*sqrt(1/(4*pi)),'AbsTol',1E-6);
        end
        function testl1mneg1(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN;1,0,0,NaN,NaN;0,0,0,0,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),sqrt(3/(4*pi))*y,'AbsTol',1E-6);
        end
        function testl1m0(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN;0,1,0,NaN,NaN;0,0,0,0,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),sqrt(3/(4*pi))*z,'AbsTol',1E-6);
        end
        function testl1m1(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN;0,0,1,NaN,NaN;0,0,0,0,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),sqrt(3/(4*pi))*x,'AbsTol',1E-6);
        end
        function testl2mneg2(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN;1,0,0,0,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),1/2*sqrt(15/pi)*x.*y,'AbsTol',1E-6);
        end
        function testl2mneg1(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN;0,1,0,0,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),1/2*sqrt(15/pi)*y.*z,'AbsTol',1E-6);
        end
        function testl2m0(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN;0,0,1,0,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),1/4*sqrt(5/pi)*(2*z.^2-x.^2-y.^2),'AbsTol',1E-6);
        end
        function testl2m1(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN;0,0,0,1,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),1/2*sqrt(15/pi)*x.*z,'AbsTol',1E-6);
        end
        function testl2m2(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN;0,0,0,0,1];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),1/4*sqrt(15/pi)*(x.^2-y.^2),'AbsTol',1E-6);
        end
        function testl3mneg3(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;1,0,0,0,0,0,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),1/4*sqrt(35/(2*pi))*y.*(3*x.^2-y.^2),'AbsTol',1E-6);
        end
        function testl3mneg2(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;0,1,0,0,0,0,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),1/2*sqrt(105/pi)*x.*y.*z,'AbsTol',1E-6);
        end
        function testl3mneg1(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;0,0,1,0,0,0,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),1/4*sqrt(21/(2*pi))*y.*(4*z.^2-x.^2-y.^2),'AbsTol',1E-6);
        end
        function testl3m0(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;0,0,0,1,0,0,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),1/4*sqrt(7/pi)*(z.*(2*z.^2-3*x.^2-3*y.^2)),'AbsTol',1E-6);
        end
        function testl3m1(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;0,0,0,0,1,0,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),1/4*sqrt(21/(2*pi))*x.*(4*z.^2-x.^2-y.^2),'AbsTol',1E-6);
        end
        function testl3m2(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;0,0,0,0,0,1,0];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),1/4*sqrt(105/pi)*z.*(x.^2-y.^2),'AbsTol',1E-6);
        end
        function testl3m3(testCase)
            shd=SphericalHarmonicDistributionReal(1/sqrt(4*pi));
            shd.coeffMat=[0,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;0,0,0,0,0,0,1];
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(shd.pdf([x;y;z]),1/4*sqrt(35/(2*pi))*x.*(x.^2-3*y.^2),'AbsTol',1E-6);
        end
        % Test conversion to ComplexSphericalHarmonic
        function testl0m0conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN;0,0,0,0,0]);
        end
        function testl1mneg1conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN;1,0,0,NaN,NaN;0,0,0,0,0]);
        end
        function testl1m0conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN;0,1,0,NaN,NaN;0,0,0,0,0]);
        end
        function testl1m1conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN;0,0,1,NaN,NaN;0,0,0,0,0]);
        end
        function testl2mneg2conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN;1,0,0,0,0]);
        end
        function testl2mneg1conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN;0,1,0,0,0]);
        end
        function testl2m0conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN;0,0,1,0,0]);
        end
        function testl2m1conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN;0,0,0,1,0]);
        end
        function testl2m2conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN;0,0,0,0,1]);
        end
        function testl3mneg3conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;1,0,0,0,0,0,0]);
        end
        function testl3mneg2conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;0,1,0,0,0,0,0]);
        end
        function testl3mneg1conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;0,0,1,0,0,0,0]);
        end
        function testl3m0conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;0,0,0,1,0,0,0]);
        end
        function testl3m1conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;0,0,0,0,1,0,0]);
        end
        function testl3m2conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;0,0,0,0,0,1,0]);
        end
        function testl3m3conversion(testCase)
            testConversion(testCase,[1,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN;0,0,0,0,0,0,1]);
        end
        function testrandomconversion(testCase)
            testConversion(testCase,rand(4,7));
        end
        function testConversionToComplexAndBack(testCase)
            warningSettings=warning('off','Normalization:notNormalized');
            rshd=SphericalHarmonicDistributionReal(rand(4,7));
            warning(warningSettings);
            cshd=rshd.toSphericalHarmonicDistributionComplex;
            rshd2=cshd.toSphericalHarmonicDistributionReal;
            testCase.verifyEqual(rshd2.coeffMat,rshd.coeffMat,'AbsTol',1E-6);
        end
        function integralAnalytical(testCase)
            unnormalizedCoeffs=rand(3,5);
            warningSettings=warning('off','Normalization:notNormalized');
            shd=SphericalHarmonicDistributionReal(unnormalizedCoeffs);
            warning(warningSettings);
            testCase.verifyEqual(shd.integralAnalytical,shd.integral,'AbsTol',1E-6);
        end
        
        function testRotationGrid(testCase)
            % Test rotation by verifying function values at identically
            % transformed points.
            rng default
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = false; end
            if enableExpensive           
                gridpoints = 7;
            else
                gridpoints = 3;
            end
            [phi,theta]=meshgrid(linspace(0,2*pi,30),linspace(-pi/2,pi/2,30));
            [x,y,z]=sph2cart(phi(:)',theta(:)',1);
            Rz=@(alpha)[cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1];
            Ry=@(beta)[cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)];
            
            [alphas,betas,gammas]=ndgrid(linspace(0,2*pi,gridpoints));
            for degree=7:10
                for i=1:numel(alphas)
                    coeffMat=rand(7,13);
                    coeffMat(1)=1/sqrt(4*pi);
                    shd=SphericalHarmonicDistributionReal(coeffMat);
                    shd2=shd.rotate(alphas(i),betas(i),gammas(i));
                    rotm=Rz(gammas(i))*Ry(betas(i))*Rz(alphas(i));
                    xyzRotated=rotm*[x;y;z];
                    testCase.verifyEqual(shd.pdf([x;y;z]),shd2.pdf(xyzRotated),'AbsTol',1E-6);
                end
            end
        end
    end
    methods
        function testConversion(testCase,coeffMat)
            warningSettings=warning('off','Normalization:notNormalized');
            rshd=SphericalHarmonicDistributionReal(coeffMat);
            warning(warningSettings);
            cshd=rshd.toSphericalHarmonicDistributionComplex;
            [phi,theta]=deal(rand(1,10)*2*pi,rand(1,10)*pi-pi/2);
            [x,y,z]=sph2cart(phi,theta,1);
            testCase.verifyEqual(cshd.pdf([x;y;z]),rshd.pdf([x;y;z]),'AbsTol',1E-6);
        end
    end
end
