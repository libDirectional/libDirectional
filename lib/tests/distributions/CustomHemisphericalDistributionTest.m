classdef CustomHemisphericalDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testSimpleDistribution2D(testCase)
            M = eye(3);
            Z = [-2 -0.5 0]';
            bd = BinghamDistribution(Z,M);
            chhd = CustomHemisphericalDistribution.fromDistribution(bd);
            % Use nonrandom numbers but do not influece seed everywhere
            % else
            rs = RandStream('twister','Seed',10);
            points = rs.randn(3,100);
            points(:,points(3,:)<0)=[];
            points = points./vecnorm(points);
            
            % Since it is symmetric, hemispherical version should have
            % twice the probability density than spherical one.
            testCase.verifyEqual(chhd.pdf(points) , 2*bd.pdf(points),'AbsTol',1e-5);
        end
        function testIntegralBinghamS2(testCase)
            M = eye(3);
            Z = [-2 -0.5 0]';
            bd = BinghamDistribution(Z,M);
            chhd = CustomHemisphericalDistribution.fromDistribution(bd);

            testCase.verifyEqual(chhd.integralNumerical,1,'AbsTol',1e-4);
        end
        
        function testWarningAsymmetric(testCase)
            vmf = VMFDistribution([0;0;1],10);
            testCase.verifyWarning(...
                @()CustomHemisphericalDistribution.fromDistribution(vmf),...
                'FromDistribution:UsePdfHypersphere');
        end
        
        function testIntegralVmfS2Manually(testCase)
            % No scaling neeeded because we only integrate over the upper
            % hemisphere
            vmf = VMFDistribution([1;0;0],1);
            hud = CustomHemisphericalDistribution(@(x)vmf.pdf(x)+vmf.pdf(-x));
            testCase.verifyEqual(hud.integralNumerical,1,'AbsTol',5e-7);
        end
        
        function testIntegralVmfS2FromDistribution(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vmf1 = VMFDistribution([0;0;1],10);
            vmf2 = VMFDistribution(1/sqrt(2)*[0;1;1],10);
            vmf3 = VMFDistribution([0;1;0],10); % This is antipodally symmetric
            fixture = testCase.applyFixture(SuppressedWarningsFixture('FromDistribution:UsePdfHypersphere'));
            chhd1 = CustomHemisphericalDistribution.fromDistribution(vmf1);
            chhd2 = CustomHemisphericalDistribution.fromDistribution(vmf2);
            chhd3 = CustomHemisphericalDistribution.fromDistribution(vmf3);
            fixture.teardown;
            
            testCase.verifyEqual(chhd1.integralNumerical,1,'AbsTol',1e-5);
            testCase.verifyEqual(chhd2.integralNumerical,1,'AbsTol',1e-5);
            testCase.verifyEqual(chhd3.integralNumerical,1,'AbsTol',1e-5);
        end
        
    end
end