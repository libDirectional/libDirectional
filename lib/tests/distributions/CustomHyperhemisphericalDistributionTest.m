classdef CustomHyperhemisphericalDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testSimpleDistribution2D(testCase)
            M = eye(3);
            Z = [-2 -0.5 0]';
            bd = BinghamDistribution(Z,M);
            chhd = CustomHyperhemisphericalDistribution.fromDistribution(bd);
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
            chhd = CustomHyperhemisphericalDistribution.fromDistribution(bd);

            testCase.verifyEqual(chhd.integralNumerical,1,'AbsTol',1e-4);
        end
        
        function testWarningAsymmetric(testCase)
            vmf = VMFDistribution([0;0;1],10);
            testCase.verifyWarning(...
                @()CustomHyperhemisphericalDistribution.fromDistribution(vmf),...
                'FromDistribution:UsePdfHypersphere');
        end
        function testIntegralVmfS2(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vmf1 = VMFDistribution([0;0;1],10);
            vmf2 = VMFDistribution(1/sqrt(2)*[0;1;1],10);
            vmf3 = VMFDistribution([0;1;0],10); % This is antipodally symmetric
            fixture = testCase.applyFixture(SuppressedWarningsFixture('FromDistribution:UsePdfHypersphere'));
            chhd1 = CustomHyperhemisphericalDistribution.fromDistribution(vmf1);
            chhd2 = CustomHyperhemisphericalDistribution.fromDistribution(vmf2);
            chhd3 = CustomHyperhemisphericalDistribution.fromDistribution(vmf3);
            fixture.teardown;
            
            testCase.verifyEqual(chhd1.integralNumerical,1,'AbsTol',1e-5);
            testCase.verifyEqual(chhd2.integralNumerical,1,'AbsTol',1e-5);
            testCase.verifyEqual(chhd3.integralNumerical,1,'AbsTol',1e-5);
        end
        
    end
end