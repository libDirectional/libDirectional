classdef HyperhemisphericalUniformDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testSimpleDistribution2D(testCase)
            hhud = HyperhemisphericalUniformDistribution(3);
            % Use nonrandom numbers but do not influece seed everywhere
            % else
            rs = RandStream('twister','Seed',10);
            points = rs.randn(3,100);
            points(:,points(3,:)<0)=[];
            points = points./vecnorm(points);
            
            testCase.verifyEqual(hhud.pdf(points), ones(1,size(points,2))*1/(2*pi),'AbsTol',1e-6);
        end
        function testIntegralS2(testCase)
            hhud = HyperhemisphericalUniformDistribution(3);
            testCase.verifyEqual(hhud.integral, 1,'AbsTol',1e-6);
        end
        function testIntegralS3(testCase)
            hhud = HyperhemisphericalUniformDistribution(4);
            testCase.verifyEqual(hhud.integral, 1,'AbsTol',1e-6);
        end
    end
end