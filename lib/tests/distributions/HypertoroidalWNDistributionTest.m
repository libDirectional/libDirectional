classdef HypertoroidalWNDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testIntegral(testCase)
            % Test integral in 3d case
            rng default
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu = 2*pi*rand(3,1);
            hwnd = HypertoroidalWNDistribution(mu,C);
            testCase.verifyEqual(integral3(@(x,y,z)reshape(hwnd.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,0,2*pi,0,2*pi),1,'AbsTol',1E-3);
        end
        
        function test2D(testCase)
            % Test against ToroidalWNDistribution in 2d case
            rng default
            C = [0.7,0.4;0.4,0.6];
            mu = 2*pi*rand(2,1);
            twnd = ToroidalWNDistribution(mu,C);
            hwnd = HypertoroidalWNDistribution(mu,C);
            testPoints = rand(2,200);
            testCase.verifyEqual(twnd.pdf(testPoints),hwnd.pdf(testPoints),'AbsTol',1E-5);
        end
        
        function test1D(testCase)
            % Test against WNDistribution in 1d case
            rng default
            sigma = 1.4;
            mu = 0.9;
            wnd = WNDistribution(mu,sigma);
            hwnd = HypertoroidalWNDistribution(mu,sigma^2);
            testPoints = rand(1,200);
            testCase.verifyEqual(wnd.pdf(testPoints),hwnd.pdf(testPoints),'AbsTol',1E-5);
        end
    end
end