classdef HypertoroidalWNDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testIntegral(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true; end
            
            rng default
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu = 2*pi*rand(3,1);
            hwnd = HypertoroidalWNDistribution(mu,C);
            
            % Test integral in 3d case
            if enableExpensive
                testCase.verifyEqual(integral3(@(x,y,z)reshape(hwnd.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,0,2*pi,0,2*pi),1,'AbsTol',1E-3);
            end
            
            function p = pdf (x, mu, C)
                p=0;
                for i=-3:3
                    for j=-3:3
                        for k=-3:3
                            p = p + mvnpdf(x,mu + 2*pi*[i;j;k], C);
                        end
                    end
                end
            end
            
            testPoints = rand(3,20);
            pdfVals = hwnd.pdf(testPoints);
            for t=1:size(testPoints,2)
                testCase.verifyEqual(pdfVals(t), pdf(testPoints(:,t), hwnd.mu, hwnd.C), 'RelTo', 1E-10);
            end
            
            testCase.verifyEqual(hwnd.toGaussian().mu, hwnd.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(hwnd.toGaussian().C, hwnd.C, 'RelTol', 1E-10); 
        end
        
        function test2D(testCase)
            % Test against ToroidalWNDistribution in 2d case
            rng default
            C = [0.7,0.4;0.4,0.6];
            mu = 2*pi*rand(2,1);
            twnd1 = ToroidalWNDistribution(mu,C);
            hwnd = HypertoroidalWNDistribution(mu,C);
            twnd2 = hwnd.toToroidalWN();
            testCase.verifyClass(twnd2, 'ToroidalWNDistribution');
            testPoints = rand(2,100);
            for twnd = [twnd1 twnd2]
                testCase.verifyEqual(twnd.pdf(testPoints),hwnd.pdf(testPoints),'AbsTol',1E-5);
                testCase.verifyEqual(twnd.trigonometricMoment(0),hwnd.trigonometricMoment(0),'RelTol',1E-10);
                testCase.verifyEqual(twnd.trigonometricMoment(1),hwnd.trigonometricMoment(1),'RelTol',1E-10);
                testCase.verifyEqual(twnd.trigonometricMoment(2),hwnd.trigonometricMoment(2),'RelTol',1E-10);
                testCase.verifyEqual(twnd.trigonometricMoment(3),hwnd.trigonometricMoment(3),'RelTol',1E-10);
                testCase.verifyEqual(twnd.toGaussian().mu, hwnd.toGaussian().mu, 'RelTol', 1E-10);
                testCase.verifyEqual(twnd.toGaussian().C, hwnd.toGaussian().C, 'RelTol', 1E-10);
            end
        end
        
        function test1D(testCase)
            % Test against WNDistribution in 1d case
            rng default
            sigma = 1.4;
            mu = 0.9;
            wnd1 = WNDistribution(mu,sigma);
            hwnd = HypertoroidalWNDistribution(mu,sigma^2);
            wnd2 = hwnd.toWN();
            testCase.verifyClass(wnd2, 'WNDistribution');
            testPoints = rand(1,100);
            for wnd = [wnd1 wnd2]
                testCase.verifyEqual(wnd.pdf(testPoints),hwnd.pdf(testPoints),'AbsTol',1E-5);
                testCase.verifyEqual(wnd.trigonometricMoment(0),hwnd.trigonometricMoment(0),'RelTol',1E-10);
                testCase.verifyEqual(wnd.trigonometricMoment(1),hwnd.trigonometricMoment(1),'RelTol',1E-10);
                testCase.verifyEqual(wnd.trigonometricMoment(2),hwnd.trigonometricMoment(2),'RelTol',1E-10);
                testCase.verifyEqual(wnd.trigonometricMoment(3),hwnd.trigonometricMoment(3),'RelTol',1E-10);
                testCase.verifyEqual(wnd.toGaussian().mu, hwnd.toGaussian().mu, 'RelTol', 1E-10);
                testCase.verifyEqual(wnd.toGaussian().C, hwnd.toGaussian().C, 'RelTol', 1E-10);
            end
        end
        
        function testSampling(testCase)
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu = 2*pi*rand(3,1);
            hwnd = HypertoroidalWNDistribution(mu,C);
            
            nSamples = 5;
            s = hwnd.sample(nSamples);
            testCase.verifyEqual(size(s,1),3);
            testCase.verifyEqual(size(s,2),nSamples);
            testCase.verifyEqual(s, mod(s,2*pi));
        end
        
        function testShift(testCase)
            twn = HypertoroidalWNDistribution([3;5;6],eye(3));
            s = [4;2;-2];
            twn2 = twn.shift(s);
            testCase.verifyClass(twn2, 'HypertoroidalWNDistribution');           
            [xTest,yTest]=meshgrid(linspace(0,2*pi,20));
            zTest = xTest;
            testCase.verifyEqual(twn2.pdf([xTest(:)';yTest(:)';zTest(:)']),twn.pdf([xTest(:)' - s(1);yTest(:)' - s(2); zTest(:)' - s(3)]),'AbsTol',1E-10);
        end
    end
end