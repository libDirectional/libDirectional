classdef GaussianDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testGaussianDistribution1d(testCase)
            mu = 3;
            sigma = 4;
            C = sigma^2;
            g = GaussianDistribution(mu,C);
            
            % test pdf for 1D case
            testCase.verifyEqual(g.pdf(mu),1/sqrt(2*pi)/sigma * exp(0), 'RelTol', 1E-10);
            testCase.verifyEqual(g.pdf(mu+2),1/sqrt(2*pi)/sigma * exp(-4/2/sigma^2), 'RelTol', 1E-10);
            xa = -5:5;
            testCase.verifyEqual(g.pdf(xa),1/sqrt(2*pi)/sigma * exp(-(xa-mu).^2/2/sigma^2), 'RelTol', 1E-10);
            testCase.verifyEqual(g.pdf(xa), mvnpdf(xa',mu,C)', 'RelTol', 1E-10);
            
            % test normalization
            gInt = integral(@(x) g.pdf(x), -Inf, Inf);
            testCase.verifyEqual(gInt, 1, 'RelTol', 1E-10);
            
            % test conversion to WN
            wn = g.toWN();
            testCase.verifyClass(wn, 'WNDistribution');
            testCase.verifyEqual(g.mu, wn.mu);
            testCase.verifyEqual(g.C, wn.sigma^2); % adjust this check, if a more sophisticated conversion is used
            
            % test sampling
            n = 10;
            s = g.sample(n);
            testCase.verifySize(s, [1 n]);        
            
            [s,w] = g.sampleDeterministic();
            n=3;
            testCase.verifySize(s, [1 n]); 
            testCase.verifySize(w, [1 n]); 
            testCase.verifyEqual(w, ones(1,n)/n);
            testCase.verifyEqual(mean(s), g.mu);
            testCase.verifyEqual(cov(s,1), g.C, 'RelTol', 1E-10);
        end
        
        function testGaussianDistribution2d(testCase)
            mu = [2,3]';
            C = [ 1.1 0.4; 0.4 0.9];
            g = GaussianDistribution(mu,C);
            
            % test pdf for 2D case
            testCase.verifyEqual(g.pdf(mu),1/2/pi/sqrt(det(C)) * exp(0), 'RelTol', 1E-10);
            x = [0.2,0.3]';
            testCase.verifyEqual(g.pdf(x),1/2/pi/sqrt(det(C))* exp(-0.5 * (x-mu)' /(C)*(x-mu) ), 'RelTol', 1E-10);
            x = [3,2.9]';
            testCase.verifyEqual(g.pdf(x),1/2/pi/sqrt(det(C))* exp(-0.5 * (x-mu)' /(C)*(x-mu) ), 'RelTol', 1E-10);
            xa = [-5:5; -3:7];
            ya = zeros(1,size(xa,2));
            for i=1:length(xa)
                ya(i) = 1/2/pi/sqrt(det(C))* exp(-0.5 * (xa(:,i)-mu)' /(C)*(xa(:,i)-mu) );
            end
            testCase.verifyEqual(g.pdf(xa), ya, 'RelTol', 1E-10);
            testCase.verifyEqual(g.pdf(xa), mvnpdf(xa',mu',C)', 'RelTol', 1E-10);
            
            % test normalization
            gInt = integral2(@(x,y) reshape(g.pdf([x(:)';y(:)']), size(x,1), size(x,2)), -20, 20, -20, 20, 'RelTol', 1E-5);
            testCase.verifyEqual(gInt, 1, 'RelTol', 1E-5);
            
            % test sampling
            n = 10;
            s = g.sample(10);
            testCase.verifySize(s, [2 n]);      
            
            [s,w] = g.sampleDeterministic();
            n=5;
            testCase.verifySize(s, [2 n]); 
            testCase.verifySize(w, [1 n]); 
            testCase.verifyEqual(w, ones(1,n)/n);
            testCase.verifyEqual(mean(s,2), g.mu);
            testCase.verifyEqual(cov(s',1), g.C, 'RelTol', 1E-10);            
        end
        
        function testGaussianDistribution3d(testCase)
            mu = [2,3,4]';
            C = [ 1.1 0.4 0; 0.4 0.9 0; 0 0 0.1];
            g = GaussianDistribution(mu,C);
            
            % test pdf for 3D case
            xa = [-5:5; -3:7; -2:8];
            testCase.verifyEqual(g.pdf(xa), mvnpdf(xa',mu',C)', 'RelTol', 1E-10);
                        
            % test sampling
            n = 10;
            s = g.sample(10);
            testCase.verifySize(s, [3 n]);
            
            [s,w] = g.sampleDeterministic();
            n=7;
            testCase.verifySize(s, [3 n]); 
            testCase.verifySize(w, [1 n]); 
            testCase.verifyEqual(w, ones(1,n)/n);
            testCase.verifyEqual(mean(s,2), g.mu);
            testCase.verifyEqual(cov(s',1), g.C, 'RelTol', 1E-10);            
        end     
        
        function testSampleDeterministic(testCase)
            mu = [2,3,4]';
            C = [ 1.1 0.4 0; 0.4 0.9 0; 0 0 0.1];
            g = GaussianDistribution(mu,C);
            
            [s,w] = g.sampleDeterministic();
            nSamples = 2*3+1;
            testCase.verifySize(s, [3 nSamples]);        
            testCase.verifySize(w, [1 nSamples]);        
            testCase.verifyEqual(sum(s.*repmat(w,3,1),2), mu, 'RelTol', 1E-10);
            testCase.verifyEqual((s-repmat(mu,1,nSamples))*diag(w)*(s-repmat(mu,1,nSamples))', C, 'RelTol', 1E-10);
        end
        
        function testMode(testCase)
            mu = [1,2,3]';
            C = [ 1.1 0.4 0; 0.4 0.9 0; 0 0 1];
            g = GaussianDistribution(mu,C);
            
            testCase.verifyEqual(g.mode, mu, 'AbsTol',1e-6);
        end
        
        function testShift(testCase)
            mu = [1,2,3]';
            C = [ 1.1 0.4 0; 0.4 0.9 0; 0 0 1];
            g = GaussianDistribution(mu,C);
            
            shiftBy = [2;-2;3];
            gShifted = g.shift(shiftBy);
            
            testCase.verifyEqual(gShifted.mode, mu+shiftBy, 'AbsTol',1e-6);
        end

        function testMarginalization(testCase)
            mu = [1,2,3]';
            C = [ 1.1 0.4 0; 0.4 0.9 0; 0 0 1];
            g = GaussianDistribution(mu,C);
            
            grid = linspace(-10,10,300);
            distMarginalized = g.marginalizeOut([2,3]);
            marginlized1DViaIntegral=@(x)arrayfun(@(xCurr)integral2(@(y,z)reshape(g.pdf([xCurr*ones(1,size(y,2));y(:)';z(:)']),size(y)),-inf,inf,-inf,inf),x);
            testCase.verifyEqual(distMarginalized.pdf(grid), marginlized1DViaIntegral(grid), 'AbsTol', 1E-10);
        end
    end
end