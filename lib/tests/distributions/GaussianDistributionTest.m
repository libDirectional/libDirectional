classdef GaussianDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testGaussianDistribution1d(testCase)
            mu = 3;
            sigma = 4;
            C = sigma^2;
            g = GaussianDistribution(mu,C);
            
            %% test pdf for 1D case
            testCase.verifyEqual(g.pdf(mu),1/sqrt(2*pi)/sigma * exp(0), 'RelTol', 1E-10);
            testCase.verifyEqual(g.pdf(mu+2),1/sqrt(2*pi)/sigma * exp(-4/2/sigma^2), 'RelTol', 1E-10);
            xa = -5:5;
            testCase.verifyEqual(g.pdf(xa),1/sqrt(2*pi)/sigma * exp(-(xa-mu).^2/2/sigma^2), 'RelTol', 1E-10);
            testCase.verifyEqual(g.pdf(xa), mvnpdf(xa',mu,C)', 'RelTol', 1E-10);
            
            %% test normalization
            gInt = integral(@(x) g.pdf(x), -Inf, Inf);
            testCase.verifyEqual(gInt, 1, 'RelTol', 1E-10);
            
            %% test conversion to WN
            wn = g.toWN();
            testCase.verifyClass(wn, 'WNDistribution');
            testCase.verifyEqual(g.mu, wn.mu);
            testCase.verifyEqual(g.C, wn.sigma^2); % adjust this check, if a more sophisticated conversion is used
            
            %% test sampling
            n = 10;
            s = g.sample(n);
            testCase.verifyEqual(size(s,1), 1);
            testCase.verifyEqual(size(s,2), n);
        end
        
        function testGaussianDistribution2d(testCase)
            mu = [2,3]';
            C = [ 1.1 0.4; 0.4 0.9];
            g = GaussianDistribution(mu,C);
            
            %% test pdf for 2D case
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
            
            %% test normalization
            gInt = integral2(@(x,y) reshape(g.pdf([x(:)';y(:)']), size(x,1), size(x,2)), -20, 20, -20, 20, 'RelTol', 1E-5);
            testCase.verifyEqual(gInt, 1, 'RelTol', 1E-5);
            
            %% test sampling
            n = 10;
            s = g.sample(10);
            testCase.verifyEqual(size(s,1), 2);
            testCase.verifyEqual(size(s,2), n);
        end
    end
end