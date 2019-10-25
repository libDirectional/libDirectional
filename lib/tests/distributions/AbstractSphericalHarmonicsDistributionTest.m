classdef AbstractSphericalHarmonicsDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testConvolutionDimAndMeanAndNormalization(testCase)
            % Test that the function is working for differing coefficient
            % matrix sizes
            coeffMat = rand(7, 13);
            coeffMat(1) = 1 / sqrt(4*pi);
            rshdRandom = SphericalHarmonicsDistributionReal(coeffMat);
            coeffMatZonal = zeros(7, 13);
            coeffMatZonal(1) = 1 / sqrt(4*pi);
            coeffMatZonal(sub2ind(size(coeffMatZonal), 2:size(coeffMatZonal, 1), 2:size(coeffMatZonal, 1))) = rand(1, size(coeffMatZonal, 1)-1);
            cshdZonal = SphericalHarmonicsDistributionComplex(coeffMatZonal);
            
            % Test both for real and complex
            shdRandoms = {rshdRandom, rshdRandom.toSphericalHarmonicsDistributionComplex};
            shdZonals = {cshdZonal.toSphericalHarmonicsDistributionReal, cshdZonal};
            for i = 1:2
                shdRandomCurr = shdRandoms{i};
                shdZonalCurr = shdZonals{i};
                shdtmp = shdRandomCurr.convolve(shdZonalCurr);
                testCase.verifySize(shdtmp.coeffMat, [7, 13]);
                testCase.verifyEqual(shdtmp.meanDirection, shdRandomCurr.meanDirection, 'AbsTol', 1E-10);
                testCase.verifyEqual(shdtmp.meanDirection, shdRandomCurr.meanDirection, 'AbsTol', 1E-4);
                
                shdtmp = shdRandomCurr.convolve(shdZonalCurr.truncate(5));
                testCase.verifySize(shdtmp.coeffMat, [6, 11]);
                testCase.verifyEqual(shdtmp.meanDirection, shdRandomCurr.meanDirection, 'AbsTol', 1E-10);
                testCase.verifyEqual(shdtmp.integralNumerical, 1, 'AbsTol', 1E-4);
                
                shdtmp = shdRandomCurr.convolve(shdZonalCurr.truncate(4));
                testCase.verifySize(shdtmp.coeffMat, [5, 9]);
                testCase.verifyEqual(shdtmp.meanDirection, shdRandomCurr.meanDirection, 'AbsTol', 1E-10);
                testCase.verifyEqual(shdtmp.integralNumerical, 1, 'AbsTol', 1E-4);
                
                shdRandomTrunc = shdRandomCurr.truncate(5);
                shdtmp = shdRandomTrunc.convolve(shdZonalCurr);
                testCase.verifySize(shdtmp.coeffMat, [6, 11]);
                testCase.verifyEqual(shdtmp.meanDirection, shdRandomCurr.meanDirection, 'AbsTol', 1E-10);
                testCase.verifyEqual(shdtmp.integralNumerical, 1, 'AbsTol', 1E-4);
                
                shdRandomTrunc = shdRandomCurr.truncate(4);
                shdtmp = shdRandomTrunc.convolve(shdZonalCurr);
                testCase.verifySize(shdtmp.coeffMat, [5, 9]);
                testCase.verifyEqual(shdtmp.meanDirection, shdRandomCurr.meanDirection, 'AbsTol', 1E-10);
                testCase.verifyEqual(shdtmp.integralNumerical, 1, 'AbsTol', 1E-4);
            end
        end
        function testConvolution(testCase)
            function p = vmfpdf(xa, mu, kappa)
                dim = 3;
                C = kappa^(dim / 2 - 1) / ((2 * pi)^(dim / 2) * besseli(dim/2-1, kappa));
                p = C * exp(kappa*mu'*xa);
            end
            order = 50;
            kappaTran = 1;
            kappaBeforePredict = 1;
            
            fBeforePredict = VMFDistribution([0; 1; 0], kappaBeforePredict);
            shdBeforePredict = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(fBeforePredict, order);
            
            fTran = @(mu, xkk)vmfpdf(xkk, mu, kappaTran);
            r = 1;
            fangles = @(phi1, phi2, xkk) fTran([ ... % xk (mu) is influenced by varying phi1 and phi2
                r .* sin(phi1) .* sin(phi2); ...
                r .* cos(phi1) .* sin(phi2); ...
                r .* cos(phi2)], xkk);
            
            g = @(phi1, phi2, xkk) reshape(fangles(phi1(:)', phi2(:)', xkk).*sin(phi2(:)')', size(phi1));
            
            fBeforePredictAng = @(phi1, phi2) fBeforePredict.pdf([ ...
                r .* sin(phi1) .* sin(phi2); ...
                r .* cos(phi1) .* sin(phi2); ...
                r .* cos(phi2)]);
            
            fBeforePredictAngArbitraryInput = @(phi1, phi2) reshape(fBeforePredictAng(phi1(:)', phi2(:)'), size(phi1));
            
            fpredscalar = @(xkk)integral2(@(phi1, phi2)g(phi1, phi2, xkk).*fBeforePredictAngArbitraryInput(phi1, phi2) ...
                , 0, 2*pi, 0, pi, 'AbsTol', 1e-5, 'RelTol', 1e-5);
            fpredArrayfun = @(xkk)arrayfun(@(x, y, z)fpredscalar([x; y; z]), xkk(1, :), xkk(2, :), xkk(3, :));
            
            chd = CustomHypersphericalDistribution(fpredArrayfun, 3);
            
            shdTrans = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(VMFDistribution([0; 0; 1], kappaTran), order);
            shdRes = shdBeforePredict.convolve(shdTrans);
            
            testCase.verifyEqual(shdRes.hellingerDistanceNumerical(chd), 0, 'AbsTol', 1E-10);
        end
    end
end
