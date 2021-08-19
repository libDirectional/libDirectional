classdef LinearMixtureTest < matlab.unittest.TestCase
    methods(Test)
        function testConstructor(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            makeMixture = @()LinearMixture({GaussianDistribution(1,1),GaussianDistribution(50,1)},[0.3,0.7]);
            testCase.verifyWarning(makeMixture,'LinearMixture:AllGaussians');
            testCase.applyFixture(SuppressedWarningsFixture('LinearMixture:AllGaussians'));
            testCase.verifyError(@()LinearMixture({GaussianDistribution(1,1)},[0.3,0.7]),'Mixture:IncompatibleDimensions');
        end
        function testPdf(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            gm1 = GaussianDistribution([1;1],diag([2,3]));
            gm2 = GaussianDistribution(-[3;1],diag([2,3]));
            makeMixture = @()LinearMixture({gm1,gm2},[0.3,0.7]);
            testCase.applyFixture(SuppressedWarningsFixture('LinearMixture:AllGaussians'));
            lm = makeMixture();
            [x,y] = meshgrid(linspace(-2,2,100));
            
            testCase.verifyEqual(lm.pdf([x(:)';y(:)']),0.3*gm1.pdf([x(:)';y(:)'])+0.7*gm2.pdf([x(:)';y(:)']),'AbsTol',1E-20);
        end
    end
end
