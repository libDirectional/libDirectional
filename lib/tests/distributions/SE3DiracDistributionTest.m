classdef SE3DiracDistributionTest < matlab.unittest.TestCase    
    methods (Test)
        function testConstructor(testCase)
            dSph = [ 1 2 3 4 5 6;
                  2 4 0 0.5 1 1;
                  5 10 20 30 40 50];
            dSph = dSph./vecnorm(dSph);
            dLin = [-5,0,5,10,15,20];
            w = [ 1 2 3 1 2 3];
            w = w/sum(w);
            testCase.verifyWarningFree(@()HypercylindricalDiracDistribution(4, [dSph;dLin], w));
        end
        
        function testFromDistribution(testCase)
            cpsd = SE3CartProdStackedDistribution({HyperhemisphericalUniformDistribution(4),GaussianDistribution([1;2;3],diag([3,2,1]))});
            testCase.verifyWarningFree(@()SE3DiracDistribution.fromDistribution(cpsd,100));
        end
    end
end