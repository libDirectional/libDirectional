classdef MvnpdffastTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testMvnpdffast(testCase)
            for n= [1 2 3 4 5 10]
                % n-D
                rng default
                mu = rand(1,n);
                C = rand(n,n);
                C=C*C';
                x = rand(100,n);

                testCase.verifyEqual(mvnpdf(x, mu, C), mvnpdffast(x, mu, C),'RelTol', 1E-9);
            end
        end
    end
end
