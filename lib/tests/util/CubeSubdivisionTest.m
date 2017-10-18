classdef CubeSubdivisionTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testTesseractSubdivision(testCase)
            for n=0:3
                m = 6*4^n+2;
                x = cubesubdivision(m-1);
                testCase.verifySize(x, [m 3]);
                testCase.verifyEqual(sum(x.^2,2), ones(m,1), 'RelTol', 1E-10); %check unit norm
                testCase.verifyEqual(sum(x,1), zeros(1,3), 'AbsTol', 1E-10); %check mean zero
                testCase.verifyEqual(x+flipud(x), zeros(m,3), 'AbsTol', 1E-10); %check symmetry
            end
        end
        
        function testCubeSubdivisionUnnormalized(testCase)
            for n=0:3
                m = 6*4^n+2;
                x = cubesubdivision(m-1, false);
                testCase.verifySize(x, [m 3]);
                testCase.verifyLessThanOrEqual(abs(x), ones(m,3)); %inside cube
                testCase.verifyEqual(max(abs(x),[],2), ones(m,1), 'AbsTol', 1E-10); %check on cube surface
                testCase.verifyEqual(sum(x,1), zeros(1,3), 'AbsTol', 1E-10); %check mean zero
                testCase.verifyEqual(x+flipud(x), zeros(m,3), 'AbsTol', 1E-10); %check symmetry
            end
        end          
    end
end
