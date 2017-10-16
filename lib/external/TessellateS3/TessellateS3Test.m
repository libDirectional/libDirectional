classdef TessellateS3Test < matlab.unittest.TestCase
    % Unit test for tessellate_S3
    % start with
    %   >> run(TessellateS3Test)
   
    properties
    end
    
    methods (Test)
        function testTessellateS3(testCase)
            for n=0:4
                m = 16*8^n;
                x = tessellate_S3(m-1)';
                testCase.verifySize(x, [m 4]);
                testCase.verifyEqual(sum(x.^2,2), ones(m,1), 'RelTol', 1E-10); %check unit norm
                testCase.verifyEqual(sum(x,1), zeros(1,4), 'AbsTol', 1E-10); %check mean zero
            end
        end        
    end
end
