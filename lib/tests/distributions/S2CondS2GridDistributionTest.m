classdef S2CondS2GridDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testWarningFreeNormalizedVMF(testCase)
            noGridPoints = 112;
            trans = @(xkk,xk)cell2mat(arrayfun(@(i)...
                VMFDistribution(xk(:,i),1).pdf(xkk),1:size(xk,2),'UniformOutput',false));
            testCase.verifyWarningFree(@()...
                S2CondS2GridDistribution.fromFunction(trans,noGridPoints,true,'eq_point_set'));
        end
        function testErrorUnnormalized(testCase)
            noGridPoints = 112;
            trans = @(xkk,xk)1./(sum((diag([0.1,0.15,1])*(xkk-xk)).^2,1)+0.01);
            testCase.verifyWarning(@()...
                S2CondS2GridDistribution.fromFunction(trans,noGridPoints,false,'eq_point_set'),...
                'Normalization:unnormalized');
        end
        function testWarningFreeCustomNormalized(testCase)
            noGridpoints = 1000;
            testCase.verifyWarningFree(@()...
                S2CondS2GridDistribution.fromFunction(@trans,noGridpoints,true,'eq_point_set'));
            function p = trans(xkk,xk)
                transUnnorm = @(xkk,xk)1./(sum((diag([0.1,0.15,0.3])*(xkk-xk)).^2,1)+0.01);
                p = NaN(size(xkk,2),size(xk,2));
                for i=1:size(xk,2)
                    chd = CustomHypersphericalDistribution(@(xkk)transUnnorm(xkk,xk(:,i)),3);
                    p(:,i) = chd.pdf(xkk)/chd.integralNumerical;
                end
            end
        end
        function testEqualWithAndWithoutCart(testCase)
            noGridpoints = 2001;
            dist = VMFDistribution([0;-1;0], 100);
            
            fTrans1 = @(xkk,xk)repmat(dist.pdf(xkk)',[1,size(xk,2)]);
            s2s2_fTrans1 = S2CondS2GridDistribution.fromFunction(fTrans1, noGridpoints, true);
            
            fTrans2 = @(xkk,xk)dist.pdf(xkk);
            s2s2_fTrans2 = S2CondS2GridDistribution.fromFunction(fTrans2, noGridpoints, false);
            
            testCase.verifyEqual(s2s2_fTrans1.getGrid(),s2s2_fTrans2.getGrid());
            testCase.verifyEqual(s2s2_fTrans1.gridValues,s2s2_fTrans2.gridValues);
        end
        function testPlottingWarningFree(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            noGridPoints = 112;
            trans=@(xkk,xk)cell2mat(arrayfun(@(i)...
                VMFDistribution(xk(:,i),1).pdf(xkk),1:size(xk,2),'UniformOutput',false));
            s2s2=S2CondS2GridDistribution.fromFunction(trans,noGridPoints,true,'eq_point_set');
            
            testCase.applyFixture(SuppressedWarningsFixture('MATLAB:hg:AutoSoftwareOpenGL'));
            testCase.verifyWarningFree(@()s2s2.plot);
            testCase.verifyWarningFree(@()s2s2.plotInterpolated);
        end
        function testFixDim(testCase)
            noGridPoints = 112;
            trans1 = @(xkk,xk)cell2mat(arrayfun(@(i)...
                VMFDistribution(xk(:,i),1).pdf(xkk),1:size(xk,2),'UniformOutput',false)); 
            s2s2 = S2CondS2GridDistribution.fromFunction(trans1,noGridPoints,true);
            
            for point = s2s2.getGrid()
                sgd = s2s2.fixDim(2,point);
                testCase.verifyEqual(sgd.meanDirection,point,'AbsTol',1E-2);
            end
        end
    end
        
end
