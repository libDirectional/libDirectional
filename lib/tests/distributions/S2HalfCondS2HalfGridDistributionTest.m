classdef S2HalfCondS2HalfGridDistributionTest < matlab.unittest.TestCase
    methods(Test)
        % Test cases adapted from S2HalfCondS2HalfGridDistributionTest
        function testClassWatsonS2(testCase)
            noGridPoints = 112;
            trans = @(xkk,xk)2*cell2mat(arrayfun(@(i)...
                WatsonDistribution(xk(:,i),1).pdf(xkk),1:size(xk,2),'UniformOutput',false));
            sd = S2HalfCondS2HalfGridDistribution.fromFunction(trans,noGridPoints,true,'eq_point_set');
            testCase.verifyClass(sd, 'S2HalfCondS2HalfGridDistribution');
        end
        function testWarningFreeNormalizedWatsonS2(testCase)
            noGridPoints = 112;
            trans = @(xkk,xk)2*cell2mat(arrayfun(@(i)...
                WatsonDistribution(xk(:,i),1).pdf(xkk),1:size(xk,2),'UniformOutput',false));
            testCase.verifyWarningFree(@()...
                S2HalfCondS2HalfGridDistribution.fromFunction(trans,noGridPoints,true,'eq_point_set'));
        end
        function testWarningUnnormalizedS2(testCase)
            noGridPoints = 112;
            trans = @(xkk,xk)1./(sum((diag([0.1,0.15,1])*(xkk-xk)).^2,1)+0.01);
            testCase.verifyWarning(@()...
                S2HalfCondS2HalfGridDistribution.fromFunction(trans,noGridPoints,false,'eq_point_set'),...
                'Normalization:unnormalized');
        end
        function testWarningFreeCustomNormalized(testCase)
            noGridpoints = 1000;
            testCase.verifyWarningFree(@()...
                S2HalfCondS2HalfGridDistribution.fromFunction(@trans,noGridpoints,true,'eq_point_set'));
            function p = trans(xkk,xk)
                transUnnorm = @(xkk,xk)2./((sum((diag([0.1,0.15,0.3])*(xkk-xk)).^2,1)+sum((diag([0.1,0.15,0.3])*(xkk+xk)).^2,1))/2+0.01);
                p = NaN(size(xkk,2),size(xk,2));
                for i=1:size(xk,2)
                    chd = CustomHypersphericalDistribution(@(xkk)transUnnorm(xkk,xk(:,i)),3);
                    p(:,i) = chd.pdf(xkk)/chd.integralNumerical;
                end
                p = 2*p;
            end
        end
        function testEqualWithAndWithoutCartS2(testCase)
            noGridpoints = 2001;
            dist = HypersphericalMixture({VMFDistribution([0;-1;0], 100),VMFDistribution([0;1;0], 100)},[0.5,0.5]);

            fTrans1 = @(xkk,xk)repmat(2*dist.pdf(xkk)',[1,size(xk,2)]);
            sdsd_fTrans1 = S2HalfCondS2HalfGridDistribution.fromFunction(fTrans1, noGridpoints, true, 'eq_point_set');
            
            fTrans2 = @(xkk,xk)2*dist.pdf(xkk);
            sdsd_fTrans2 = S2HalfCondS2HalfGridDistribution.fromFunction(fTrans2, noGridpoints, false, 'eq_point_set');
            
            testCase.verifyEqual(sdsd_fTrans1.getGrid(),sdsd_fTrans2.getGrid());
            testCase.verifyEqual(sdsd_fTrans1.gridValues,sdsd_fTrans2.gridValues);
        end
        function testIndepdenentS2(testCase)
            noGridpoints = 100;
            % Prediction with Bingham
            Z = [-7.44;-2.84;0];
            M = [-0.18 -0.172 -0.969;-0.738 -0.627 0.248;-0.65 0.76 -0.0141];
            distPredictionResult = BinghamDistribution(Z,M);
            distPredictionResult.F = distPredictionResult.F*distPredictionResult.integralNumerical;
            
            fTrans = @(xkk,xk)2*repmat(distPredictionResult.pdf(xkk)',[1,size(xk,2)]);
            sdsd_fTrans = S2HalfCondS2HalfGridDistribution.fromFunction(fTrans, noGridpoints, true, 'eq_point_set');
            testCase.verifyEqual(diff(sdsd_fTrans.gridValues,1,2),zeros(100,99));
        end
        function testPlottingWarningFree(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            noGridPoints = 112;
            trans = @(xkk,xk)2*cell2mat(arrayfun(@(i)...
                WatsonDistribution(xk(:,i),1).pdf(xkk),1:size(xk,2),'UniformOutput',false));
            sdsd = S2HalfCondS2HalfGridDistribution.fromFunction(trans,noGridPoints,true,'eq_point_set');
            
            testCase.applyFixture(SuppressedWarningsFixture('MATLAB:hg:AutoSoftwareOpenGL'));
            testCase.verifyWarningFree(@()sdsd.plotInterpolated);
            testCase.verifyWarningFree(@()sdsd.plotFullSphere);
            testCase.verifyWarningFree(@()sdsd.plotInterpolatedFullSphere);
        end
        function testFixDimS2(testCase)
            noGridPoints = 112;
            trans = @(xkk,xk)2*cell2mat(arrayfun(@(i)...
                WatsonDistribution(xk(:,i),1).pdf(xkk),1:size(xk,2),'UniformOutput',false));
            sdsd = S2HalfCondS2HalfGridDistribution.fromFunction(trans,noGridPoints,true,'eq_point_set');
            
            grid = sdsd.getGrid();
            for i = 1:size(grid,2)
                point = grid(:,i);
                sgd = sdsd.fixDim(2,point);
                testCase.verifyEqual(sgd.gridValues, sdsd.gridValues(:,i));
                
                sgd = sdsd.fixDim(1,point);
                testCase.verifyEqual(sgd.gridValues', sdsd.gridValues(i,:));
            end
        end
    end
        
end
