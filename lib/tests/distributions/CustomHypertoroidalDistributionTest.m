classdef CustomHypertoroidalDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testSimpleDistribution2D(testCase)
            twn = ToroidalWNDistribution([1;2],eye(2));
            chd = CustomHypertoroidalDistribution(@(xa)twn.pdf(xa),2);
            ct1 = CustomToroidalDistribution(@(xa)twn.pdf(xa));
            ct2 = chd.toCustomToroidal();
            [xTest,yTest]=meshgrid(linspace(0,2*pi,50));
            testCase.verifyEqual(chd.pdf([xTest(:)';yTest(:)']),twn.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-10);
            testCase.verifyEqual(ct1.pdf([xTest(:)';yTest(:)']),twn.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-10);
            testCase.verifyEqual(ct2.pdf([xTest(:)';yTest(:)']),twn.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-10);
        end
        
        function testSimpleDistribution1D(testCase)
            wn = WNDistribution(2, 0.7);
            chd = CustomHypertoroidalDistribution(@(xa)wn.pdf(xa),1);
            ct1 = CustomCircularDistribution(@(xa)wn.pdf(xa));
            ct2 = chd.toCustomCircular();
            xTest = linspace(0,2*pi,50);
            testCase.verifyEqual(chd.pdf(xTest),wn.pdf(xTest),'AbsTol',1E-10);
            testCase.verifyEqual(ct1.pdf(xTest),wn.pdf(xTest),'AbsTol',1E-10);
            testCase.verifyEqual(ct2.pdf(xTest),wn.pdf(xTest),'AbsTol',1E-10);
        end        
        
        function testFromDistribution(testCase)
            twn=ToroidalWNDistribution([1;2],eye(2));
            chd=CustomHypertoroidalDistribution.fromDistribution(twn);
            [xTest,yTest]=meshgrid(linspace(0,2*pi,100));
            testCase.verifyEqual(chd.pdf([xTest(:)';yTest(:)']),twn.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-10);
        end
        
        function testShifting(testCase)
            twn = ToroidalWNDistribution([3;5],eye(2));
            chd = CustomHypertoroidalDistribution.fromDistribution(twn);
            s = [1;2];
            chd2 = chd.shift(s);
            twn2 = twn.shift(s);
            [xTest,yTest]=meshgrid(linspace(0,2*pi,100));
            testCase.verifyEqual(chd2.pdf([xTest(:)';yTest(:)']),twn.pdf([xTest(:)' - s(1);yTest(:)' - s(2)]),'AbsTol',1E-10);
            testCase.verifyEqual(chd2.pdf([xTest(:)';yTest(:)']),twn2.pdf([xTest(:)';yTest(:)']),'AbsTol',1E-10);
        end
        
        function testScaling(testCase)
            twn = ToroidalWNDistribution([3;5],eye(2));
            chd = CustomHypertoroidalDistribution.fromDistribution(twn);
            
            testCase.verifyWarningFree(@()chd.scale(2));
            testCase.verifyWarningFree(@()chd.scale(0.5));
            testCase.verifyWarning(@()chd.scale(-0.5),'CustomHypertoroidalDistribution:NegativeScaling');
            
            chdScaled = chd.scale(2);
            [xTest,yTest]=meshgrid(linspace(0,2*pi,100));
            testCase.verifyEqual(chdScaled.pdf([xTest(:)';yTest(:)']),2*twn.pdf([xTest(:)';yTest(:)']));
        end
        
        function testSlicingAndConditioning2DTo1D(testCase)
            twn = ToroidalWNDistribution([3;4],[1,0.8;0.8,1]);
            chd = CustomHypertoroidalDistribution(@(xs)twn.pdf(xs),2);
            for z=[[1;3],zeros(2,1)] % Also test case 0 because it is handled differently
                z1 = z(1);
                z2 = z(2);
                
                cdUnnormDim1 = chd.sliceAt(1,z1);
                cdUnnormDim2 = chd.sliceAt(2,z2);

                grid = linspace(0,2*pi,100);
                testCase.verifyEqual(cdUnnormDim1.pdf(grid),twn.pdf([z1*ones(1,100);grid]));
                testCase.verifyEqual(cdUnnormDim2.pdf(grid),twn.pdf([grid;z2*ones(1,100)]));

                cdNormDim1 = chd.conditionOn(1,z1);
                cdNormDim2 = chd.conditionOn(2,z2);
                normConst1 = 1/integral(@(x)reshape(twn.pdf([z1*ones(1,numel(x));x(:)']),size(x)),0,2*pi);
                normConst2 = 1/integral(@(x)reshape(twn.pdf([x(:)';z2*ones(1,numel(x))]),size(x)),0,2*pi);

                testCase.verifyEqual(cdNormDim1.pdf(grid),normConst1*twn.pdf([z1*ones(1,100);grid]));
                testCase.verifyEqual(cdNormDim2.pdf(grid),normConst2*twn.pdf([grid;z2*ones(1,100)]));
            end
        end
        
        function testSlicingAndConditioning2DTo1DShifted(testCase)
            twn = ToroidalWNDistribution([3;4],[1,0.8;0.8,1]);
            chd = CustomHypertoroidalDistribution(@(xs)twn.pdf(xs),2);
            z1 = 1;
            z2 = 3;
            cdSliceDim1 = chd.shift([-z1;0]).sliceAt(1,0);
            cdSliceDim2 = chd.shift([0;-z2]).sliceAt(2,0);
            
            grid = linspace(0,2*pi,100);
            testCase.verifyEqual(cdSliceDim1.pdf(grid),twn.pdf([z1*ones(1,100);grid]));
            testCase.verifyEqual(cdSliceDim2.pdf(grid),twn.pdf([grid;z2*ones(1,100)]));
            
            cdNormDim1 = chd.shift([-z1;0]).conditionOn(1,0);
            cdNormDim2 = chd.shift([0;-z2]).conditionOn(2,0);
            normConst1 = 1/integral(@(x)reshape(twn.pdf([z1*ones(1,numel(x));x(:)']),size(x)),0,2*pi);
            normConst2 = 1/integral(@(x)reshape(twn.pdf([x(:)';z2*ones(1,numel(x))]),size(x)),0,2*pi);
            
            testCase.verifyEqual(cdNormDim1.pdf(grid),normConst1*twn.pdf([z1*ones(1,100);grid]));
            testCase.verifyEqual(cdNormDim2.pdf(grid),normConst2*twn.pdf([grid;z2*ones(1,100)]));
        end
        
        function testSlicingAndConditioning2DTo1DScaled(testCase)
            twn = ToroidalWNDistribution([3;4],[1,0.8;0.8,1]);
            chd = CustomHypertoroidalDistribution(@(xs)twn.pdf(xs),2);
            z1 = 1;
            z2 = 3;
            scalingFactor1 = 2;
            scalingFactor2 = 0.5;
            cdSliceDim1 = chd.scale(scalingFactor1).sliceAt(1,z1);
            cdSliceDim2 = chd.scale(scalingFactor2).sliceAt(2,z2);
            
            grid = linspace(0,2*pi,100);
            testCase.verifyEqual(cdSliceDim1.pdf(grid),scalingFactor1*twn.pdf([z1*ones(1,100);grid]));
            testCase.verifyEqual(cdSliceDim2.pdf(grid),scalingFactor2*twn.pdf([grid;z2*ones(1,100)]));
        end
        
        function testSlicingAndConditioning3DTo2D(testCase)
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            chd = CustomHypertoroidalDistribution(@(xs)hwn.pdf(xs),hwn.dim);
            z1 = 1;
            z2 = 3;
            z3 = 4;
            cdUnnormDim1 = chd.sliceAt(1,z1);
            cdUnnormDim2 = chd.sliceAt(2,z2);
            cdUnnormDim3 = chd.sliceAt(3,z3);
            
            [mesh1,mesh2] = meshgrid(linspace(0,2*pi,20));
            testCase.verifyEqual(cdUnnormDim1.pdf([mesh1(:)';mesh2(:)']),hwn.pdf([z1*ones(1,numel(mesh1));mesh1(:)';mesh2(:)']));
            testCase.verifyEqual(cdUnnormDim2.pdf([mesh1(:)';mesh2(:)']),hwn.pdf([mesh1(:)';z2*ones(1,numel(mesh1));mesh2(:)']));
            testCase.verifyEqual(cdUnnormDim3.pdf([mesh1(:)';mesh2(:)']),hwn.pdf([mesh1(:)';mesh2(:)';z3*ones(1,numel(mesh1))]));
            
            cdNormDim1 = chd.conditionOn(1,z1);
            cdNormDim2 = chd.conditionOn(2,z2);
            cdNormDim3 = chd.conditionOn(3,z3);
            normConst1 = 1/integral2(@(x,y)reshape(hwn.pdf([ones(1,numel(x));x(:)';y(:)']),size(x)),0,2*pi,0,2*pi);
            normConst2 = 1/integral2(@(x,y)reshape(hwn.pdf([x(:)';z2*ones(1,numel(x));y(:)']),size(x)),0,2*pi,0,2*pi);
            normConst3 = 1/integral2(@(x,y)reshape(hwn.pdf([x(:)';y(:)';z3*ones(1,numel(x))]),size(x)),0,2*pi,0,2*pi);
            
            testCase.verifyEqual(cdNormDim1.pdf([mesh1(:)';mesh2(:)']),normConst1*hwn.pdf([z1*ones(1,numel(mesh1));mesh1(:)';mesh2(:)']));
            testCase.verifyEqual(cdNormDim2.pdf([mesh1(:)';mesh2(:)']),normConst2*hwn.pdf([mesh1(:)';z2*ones(1,numel(mesh1));mesh2(:)']));
            testCase.verifyEqual(cdNormDim3.pdf([mesh1(:)';mesh2(:)']),normConst3*hwn.pdf([mesh1(:)';mesh2(:)';z3*ones(1,numel(mesh1))]));
        end
        
        function testSlicingAndConditioning3DTo1D(testCase)
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            chd = CustomHypertoroidalDistribution(@(xs)hwn.pdf(xs),hwn.dim);
            z1 = [1;5];
            z2 = [2;4];
            z3 = [3;6];
            cdUnnormDim1 = chd.sliceAt([1,2],z1);
            cdUnnormDim2 = chd.sliceAt([1,3],z2);
            cdUnnormDim3 = chd.sliceAt([2,3],z3);
            
            grid = linspace(0,2*pi,100);
            testCase.verifyEqual(cdUnnormDim1.pdf(grid),hwn.pdf([z1(1)*ones(1,numel(grid));z1(2)*ones(1,numel(grid));grid]));
            testCase.verifyEqual(cdUnnormDim2.pdf(grid),hwn.pdf([z2(1)*ones(1,numel(grid));grid;z2(2)*ones(1,numel(grid))]));
            testCase.verifyEqual(cdUnnormDim3.pdf(grid),hwn.pdf([grid;z3(1)*ones(1,numel(grid));z3(2)*ones(1,numel(grid))]));
            
            cdNormDim1 = chd.conditionOn([1,2],z1);
            cdNormDim2 = chd.conditionOn([1,3],z2);
            cdNormDim3 = chd.conditionOn([2,3],z3);
            
            normConst1 = 1/integral(@(x)reshape(hwn.pdf([z1(1)*ones(1,numel(x));z1(2)*ones(1,numel(x));x(:)']),size(x)),0,2*pi);
            normConst2 = 1/integral(@(x)reshape(hwn.pdf([z2(1)*ones(1,numel(x));x(:)';z2(2)*ones(1,numel(x))]),size(x)),0,2*pi);
            normConst3 = 1/integral(@(x)reshape(hwn.pdf([x(:)';z3(1)*ones(1,numel(x));z3(2)*ones(1,numel(x))]),size(x)),0,2*pi);
            
            testCase.verifyEqual(cdNormDim1.pdf(grid),normConst1*hwn.pdf([z1(1)*ones(1,numel(grid));z1(2)*ones(1,numel(grid));grid]));
            testCase.verifyEqual(cdNormDim2.pdf(grid),normConst2*hwn.pdf([z2(1)*ones(1,numel(grid));grid;z2(2)*ones(1,numel(grid))]));
            testCase.verifyEqual(cdNormDim3.pdf(grid),normConst3*hwn.pdf([grid;z3(1)*ones(1,numel(grid));z3(2)*ones(1,numel(grid))]));
        end
        
        function testMarginalizeOut(testCase)        
            twn = ToroidalWNDistribution([3;4],2*[1,0.8;0.8,1]);
            tfd = ToroidalFourierDistribution.fromDistribution(twn,[31,31],'identity');
            chd = CustomHypertoroidalDistribution.fromDistribution(tfd);
            grid = linspace(-pi,3*pi,300);
            
            fdMarg = tfd.marginalizeOut(1);
            cdMarg = chd.marginalizeOut(1);
            testCase.verifyEqual(cdMarg.pdf(grid),fdMarg.pdf(grid),'AbsTol',5e-16);
            
            fdMarg = tfd.marginalizeOut(2);
            cdMarg = chd.marginalizeOut(2);
            testCase.verifyEqual(cdMarg.pdf(grid),fdMarg.pdf(grid),'AbsTol',5e-16);
        end
        
        function testMarginalize3Dto2D(testCase)
            twn = HypertoroidalWNDistribution([3;4;6],2*[1,0.8,0.3;0.8,1,0.5;0.3,0.5,2]);
            hfd = HypertoroidalFourierDistribution.fromDistribution(twn, [13, 13, 13], 'identity');
            chd = CustomHypertoroidalDistribution.fromDistribution(hfd);
            [mesh1,mesh2] = meshgrid(linspace(-pi,3*pi,25));
            
            fdMarg = hfd.marginalizeOut(1);
            cdMarg = chd.marginalizeOut(1);
            testCase.verifyEqual(cdMarg.pdf([mesh1(:)';mesh2(:)']),fdMarg.pdf([mesh1(:)';mesh2(:)']),'AbsTol',5e-16);
            
            fdMarg = hfd.marginalizeOut(2);
            cdMarg = chd.marginalizeOut(2);
            testCase.verifyEqual(cdMarg.pdf([mesh1(:)';mesh2(:)']),fdMarg.pdf([mesh1(:)';mesh2(:)']),'AbsTol',5e-16);
            
            fdMarg = hfd.marginalizeOut(3);
            cdMarg = chd.marginalizeOut(3);
            testCase.verifyEqual(cdMarg.pdf([mesh1(:)';mesh2(:)']),fdMarg.pdf([mesh1(:)';mesh2(:)']),'AbsTol',5e-16);
        end
        
        function testMarginalize3Dto1D(testCase)
            twn = HypertoroidalWNDistribution([3;4;6],2*[1,0.8,0.3;0.8,1,0.5;0.3,0.5,2]);
            hfd = HypertoroidalFourierDistribution.fromDistribution(twn, [13, 13, 13], 'identity');
            chd = CustomHypertoroidalDistribution.fromDistribution(hfd);
            grid = linspace(-pi,3*pi,100);
            
            fdMarg = hfd.marginalizeOut([2,3]);
            cdMarg = chd.marginalizeOut([2,3]);
            testCase.verifyEqual(cdMarg.pdf(grid),fdMarg.pdf(grid),'RelTol',5e-7);
            
            fdMarg = hfd.marginalizeOut([1,3]);
            cdMarg = chd.marginalizeOut([1,3]);
            testCase.verifyEqual(cdMarg.pdf(grid),fdMarg.pdf(grid),'AbsTol',5e-7);
            
            fdMarg = hfd.marginalizeOut([1,2]);
            cdMarg = chd.marginalizeOut([1,2]);
            testCase.verifyEqual(cdMarg.pdf(grid),fdMarg.pdf(grid),'AbsTol',5e-7);
        end
        
        function testLikelihood2DTo1D(testCase)
            twn = ToroidalWNDistribution([3;4],[1,0.8;0.8,1]);
            chd = CustomHypertoroidalDistribution(@(xs)twn.pdf(xs),2);
            xFix = 1;
            yFix = 3;
            likelihoodDim1 = chd.likelihood(1,xFix);
            likelihoodDim2 = chd.likelihood(2,yFix);
            
            funSlicedDim1=@(y)twn.pdf([xFix*ones(1,size(y,2));y]); % Sliced
            funMarginalDim1=@(y)arrayfun(@(yCurr)integral(@(x)reshape(twn.pdf([x(:)';yCurr*ones(1,numel(x))]),size(x)),0,2*pi),y); % Marginal
            likelihoodDim1Fun = @(x)funSlicedDim1(x)./funMarginalDim1(x);
            funSlicedDim2=@(x)twn.pdf([x;yFix*ones(1,size(x,2))]); % Sliced
            funMarginalDim2=@(x)arrayfun(@(xCurr)integral(@(y)reshape(twn.pdf([xCurr*ones(1,numel(y));y(:)']),size(y)),0,2*pi),x); % Marginal
            likelihoodDim2Fun = @(x)funSlicedDim2(x)./funMarginalDim2(x);
            
            grid = linspace(-pi,3*pi,100);
            testCase.verifyEqual(likelihoodDim1.pdf(grid),likelihoodDim1Fun(grid));
            testCase.verifyEqual(likelihoodDim2.pdf(grid),likelihoodDim2Fun(grid));
        end
        
        function testLikelihood3DTo2D(testCase)
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            chd = CustomHypertoroidalDistribution(@(xs)hwn.pdf(xs),hwn.dim);
            xFix = 1;
            yFix = 3;
            zFix = 4;
            chdLikelihoodDim1 = chd.likelihood(1,xFix);
            chdLikelihoodDim2 = chd.likelihood(2,yFix);
            chdLikelihoodDim3 = chd.likelihood(3,zFix);
            
            funSlicedDim1=@(y,z)hwn.pdf([xFix*ones(1,size(y,2));y;z]); % Sliced
            funSlicedDim2=@(x,z)hwn.pdf([x;yFix*ones(1,size(x,2));z]); % Sliced
            funSlicedDim3=@(x,y)hwn.pdf([x;y;zFix*ones(1,size(x,2))]); % Sliced
            funMarginalDim1=@(y,z)arrayfun(@(yCurr,zCurr)integral(@(x)reshape(hwn.pdf([x(:)';yCurr*ones(1,size(x,2));zCurr*ones(1,size(x,2))]),size(x)),0,2*pi),y,z); % Marginal
            funMarginalDim2=@(x,z)arrayfun(@(xCurr,zCurr)integral(@(y)reshape(hwn.pdf([xCurr*ones(1,size(y,2));y(:)';zCurr*ones(1,size(y,2))]),size(y)),0,2*pi),x,z); % Marginal
            funMarginalDim3=@(x,y)arrayfun(@(xCurr,yCurr)integral(@(z)reshape(hwn.pdf([xCurr*ones(1,size(z,2));yCurr*ones(1,size(z,2));z(:)']),size(z)),0,2*pi),x,y); % Marginal
            likelihoodDim1Fun = @(y,z)funSlicedDim1(y,z)./funMarginalDim1(y,z);
            likelihoodDim2Fun = @(x,z)funSlicedDim2(x,z)./funMarginalDim2(x,z);
            likelihoodDim3Fun = @(x,y)funSlicedDim3(x,y)./funMarginalDim3(x,y);
            
            [mesh1,mesh2] = meshgrid(linspace(0,2*pi,20));
            testCase.verifyEqual(chdLikelihoodDim1.pdf([mesh1(:)';mesh2(:)']),likelihoodDim1Fun(mesh1(:)',mesh2(:)'));
            testCase.verifyEqual(chdLikelihoodDim2.pdf([mesh1(:)';mesh2(:)']),likelihoodDim2Fun(mesh1(:)',mesh2(:)'));
            testCase.verifyEqual(chdLikelihoodDim3.pdf([mesh1(:)';mesh2(:)']),likelihoodDim3Fun(mesh1(:)',mesh2(:)'));
        end
        
        function testLikelihood3DTo1D(testCase)
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            chd = CustomHypertoroidalDistribution(@(xs)hwn.pdf(xs),hwn.dim);
            xyFix = [1;5];
            xzFix = [2;4];
            yzFix = [3;6];
            
            chdLikelihoodDim12 = chd.likelihood([1,2],xyFix);
            chdLikelihoodDim13 = chd.likelihood([1,3],xzFix);
            chdLikelihoodDim23 = chd.likelihood([2,3],yzFix);
            
            funSlicedDim12=@(z)hwn.pdf([xyFix.*ones(1,size(z,2));z]); % Sliced
            funSlicedDim13=@(y)hwn.pdf([xzFix(1)*ones(1,size(y,2));y;xzFix(2)*ones(1,size(y,2))]); % Sliced
            funSlicedDim23=@(x)hwn.pdf([x;yzFix.*ones(1,size(x,2))]); % Sliced
            funMarginalDim12=@(z)arrayfun(@(zCurr)integral2(@(x,y)reshape(hwn.pdf([x(:)';y(:)';zCurr*ones(1,numel(x))]),size(x)),0,2*pi,0,2*pi),z);
            funMarginalDim13=@(y)arrayfun(@(yCurr)integral2(@(x,z)reshape(hwn.pdf([x(:)';yCurr*ones(1,numel(x));z(:)']),size(x)),0,2*pi,0,2*pi),y);
            funMarginalDim23=@(x)arrayfun(@(xCurr)integral2(@(y,z)reshape(hwn.pdf([xCurr*ones(1,numel(y));y(:)';z(:)']),size(y)),0,2*pi,0,2*pi),x);
            likelihoodDim12Fun = @(z)funSlicedDim12(z)./funMarginalDim12(z);
            likelihoodDim13Fun = @(y)funSlicedDim13(y)./funMarginalDim13(y);
            likelihoodDim23Fun = @(x)funSlicedDim23(x)./funMarginalDim23(x);
            
            grid = linspace(-pi,3*pi,100);
            testCase.verifyEqual(chdLikelihoodDim12.pdf(grid),likelihoodDim12Fun(grid));
            testCase.verifyEqual(chdLikelihoodDim13.pdf(grid),likelihoodDim13Fun(grid));
            testCase.verifyEqual(chdLikelihoodDim23.pdf(grid),likelihoodDim23Fun(grid));
        end
    end
end
