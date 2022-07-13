classdef HypercylindricalStateSpaceSubdivisionGaussianDistributionTest < matlab.unittest.TestCase
    methods (Test)                                
        function testFromDistributionUncorrelated(testCase)
            n = 100;
            vm = VMDistribution(0,1);
            gauss = GaussianDistribution(0,1);
            customUncorrelated = CustomHypercylindricalDistribution(...
                @(x)vm.pdf(x(1,:)).*gauss.pdf(x(2,:)),1,1);
            
            hcrbd = HypercylindricalStateSpaceSubdivisionGaussianDistribution.fromDistribution(...
                customUncorrelated,n);
            
            [xMesh, yMesh] = meshgrid(linspace(-2,2,100),linspace(-pi,3*pi,100));
            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)']),...
                customUncorrelated.pdf([xMesh(:)';yMesh(:)']),'AbsTol',1e-13);            
        end
        
        function testFromDistributionWNLittleCorrelation(testCase)
            n = 500;
            hcwn = HypercylindricalWNDistribution([1;1],[5,0.1;0.1,2],1);
            % This approximation is only good because there is little
            % correlation (otherwise, the tails of the Gaussian could wrap
            % around the cylinder
            hcrbd = HypercylindricalStateSpaceSubdivisionGaussianDistribution.fromDistribution(hcwn,n);
            
            [xMesh, yMesh] = meshgrid(linspace(-2,2,100),linspace(-pi,3*pi,100));
            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)']),...
                hcwn.pdf([xMesh(:)';yMesh(:)']),'AbsTol',5e-6);
        end
        
        function testFromDistributionWNHigherCorrelation(testCase)
            n = 500;
            hcwn = HypercylindricalWNDistribution([1;1],[1,0.2;0.2,2],1);
            % This get more and more tricky sind the wrap around part will
            % become significant, making it impossible to approximate by a
            % Gaussian
            hcrbd = HypercylindricalStateSpaceSubdivisionGaussianDistribution.fromDistribution(hcwn,n);
            
            [xMesh, yMesh] = meshgrid(linspace(-2,2,100),linspace(-pi,3*pi,100));
            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)']),...
                hcwn.pdf([xMesh(:)';yMesh(:)']),'AbsTol',5e-5);
        end
        
        function testLinearMean(testCase)
            n = 101;
            rbd = HypercylindricalStateSpaceSubdivisionGaussianDistribution(...
                FIGDistribution.fromDistribution(...
                CircularUniformDistribution,n),...
                arrayfun(@(i)GaussianDistribution(i,1),(0:n-1)'));
            % Because we are only testing the linear mean and all are
            % equally weighted, we should land on 50!
            testCase.verifyEqual(rbd.linearMean,50,'AbsTol',1E-13);
        end
        
        function testMultiplyS1xR1hwnWithVisualization(testCase)
            n = 500;
            visualizeDensityAndDifference = false;

            hwn1 = HypercylindricalWNDistribution(pi+[0;1],[1,-0.2;-0.2,2],1);
            hwn2 = HypercylindricalWNDistribution([1;0],[1,0.2;0.2,2],1);

            chdFus = CustomHypercylindricalDistribution(@(x)hwn1.pdf(x).*hwn2.pdf(x),1,1);
            chdFus = chdFus.normalize();
            chdFus = chdFus.normalize();

            apd1 = HypercylindricalStateSpaceSubdivisionGaussianDistribution.fromDistribution(hwn1,n);
            apd2 = HypercylindricalStateSpaceSubdivisionGaussianDistribution.fromDistribution(hwn2,n);
            apdUp = apd1.multiply(apd2);

            [xMesh, yMesh] = meshgrid(linspace(-2,2,100),linspace(-pi,3*pi,100));
            testCase.verifyEqual(apd1.pdf([xMesh(:)';yMesh(:)']),hwn1.pdf([xMesh(:)';yMesh(:)']),'AbsTol',5e-5);
            testCase.verifyEqual(apd2.pdf([xMesh(:)';yMesh(:)']),hwn2.pdf([xMesh(:)';yMesh(:)']),'AbsTol',5e-5);
            
            if visualizeDensityAndDifference
                chdtmp1 = CustomHypercylindricalDistribution(@(x)hwn1.pdf(x)-apd1.pdf(x),1,1); %#ok<UNRCH>
                chdtmp2 = CustomHypercylindricalDistribution(@(x)hwn2.pdf(x)-apd2.pdf(x),1,1);
                chdtmpUp = CustomHypercylindricalDistribution(@(x)chdFus.pdf(x)-apdUp.pdf(x),1,1);

                figure(1)
                hwn1.plotCylinder()
                figure(2)
                hwn2.plotCylinder()
                figure(3)
                chdFus.plotCylinder()

                figure(11)
                apd1.plotCylinder()
                zlim1=zlim();
                figure(12)
                apd2.plotCylinder()
                zlim2=zlim();
                figure(13)
                apdUp.plotCylinder()
                zlim3=zlim();

                figure(21)
                chdtmp1.plotCylinder(zlim1)
                figure(22)
                chdtmp2.plotCylinder(zlim2)
                figure(23)
                chdtmpUp.plotCylinder(zlim3)
            end
            
            testCase.verifyEqual(apdUp.pdf([xMesh(:)';yMesh(:)']),chdFus.pdf([xMesh(:)';yMesh(:)']),'AbsTol',0.0003);
        end
        
    end
   
end
