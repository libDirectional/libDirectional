classdef SphericalHarmonicFilterTest < matlab.unittest.TestCase
   
    methods (Test)
        function testUpdateIdentity(testCase)
            shdFilter=SphericalHarmonicFilter(30);
            vmfFilter=VMFFilter();
            
            vmf1=VMFDistribution([0;1;0],1);
            vmf2=VMFDistribution([0;0;1],0.1);
            shd1=SphericalHarmonicDistributionComplex.fromDistributionFast(vmf1,30);
            shd2=SphericalHarmonicDistributionComplex.fromDistributionFast(vmf2,30);
            
            vmfFilter.setState(vmf1);
            vmfFilter.updateIdentity(vmf2,[1;0;0]);

            shdFilter.setState(shd1);
            shdFilter.updateIdentity(shd2,[1;0;0]);
            
            testCase.verifyEqual(vmfFilter.getEstimateMean,shdFilter.getEstimateMean,'AbsTol',1E-10);
        end
        function testUpdateUsingLikelihood(testCase)
            shFilter=SphericalHarmonicFilter(30); % Initialize with uniform
            rng(1)
            posTrue=-1/sqrt(3)*ones(3,1);
            %% generate measurements according to truncated gaussian along x and y axis
            measX=NaN(1,5);
            sigmaX=.3;
            measY=NaN(1,5);
            sigmaY=.3; 
            measZ=NaN(1,5);
            sigmaZ=.3; 

            for i=1:numel(measX)
                measX(i)=normrnd(posTrue(1),sigmaX);
                if measX(i)>1 || measX(i)<-1
                    i=i-1; %#ok<FXSET>
                end 
            end
            for i=1:numel(measY)
                measY(i)=normrnd(posTrue(2),sigmaY);
                if measY(i)>1 || measY(i)<-1
                    i=i-1; %#ok<FXSET>
                end 
            end
            for i=1:numel(measZ)
                measZ(i)=normrnd(posTrue(3),sigmaZ);
                if measZ(i)>1 || measZ(i)<-1
                    i=i-1; %#ok<FXSET>
                end 
            end
            %% Update first using measurements of the x-axis, then of the y-axis and then of the z-axis
            for i=1:numel(measX)
                shFilter.updateNonlinearUsingLikelihood(@(z,x)normpdf(z(1),x(1,:),sigmaX),[measX(i);0;0])
            end
            for i=1:numel(measY)
                shFilter.updateNonlinearUsingLikelihood(@(z,x)normpdf(z(2),x(2,:),sigmaY),[0;measY(i);0])
            end
            for i=1:numel(measZ)
                shFilter.updateNonlinearUsingLikelihood(@(z,x)normpdf(z(3),x(3,:),sigmaZ),[0;0;measZ(i)])
            end
            testCase.verifyEqual(shFilter.getEstimateMean,posTrue,'AbsTol',.3);
        end
        function testUpdateUsingLikelihoodMultiple(testCase)
            % Test if approximately equal to update sequentially or at
            % once.
            shFilter1=SphericalHarmonicFilter(10);
            shFilter2=SphericalHarmonicFilter(10);

            sigmaX=0.3;sigmaY=0.3;sigmaZ=0.3;
            shFilter1.updateNonlinearUsingLikelihood(@(z,x)normpdf(z(1),x(1,:),sigmaX),[-1/sqrt(3);0;0]);
            shFilter1.updateNonlinearUsingLikelihood(@(z,x)normpdf(z(2),x(2,:),sigmaX),[0;-1/sqrt(3);0]);
            shFilter1.updateNonlinearUsingLikelihood(@(z,x)normpdf(z(3),x(3,:),sigmaX),[0;0;-1/sqrt(3)]);

            shFilter2.updateNonlinearUsingLikelihood({@(z,x)normpdf(z(1),x(1,:),sigmaX)...
                @(z,x)normpdf(z(2),x(2,:),sigmaY),@(z,x)normpdf(z(3),x(3,:),sigmaZ)},...
                {[-1/sqrt(3);0;0],[0;-1/sqrt(3);0],[0;0;-1/sqrt(3)]});
            testCase.verifyEqual(shFilter2.getEstimateMean,shFilter1.getEstimateMean,'AbsTol',1E-5);
        end
    end
end

