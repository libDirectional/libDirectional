classdef HypersphericalMixtureTest < matlab.unittest.TestCase
    properties
    end
    
    methods (Test)                                
        function testPdf3D(testCase)
            wad=WatsonDistribution([0;0;1],-10);
            vmf=VMFDistribution([0;0;1],1);
            w=[0.3,0.7];
            smix=HypersphericalMixture({wad,vmf},w);
            
            [phi,theta]=meshgrid(linspace(0,2*pi,10),linspace(-pi/2,pi/2,10));
            [x,y,z]=sph2cart(phi(:)',theta(:)',1);
            testCase.verifyEqual(smix.pdf([x;y;z]),w(1)*wad.pdf([x;y;z])+w(2)*vmf.pdf([x;y;z]),'AbsTol',1E-10);
        end
        
        function testPdf4D(testCase)
            wad=WatsonDistribution([0;0;0;1],-10);
            vmf=VMFDistribution([0;1;0;0],1);
            w=[0.3,0.7];
            smix=HypersphericalMixture({wad,vmf},w);
            
            [a,b,c,d]=ndgrid(linspace(-1,1,4));
            points=[a(:)';b(:)';c(:)';d(:)'];
            points=points./repmat(sqrt(sum(points.^2,1)),4,1);
            testCase.verifyEqual(smix.pdf(points),w(1)*wad.pdf(points)+w(2)*vmf.pdf(points),'AbsTol',1E-10);
        end
    end    
   
end
