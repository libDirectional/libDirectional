classdef CustomHypercylindricalDistributionTest < matlab.unittest.TestCase
    methods (Test)                                
        function testConstructor(testCase)
            mat = 0.001*(magic(7)*magic(7)');
            mat = mat+mat';
            rbd = HypercylindricalWNDistribution([2;3;4;5;6;7],mat(2:end,2:end),3);
                
            chd = CustomHypercylindricalDistribution(@(x)rbd.pdf(x),3,3);
            [a,b,c,d,e,f] = ndgrid(-3:3,-3:3,-2:2,-2:2,-2:2,-2:2);
            testCase.verifyEqual(rbd.pdf([a(:)';b(:)';c(:)';d(:)';e(:)';f(:)']),...
                chd.pdf([a(:)';b(:)';c(:)';d(:)';e(:)';f(:)']));
        end
        
        function testFromDistribution(testCase)
            mat = 0.001*(magic(7)*magic(7)');
            mat = mat+mat';
            rbd = HypercylindricalWNDistribution([2;3;4;5;6;7],mat(2:end,2:end),3);
                
            chd = CustomHypercylindricalDistribution.fromDistribution(rbd);
            [a,b,c,d,e,f] = ndgrid(-3:3,-3:3,-2:2,-2:2,-2:2,-2:2);
            testCase.verifyEqual(rbd.pdf([a(:)';b(:)';c(:)';d(:)';e(:)';f(:)']),...
                chd.pdf([a(:)';b(:)';c(:)';d(:)';e(:)';f(:)']));
        end
        
        function testConditionOnLinear(testCase)
            vm = VMDistribution(0,1);
            gauss = GaussianDistribution([1;2],eye(2));
            fun = @(x)reshape(vm.pdf(x(1,:)).*gauss.pdf([x(2,:);x(3,:)]),[1,size(x,2)]);
            
            chcd = CustomHypercylindricalDistribution(fun,1,2);
            chtd = chcd.conditionOnLinear([2;1]);
            
            x = linspace(0,2*pi,100);
            testCase.verifyEqual(chtd.pdf(x), vm.pdf(x),'AbsTol',1e-15);
        end
        
        function testConditionOnPeriodic(testCase)
            vm = VMDistribution(0,1);
            gauss = GaussianDistribution([1;2],eye(2));
            fun = @(x)reshape(vm.pdf(x(1,:)).*gauss.pdf([x(2,:);x(3,:)]),[1,size(x,2)]);
            
            chcd = CustomHypercylindricalDistribution(fun,1,2);
            cld = chcd.conditionOnPeriodic(1);
            
            [a,b] = ndgrid(-3:3,-3:3);
            testCase.verifyEqual(cld.pdf([a(:)';b(:)']),gauss.pdf([a(:)';b(:)']),'AbsTol',1e-10);
        end
    end
   
end