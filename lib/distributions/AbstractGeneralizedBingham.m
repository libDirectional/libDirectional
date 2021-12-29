classdef (Abstract) AbstractGeneralizedBingham < AbstractLinPeriodicDistribution
    % Distribution on S^(d-1) x R^n.
    %  
    
    methods
        function h = plot(this)
            assert(this.linD == 1);
            assert(this.boundD == 2);
            
            phi = linspace(0,2*pi,100);
            l = linspace(-5,5, 100);
            [Phi,L] = meshgrid(phi,l);
            C = this.pdf([cos(Phi(:))'; sin(Phi(:))'; L(:)']);
            h = surf(cos(Phi),sin(Phi),L, reshape(C, size(Phi)));
            shading interp
        end
        
        function result = integralNumerical(this)
            assert(this.linD == 1);
            assert(this.boundD == 2);
            
            f = @(x,y) reshape(this.pdf([cos(x(:)');sin(x(:)');y(:)']), size(x));
            result = integral2(f, 0, 2*pi, -10, 10);            
        end
        
        function l = logLikelihood(this, samples, weights)
            % Calculates the log-likelihood of the given samples
            %
            % Parameters:
            %   samples (dim x n)
            %       n samples on the torus
            %   weights (1 x n)
            %       weight for each sample (by default, all weights are 1)
            % Returns:
            %   l (scalar)
            %       log-likelihood of obtaining the given samples
            assert(size(samples,1)==this.boundD + this.linD);
            assert(size(samples,2)>=1);
            
            if nargin>2
                assert(size(weights,1)==1);
                assert(size(weights,2)==size(samples,2));
                l = sum(weights.*log(this.pdf(samples)));
            else
                l = sum(log(this.pdf(samples)));
            end            
        end        
    end
    
end

