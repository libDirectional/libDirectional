classdef (Abstract) AbstractComplexHypersphericalDistribution < AbstractDistribution
    % Abstract base class for distributions on the hypershere (Sd)
    
    methods
        function m = mean(~) %#ok<STOUT> 
            error('Mean currently not supported.')
        end
        function i = integral(this, n)
            % Calculate integral of pdf to check normalization
            % (should always be 1)
            % Returns:
            %   i (scalar)
            %       integral over hypersphere surface of pdf (uses
            %       approximation, not very accurate for higher dimensions)
            
            if nargin <2
                %use monte carlo integration
                n = 100000; %number of samples for integration
            end
            Z = complex(randn(this.dim, n), randn(this.dim, n));
            Z = bsxfun(@rdivide, Z, sqrt(sum(Z .* conj(Z), 1)));
            p = this.pdf(Z);
            Sd = AbstractComplexHypersphericalDistribution.computeUnitSphereSurface(this.dim);
            i = sum(p)/n * Sd; %average value of pdf times surface area of unit sphere
        end
        
    end
    
    methods (Static)
        function surfaceArea = computeUnitSphereSurface(dimension)
            % Computes surface area of (d-1)-sphere
            % see http://en.wikipedia.org/wiki/N-sphere#Volume_and_surface_area
            % dimension = 2 => circle
            % dimension = 3 => sphere
            % ...
            %
            % Parameters:
            %   dimension (scalar)
            %       the dimension of the sphere
            % Returns:
            %   surfaceArea (scalar)
            %       the surface area of the sphere
            surfaceArea = 2*pi^dimension / factorial(dimension-1);
        end
        
        function [h1, h2] = scatter(Z, labels, normalized, varargin)
            if normalized
                Z = bsxfun(@times, Z, exp(-1i * angle(Z(1, :))));
            end
            
            N = size(Z, 2);
            
            clf('reset'); hold on;
            
            colors = get(gca, 'ColorOrder');
            
            h1 = zeros(1,N);
            h2 = zeros(1,N);
            
            for i = 1:N
                if isempty(labels)
                    color = colors(mod(i - 1, size(colors, 1)) + 1, :);
                else
                    color = colors(mod(labels(i) - 1, size(colors, 1)) + 1, :);
                end
                h1(i) = plot(real(Z(:, i)), imag(Z(:, i)), 'color', color, varargin{:});
                h2(i) = plot(real(Z(1, i)), imag(Z(1, i)), 'o', 'color', color, varargin{:});
            end
            
            axis equal; grid on; hold off; box on;
            
            xlabel('Real(Z)');
            ylabel('Imag(Z)');
        end
    end
end