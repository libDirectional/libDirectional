classdef HypersphericalMixture < AbstractHypersphericalDistribution
    % Mixture of multiple hyperspherical distributions. The distributions may belong
    % to different classes.
    
    properties
        dists     % Cell array of hyperspherical distributions
        w         % Weights
    end
    
    methods
        function this = HypersphericalMixture(dists, w)
            % Constructor
            assert(isa(dists,'cell') && all(cellfun(@(dist)isa(dist,'AbstractHypersphericalDistribution'),dists)),...
                'dists must be a cell array of hypertoroidal distributions');
            assert(all(size(dists) == size(w)),'size of dists and w must be equal');
            assert(all(dists{1}.dim==cellfun(@(dist)dist.dim,dists))); % Ensure equal dimensions
            assert(all(w >= 0), 'weights must be nonnegative');
            if all(cellfun(@(dist)isa(dist,'AbstractSphericalHarmonicDistribution'),dists))
                warning('Creating a mixture of Spherical Harmonics may not be necessary.');
            end
            this.dim=dists{1}.dim;
            this.dists=dists;
            this.w=w/sum(w);
        end
        
        function p = pdf(this, xa)
            % Evaluate pdf at each column of xa
            %
            % Parameters:
            %   xa (dim x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            assert(size(xa,1)==this.dim);
            p = zeros(1, size(xa,2));
            for i=1:length(this.dists);
                p = p + this.w(i)*this.dists{i}.pdf(xa); % Calculate pdf using individual pdfs
            end
        end
    end
    methods (Sealed)
        function s = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       n samples on the dim-dimensional hypersphere
            
            % Sample component first, then sample from the chosen component
            d = discretesample(this.w,n);
            s = zeros(this.dim,n);
            occurrences=histc(d,1:length(this.dists));
            count=1;
            for i=1:length(this.dists)
                s(:,count:count+occurrences(i)-1) = this.dists{i}.sample(occurrences(i));
                count=count+occurrences(i);
            end
            [~,order]=sort(d);
            s(:,order)=s;
        end
    end
end

