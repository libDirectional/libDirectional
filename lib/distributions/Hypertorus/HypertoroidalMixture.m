classdef HypertoroidalMixture < AbstractHypertoroidalDistribution
    % Mixture of multiple hypertoroidal distributions. The distributions may belong
    % to different classes.
    
    properties
        dists     % Cell array of hypertoroidal distributions
        w         % Weights
    end
    
    methods
        function this = HypertoroidalMixture(dists, w)
            % Constructor
            assert(isa(dists,'cell') && all(cellfun(@(dist)isa(dist,'AbstractHypertoroidalDistribution'),dists)),...
                'dists must be a cell array of hypertoroidal distributions');
            assert(all(size(dists) == size(w)),'size of cds and w must be equal');
            assert(all(dists{1}.dim==cellfun(@(dist)dist.dim,dists))); % Ensure equal dimensions
            assert(all(w >= 0), 'weights must be nonnegative');
            if all(cellfun(@(cd)isa(cd,'HypertoroidalFourierDistribution'),dists))
                warning('Mixtures of HypertoroidalFourierDistributions can be built by combining the Fourier coefficients so using a mixture may not be necessary.');
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
        
        function m = trigonometricMoment(this,n)
            % Calculate n-th trigonometric moment analytically
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (scalar)
            %       n-th trigonometric moment (complex number)
            m = zeros(this.dim,1);
            for i=1:length(this.dists);
                m = m + this.w(i)*this.dists{i}.trigonometricMoment(n); % Calculate moments using moments of each component
            end
        end
        
        function cm = toCircularMixture(this)
            % Convert to a circular mixture (only in 1D case)
            %
            % Returns:
            %   cm (CircularMixture)
            %       CircularMixture with same parameters
            assert(this.dim == 1);
            % This also requires that all mixture components are circular
            % distributions rather than hypertoroidal distributions!
            cm = CircularMixture(this.dists, this.w);
        end
        
        function tm = toToroidalMixture(this)
            % Convert to a toroidal mixture (only in 2D case)
            %
            % Returns:
            %   tm (ToroidalMixture)
            %       ToroidalMixture with same parameters
            assert(this.dim == 2);
            % This also requires that all mixture components are toroidal
            % distributions rather than hypertoroidal distributions!
            tm = ToroidalMixture(this.dists, this.w);
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
            %       n samples on the dim-dimensional hypertorus
            
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

