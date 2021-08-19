classdef HypertoroidalMixture < AbstractHypertoroidalDistribution & AbstractMixture
    % Mixture of multiple hypertoroidal distributions. The distributions may belong
    % to different classes.
    
    methods
        function this = HypertoroidalMixture(dists, w)
            % Constructor
            assert(isa(dists,'cell') && all(cellfun(@(dist)isa(dist,'AbstractHypertoroidalDistribution'),dists)),...
                'dists must be a cell array of hypertoroidal distributions');
            if all(cellfun(@(dist)isa(dist,'HypertoroidalFourierDistribution'),dists))
                warning('Mixtures of HypertoroidalFourierDistributions can be built by combining the Fourier coefficients so using a mixture may not be necessary.');
            end
            this = this@AbstractMixture(dists, w);
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
            for i=1:length(this.dists)
                m = m + this.w(i)*this.dists{i}.trigonometricMoment(n); % Calculate moments using moments of each component
            end
        end
        
        function hd = shift(this, shiftAngles)
            % Shifts the distribution by shiftAngles.
            %
            % Parameters:
            %   shiftAngles (dim x 1 column vector)
            %       angles to shift by
            % Return:
            %   hd (HypertoroidalMixture)
            %       shifted distribution
            arguments
                this (1,1) HypertoroidalMixture
                shiftAngles (:,1) double
            end
            assert(size(shiftAngles,1)==this.dim);
            hd = this;
            hd.dists = cellfun(@(dist){dist.shift(shiftAngles)},hd.dists);
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
end

