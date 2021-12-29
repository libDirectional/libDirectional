classdef HyperhemisphericalDiracDistribution < AbstractHypersphereSubsetDiracDistribution & AbstractHyperhemisphericalDistribution
    % Dirac distribution on the hypersphere with dirac positions d and
    % weights w.
    
    methods
        function this = HyperhemisphericalDiracDistribution(d_, w_)
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (dim x L)
            %       Dirac locations on the unit hypersphere
            %   w_ (1 x L)
            %       weights for each Dirac
            arguments
                d_ (:,:) double {mustBeNonempty}
                w_ (1,:) double = ones(1,size(d_,2))/size(d_,2);
            end
            assert(all(d_(end,:)>=0), 'DiracLocation:notOnUpperHemisphere'); % Assert there are only Diracs on the upper hemisphere
            this@AbstractHypersphereSubsetDiracDistribution(d_, w_);
        end      
        
        function C = moment(this)
            C = this.d.*this.w*this.d';
        end
        
        function p = plot(this, varargin)
            % Create an appropriate plot
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to plot/surf command
            % Returns:
            %   p (scalar)
            %       plot handle
            switch this.dim
                case 2
                    % use polar coordinates, plot angle->pdf(angle)
                    p = stem(atan2(this.d(2,:), this.d(1,:)), this.w, varargin{:});
                case 3
                    % plot points on hemisphere
                    p = scatter3(this.d(1,:),this.d(2,:),this.d(3,:), this.w*size(this.d,2)*20, varargin{:});
                otherwise
                    error('Plotting for this dimension is currently not supported');
            end                        
        end
    end
    methods (Static)
        function f = fromDistribution(distribution, noOfSamples)
            arguments
                distribution (1,1) AbstractHyperhemisphericalDistribution
                noOfSamples (1,1) {mustBePositive, mustBeInteger}
            end
            f = HyperhemisphericalDiracDistribution(...
                distribution.sample(noOfSamples),ones(1,noOfSamples)/noOfSamples);
        end
    end
end

