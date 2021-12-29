classdef LinearDiracDistribution < AbstractLinearDistribution & AbstractDiracDistribution
    
    methods
        function this = LinearDiracDistribution(d_, w_)
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (1 x L)
            %       Dirac locations in [0,2pi)
            %   w_ (1 x L)
            %       weights for each Dirac
            arguments
                d_ (:,:) double
                w_ (1,:) double = ones(1,size(d_,2))/size(d_,2);
            end
            this@AbstractDiracDistribution(d_, w_);
        end
        
        function mu = mean(this)
            mu = this.d*this.w';
        end
        
        function C = covariance(this)
            [~,C] = LinearDiracDistribution.weightedSamplesToMeanAndCov(this.d,this.w);
        end

        function result = integral(this)
            arguments
                this (1,1) HypersphericalDiracDistribution
            end
            result = integral@AbstractDiracDistribution(this);
        end
        
        function h = plot(this, varargin)
            switch this.dim
                case 1
                    h = stem(this.d, this.w, varargin{:});
                case 2
                    h = scatter(this.d(1,:),this.d(2,:), this.w/max(this.w)*100, varargin{:});
                case 3
                    h = scatter3(this.d(1,:),this.d(2,:), this.w/max(this.w)*100, varargin{:});
                otherwise
                    error('Plotting not supported for this dimension');
            end
        end
    end
    
    methods(Static)
        function f = fromDistribution(distribution, noOfSamples)
            arguments
                distribution (1,1) AbstractLinearDistribution
                noOfSamples (1,1) {mustBePositive, mustBeInteger}
            end
            f = LinearDiracDistribution(...
                    distribution.sample(noOfSamples),ones(1,noOfSamples)/noOfSamples);
        end
        
        function [mu,C] = weightedSamplesToMeanAndCov(samples, w)
            arguments
                samples (:,:) double
                w (1,:) double = ones(1,size(samples,2))/size(samples,2)
            end
            mu = samples*w';
            if nargout==2
                samplesZMWeighted = (samples-mu).*sqrt(w);
                C = samplesZMWeighted*samplesZMWeighted';
            end
        end
    end
end

