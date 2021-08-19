classdef GaussianMixtureDistribution < AbstractMixture & AbstractLinearDistribution
    % Chosen this name to prevenet collision with GaussianMixture of the
    % Nonlinear Estimation Toolbox
    methods
        function this = GaussianMixtureDistribution(dists, w)
            % Constructor
            assert(isa(dists,'cell') && all(cellfun(@(dist)isa(dist,'AbstractLinearDistribution'),dists)),...
                'dists must be a cell array of Gaussian distributions');
            this = this@AbstractMixture(dists, w);
        end
        
        function gauss = toGaussian(this)
            arguments
                this (1,1) GaussianMixtureDistribution
            end
            gaussArray = [this.dists{:}];
            [mu, C] = mixtureParametersToGaussianParamters([gaussArray.mu],cat(3,gaussArray.C),this.w);
            gauss = GaussianDistribution(mu,C);
        end
        
        function mu = mean(this)
            arguments
                this (1,1) GaussianMixtureDistribution
            end
            gaussArray = [this.dists{:}];
            mu = [gaussArray.mu]*this.w';
        end
        
        function C = covariance(this)
            arguments
                this (1,1) GaussianMixtureDistribution
            end
            gaussArray = [this.dists{:}];
            [~, C] = GaussianMixtureDistribution.mixtureParametersToGaussianParamters([gaussArray.mu],cat(3,gaussArray.C),this.w);
        end
        
    end
    
    methods (Static)
        function [mu, C] = mixtureParametersToGaussianParamters(means, covarianceMatrices, weights)
            arguments
                means (:,:) double
                covarianceMatrices (:,:,:) double
                weights = ones(1,size(means,2))/size(means,2)
            end
            if nargout==1
                mu = means*weights';
            else 
                CfromCov = sum(covarianceMatrices.*reshape(weights, [1,1,size(means,2)]),3);
                [mu,CfromMeans] = LinearDiracDistribution.weightedSamplesToMeanAndCov(means, weights);
                C = CfromCov + CfromMeans;
            end
        end
    end
end