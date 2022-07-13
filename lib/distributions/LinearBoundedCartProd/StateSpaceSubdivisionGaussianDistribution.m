classdef StateSpaceSubdivisionGaussianDistribution < StateSpaceSubdivisionDistribution
	
    methods
        function this = StateSpaceSubdivisionGaussianDistribution(gd_, gaussians)
            arguments
                gd_ (1,1) AbstractGridDistribution
                gaussians (:,1) GaussianDistribution
            end
            this@StateSpaceSubdivisionDistribution(gd_, gaussians);
        end      
        
        function d = marginalizePeriodic(this)
            d = GaussianMixtureDistribution(mat2cell(this.linearDistributions,ones(1,numel(this.linearDistributions)),1), this.gd.gridValues/sum(this.gd.gridValues)); %#ok<MMTC> 
        end

        function sed = multiply(this, other)
            arguments
                this (1,1) StateSpaceSubdivisionGaussianDistribution
                other (1,1) StateSpaceSubdivisionGaussianDistribution
            end
            assert(isequal(this.gd.getGrid(),other.gd.getGrid()),...
                'SEdRBDistribution:IncompatibleGrid','Can only multiply for equal grids.');
            sed = this;
            % We need the likelihoods here
            factorsLinear = mvnpdf([this.linearDistributions.mu]',[other.linearDistributions.mu]',cat(3,this.linearDistributions.C)+cat(3,other.linearDistributions.C));
            sed.gd.gridValues = this.gd.gridValues.*other.gd.gridValues.*factorsLinear;
            sed.gd = sed.gd.normalize(warnUnnorm=false);
            for i = 1:numel(this.linearDistributions)
                sed.linearDistributions(i) = this.linearDistributions(i).multiply(other.linearDistributions(i));
            end
        end
        
        function muLin = linearMean(this)
            arguments
                this (1,1) StateSpaceSubdivisionGaussianDistribution
            end
            % Use GM reduction for the linear part
            muLin = GaussianMixtureDistribution.mixtureParametersToGaussianParamters(...
                [this.linearDistributions.mu],cat(3,this.linearDistributions.C),this.gd.gridValues'/sum(this.gd.gridValues));
        end
        
        function C = linearCovariance(this)
            % Computes covariance of linear dimensions
            %
            % Returns:
            %   C (linD x linD)
            %       covariance matrix
            %
            % Reinterpret as GM to get linear covariance.
            arguments
                this (1,1) StateSpaceSubdivisionGaussianDistribution
            end
            [~,C] = GaussianMixtureDistribution.mixtureParametersToGaussianParamters(...
                [this.linearDistributions.mu],cat(3,this.linearDistributions.C),this.gd.gridValues'/sum(this.gd.gridValues));
        end
        
        function rbd = convolve(~, ~)
            error('Not supported');
        end
        
        function m = mode(this)
            arguments
                this (1,1) StateSpaceSubdivisionGaussianDistribution
            end
            % We have to check all the modes from the linear distributions
            % and multiply them by the weights of the grid points. The
            % highest value is then the max (not interpolated).
            allC = cat(3,this.linearDistributions.C);
            maxfunValsCond = mvnpdf(zeros(1,size(allC,1)),zeros(1,size(allC,1)),allC); % The pdf value at the mean only depends on the covariance matrix
            funValsJoint = maxfunValsCond.*this.gd.gridValues;
            [maxVal, index] = max(funValsJoint);
            funValsJoint(index)=[];
            if any(((maxVal-funValsJoint) < 1e-15)) || any((maxVal./funValsJoint)<1.001)
                warning('Mode:PotentiallyMultimodal','Density may not be unimodal. However, this can also be caused by a high grid resolution and thus very similar function values at the grid points.')
            end
            mPeriodic = this.gd.getGridPoint(index);
            m = [mPeriodic;this.linearDistributions(index).mu];
        end
    end
end

