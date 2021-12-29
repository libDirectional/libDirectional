classdef HypercylindricalStateSpaceSubdivisionGaussianDistribution < StateSpaceSubdivisionGaussianDistribution & HypercylindricalStateSpaceSubdivisionDistribution
    methods
        function this = HypercylindricalStateSpaceSubdivisionGaussianDistribution(gd_, gaussians)
            % In this class, the grid must be on the hypertorus.
            arguments
                gd_ (1,1) HypertoroidalGridDistribution
                gaussians (:,1) GaussianDistribution
            end
            this@HypercylindricalStateSpaceSubdivisionDistribution(gd_, gaussians);
            this@StateSpaceSubdivisionGaussianDistribution(gd_, gaussians);
        end
        
        function m = mode(this)
            arguments
                this (1,1) HypercylindricalStateSpaceSubdivisionGaussianDistribution
            end
            m = mode@StateSpaceSubdivisionGaussianDistribution(this);
        end
        
        function muLin = linearMean(this)
            arguments
                this (1,1) HypercylindricalStateSpaceSubdivisionGaussianDistribution
            end
            % Use GM reduction for the linear part
            muLin = linearMean@StateSpaceSubdivisionGaussianDistribution(this);
        end
        
        function CLin = linearCovariance(this)
            arguments
                this (1,1) HypercylindricalStateSpaceSubdivisionGaussianDistribution
            end
            % Use GM reduction for the linear part
            CLin = linearCovariance@StateSpaceSubdivisionGaussianDistribution(this);
        end
        
        function mu = hybridMean(this)
            arguments
                this (1,1) HypercylindricalStateSpaceSubdivisionGaussianDistribution
            end
            mu = hybridMean@StateSpaceSubdivisionGaussianDistribution(this);
        end
        
        function m = hybridMoment(this)
            arguments
                this (1,1) HypercylindricalStateSpaceSubdivisionGaussianDistribution
            end
            trigMomComplex = this.gd.trigonometricMoment(1);
            trigMomReal = [real(trigMomComplex)';imag(trigMomComplex)'];
            m = [trigMomReal(:); this.linearMean()];
        end
    end
    
    methods (Static)
        function hcrbd = fromDistribution(dist, noOfGridPoints, gridType)
            arguments
                dist AbstractHypercylindricalDistribution
                noOfGridPoints {mustBeInteger,mustBePositive}
                gridType char = 'CartesianProd'
            end
            % Even for WN, the conditional is a HypercylindricalWN divided
            % by a WN. No analytical formula is known to me
            hcrbdNonGauss = HypercylindricalStateSpaceSubdivisionDistribution.fromFunction(@(x)dist.pdf(x),...
                noOfGridPoints, dist.linD, dist.boundD, gridType);
            % Superclass generates CustomLinearDistributions. Convert.
            linDists = arrayfun(@(dist)GaussianDistribution.fromDistribution(dist),hcrbdNonGauss.linearDistributions);
            hcrbd = HypercylindricalStateSpaceSubdivisionGaussianDistribution(hcrbdNonGauss.gd,linDists);
        end
        
        function hcrbd = fromFunction(fun, noOfGridPoints, dimLin, dimBound, gridType, intRange)
            arguments
                fun (1,1) function_handle
                noOfGridPoints (1,1) {mustBeInteger,mustBePositive}
                dimLin (1,1) {mustBeInteger,mustBePositive}
                dimBound (1,1) {mustBeInteger,mustBePositive} = 1
                gridType char = 'CartesianProd'
                intRange (1,:) double = [-inf,inf]
            end
            hcrbdNonGauss = HypercylindricalStateSpaceSubdivisionDistribution.fromFunction(...
                fun, noOfGridPoints, dimLin, dimBound, gridType, intRange);
            hcrbd = HypercylindricalStateSpaceSubdivisionGaussianDistribution(hcrbdNonGauss.gd,linDists);
        end
    end
end

