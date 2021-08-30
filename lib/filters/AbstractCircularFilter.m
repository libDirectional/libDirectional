classdef (Abstract) AbstractCircularFilter < AbstractHypertoroidalFilter
    % Abstract base class for filters on the circle (S1, SO(2))
    % parameterized as [0, 2pi)
    methods
        function likelihoodVal = associationLikelihood(this, likelihood)
            % Provides the association likelihood
            likelihoodVal = associationLikelihoodNumerical(this, likelihood);
        end
        function likelihoodVal = associationLikelihoodNumerical(this, likelihood)
            likelihoodVal = integral(@(x)reshape( ...
                this.getEstimate.pdf(x(:)').*likelihood.pdf(x(:)'), size(x)), 0, 2*pi);
        end
    end
end
