classdef CustomSE2Distribution < CustomHypercylindricalDistribution & AbstractSE2Distribution
    methods
        function this = CustomSE2Distribution(mu_, C_)
            arguments
                mu_ (:,1) double
                C_ (:,:) double
            end
            this@CustomHypercylindricalDistribution(mu_, C_, 1);
        end
    end
end