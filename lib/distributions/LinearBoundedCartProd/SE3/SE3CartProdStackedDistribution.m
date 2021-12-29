classdef SE3CartProdStackedDistribution < CartProdStackedDistribution & AbstractSE3Distribution
    methods
        function this = SE3CartProdStackedDistribution(dists)
            arguments
                dists (2,1) cell
            end
            assert(dists{1}.dim==4);
            assert(isa(dists{1},'AbstractHyperhemisphericalDistribution'));
            assert(dists{2}.dim==3);
            assert(isa(dists{2},'AbstractLinearDistribution'));
            this@CartProdStackedDistribution(dists);
            this.boundD = dists{1}.dim;
            this.linD = dists{2}.dim;
            this.periodicManifoldType = 'hyperhemisphere';
        end
        function dist = marginalizeLinear(this)
            arguments
                this (1,1) SE3CartProdStackedDistribution
            end
            dist = this.dists{1};
        end
        function dist = marginalizePeriodic(this)
            arguments
                this (1,1) SE3CartProdStackedDistribution
            end
            dist = this.dists{2};
        end
    end
end

