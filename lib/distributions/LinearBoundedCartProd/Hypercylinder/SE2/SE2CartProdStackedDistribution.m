classdef SE2CartProdStackedDistribution < CartProdStackedDistribution & AbstractSE2Distribution
    methods
        function this = SE2CartProdStackedDistribution(dists)
            arguments
                dists (2,1) cell
            end
            assert(dists{1}.dim==1);
            assert(isa(dists{1},'AbstractHypertoroidalDistribution'));
            assert(dists{2}.dim==2);
            assert(isa(dists{2},'AbstractLinearDistribution'));
            this@CartProdStackedDistribution(dists);
            this.boundD = dists{1}.dim;
            this.linD = dists{2}.dim;
        end
        function stackedDist = shift(this, shiftVals)
            arguments
                this (1,1) SE2CartProdStackedDistribution
                shiftVals (3,1) double
            end
            stackedDist = this;
            stackedDist.dists{1} = stackedDist.dists{1}.shift(shiftVals(1));
            stackedDist.dists{2} = stackedDist.dists{2}.shift(shiftVals(2:3));
        end
    end
end

