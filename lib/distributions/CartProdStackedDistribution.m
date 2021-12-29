classdef CartProdStackedDistribution < AbstractCartProdDistribution
    % This is for arbitrary Cartesian products of indepedent (!)
    % distributions. Their states are assumed to be simply stacked
    properties
        dists (:,1) cell
    end
    methods
        function this = CartProdStackedDistribution(dists)
            arguments
                dists (:,1) cell
            end
            this.dists = dists;
            this.dim = sum(cellfun(@(dist)dist.dim,dists));
        end
        
        function s = sample(this, n)
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       one sample per column
            arguments
                this (1,1) CartProdStackedDistribution
                n (1,1) {mustBePositive,mustBeInteger}
            end
            s = cell2mat(cellfun(@(dist){dist.sample(n)},this.dists));
        end
        
        function p = pdf(this, xa)
            ps = NaN(numel(this.dists),size(xa,2));
            nextDim = 1;
            for i=1:numel(this.dists)
                ps(i,:) = this.dists{i}.pdf(xa(nextDim:nextDim+this.dists{i}.dim-1,:));
                nextDim = nextDim+this.dists{i}.dim;
            end
            p = prod(ps,1);
        end

        function stackedDists = shift(this, offsets)
            % This will fail if not all can be shifted
            arguments
                this (1,1) CartProdStackedDistribution
                offsets (:,1) double
            end
            assert(numel(offsets)==this.dim);
            stackedDists = this;
            currDim = 1;
            for i = 1:numel(this.dists)
                stackedDists.dists{i} = this.dists{i}.shift(offsets(currDim:currDim+this.dists{i}.dim-1));
                currDim = currDim + this.dists{i}.dim;
            end
        end

        function dist = setMode(this, newMode)
            arguments
                this (1,1) CartProdStackedDistribution
                newMode (:,1) double
            end
            newDists = this.dists;
            currInd = 1;
            for i = 1:numel(this.dists)
                newDists{i} = newDists{i}.setMode(newMode(currInd:currInd+newDists{i}.dim-1));
                currInd = currInd + newDists{i}.dim;
            end
            dist = this;
            dist.dists = newDists;
        end
    end
    methods (Sealed)
        function m = hybridMean(this)
            arguments
                this (1,1) CartProdStackedDistribution
            end
            meanCell = cellfun(@(d){d.mean},this.dists);
            m = cat(1,meanCell{:});
        end
        function m = mean(this)
            arguments
                this (1,1) CartProdStackedDistribution
            end
            m = this.hybridMean();
        end
        function m = mode(this)
            arguments
                this (1,1) CartProdStackedDistribution
            end
            muCell = cell(numel(this.dists),1);
            for i = 1:numel(this.dists)
                muCell{i} = this.dists{i}.mode();      
            end
            m = cat(1,muCell{:});
        end
    end
end

