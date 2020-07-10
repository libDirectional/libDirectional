classdef (Abstract) AbstractMixture < AbstractDistribution
    properties
        dists cell     % Cell array of hypertoroidal distributions
        w double {mustBePositive}         % Weights
    end
    
    methods
        function this = AbstractMixture(dists, w)
            arguments
                dists cell     % Cell array of hypertoroidal distributions
                w double {mustBeNonnegative}         % Weights
            end
            % Constructor
            assert(all(size(dists) == size(w)),'size of dists and w must be equal');
            assert(all(dists{1}.dim==cellfun(@(dist)dist.dim,dists))); % Ensure equal dimensions
            this.dim=dists{1}.dim;
            this.dists=dists;
            this.w=w/sum(w);
        end
        
        function p = pdf(this,xa)
            % Evaluate pdf at each column of xa.
            %
            % Parameters:
            %   xa (dim x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location            
            assert(size(xa,1)==this.dim);
            p = zeros(1, size(xa,2));
            for i=1:length(this.dists)
                p = p + this.w(i)*this.dists{i}.pdf(xa); % Calculate pdf using individual pdfs
            end
        end 
    end
end