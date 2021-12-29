classdef (Abstract) AbstractMixture < AbstractDistribution
    properties
        dists (1,:) cell     % Cell array of hypertoroidal distributions
        w (1,:) double {mustBePositive}         % Weights
    end
    
    methods
        function this = AbstractMixture(dists, w)
            arguments
                dists (1,:) cell {mustBeNonempty}  % Cell array of hypertoroidal distributions
                w (1,:) double {mustBeNonnegative,mustBeNonempty}         % Weights
            end
            % Constructor
            assert(all(size(dists) == size(w)),'Mixture:IncompatibleDimensions','size of dists and w must be equal');
            assert(all(dists{1}.dim==cellfun(@(dist)dist.dim,dists))); % Ensure equal dimensions
            this.dim=dists{1}.dim;
            if any(w==0)
                warning('Pruning elements in mixture with weight zero.');
                dists = dists(w~=0);
                w = w(w~=0);
            end
            this.dists=dists;
            if abs(sum(w)-1)>1e-10
                warning('Weights of mixture do not sum to one.');
                this.w = w/sum(w);
            else
                this.w = w;
            end
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
    methods (Sealed)
        function s = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       n samples on the dim-dimensional hypertorus
            
            % Sample component first, then sample from the chosen component
            d = discretesample(this.w,n);
            if isa(this.dists{1},'SE2BinghamDistribution')
              s = zeros(this.dim+1,n);
            else
              s = zeros(this.dim,n);
            end
            occurrences=histc(d,1:length(this.dists));
            count=1;
            indRange = 1:length(this.dists);
            for i=indRange(occurrences~=0) % Skip all for which zero samples are to be drawn
                s(:,count:count+occurrences(i)-1) = this.dists{i}.sample(occurrences(i));
                count=count+occurrences(i);
            end
            [~,order]=sort(d);
            s(:,order)=s;
        end
    end
end