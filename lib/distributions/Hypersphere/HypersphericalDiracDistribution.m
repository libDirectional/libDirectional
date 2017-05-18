classdef HypersphericalDiracDistribution < AbstractHypersphericalDistribution
    % Dirac distribution on the hypersphere with dirac positions d and
    % weights w.
    
    properties
        d
        w
    end
    
    methods
        function this = HypersphericalDiracDistribution(d_, w_)
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (dim x L)
            %       Dirac locations on the unit hypersphere
            %   w_ (1 x L)
            %       weights for each Dirac            
            this.dim = size(d_,1);
            
            % check normalization
            assert ( max(abs(sum(d_.^2)) - ones(1, size(d_,2))) < 1E-5);
            this.d = d_;
            
            if (nargin<2)
                %all diracs have equal weights
                this.w = ones(1,size(this.d,2))/size(this.d,2);
            else
                assert(size(w_,1) == 1 );
                assert(size(d_,2) == size(w_,2));
                this.w = w_/sum(w_);
            end
        end      
        
        function p = plot(this, varargin)
            % Create an appropriate plot
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to plot/surf command
            % Returns:
            %   p (scalar)
            %       plot handle
            switch this.dim
                case 2
                    % use polar coordinates, plot angle->pdf(angle)
                    p = stem(atan2(this.d(2,:), this.d(1,:)), this.w, varargin{:});
                case 3
                    % plot points on sphere
                    p = scatter3(this.d(1,:),this.d(2,:),this.d(3,:), this.w*size(this.d,2)*20, varargin{:});
                otherwise
                    error('Plotting for this dimension is currently not supported');
            end                        
        end
        
        function result = integral(this)
            % Calculate integral to check normalization
            % (should always be 1)
            % Returns:
            %   result (scalar)
            %       integral over hypersphere surface of pdf (uses
            %       approximation, not very accurate for higher dimensions)            
            result = sum(this.w);
        end
        
        function integralNumerical(~, ~)
            error('PDF:UNDEFINED', 'not supported');
        end
                
        function hdd = applyFunction(this,f)
            % Apply a function f(x) to each Dirac component and obtain its new position
            %
            % Parameters:
            %   f (function handle)
            %       function from S^dim to S^dim
            % Returns:
            %   hdd (HypersphericalDiracDistribution)
            %       distribution with new Dirac locations (and same
            %       weights as before)
            assert(isa(f,'function_handle'));
            
            d_ = zeros(size(this.d));
            for i=1:size(this.d,2)
                d_(:,i) = f(this.d(:,i));
            end
            % check normalization
            assert(max(abs(sum(d_.^2)) - ones(1, size(d_,2))) < 1E-5);
            
            hdd = this;
            hdd.d = d_;
        end
        
        function hdd = reweigh(this, f)
            % Uses a function f(x) to calculate the weight of each Dirac
            % component. The new weight is given by the product of the old 
            % weight and the weight obtained with f. Restores normalization
            % afterwards.
            %
            % Parameters:
            %   f (function handle)
            %       function from S^dim to [0, infinity)
            %       (needs to support vectorized inputs, i.e., dim x n matrices)
            % Returns:
            %   hdd (HypersphericalDiracDistribution)
            %       distribution with new weights and same Dirac locations
            assert(isa(f,'function_handle'));
            
            wNew = f(this.d);
            assert(all(wNew >= 0));
            assert(sum(wNew) > 0);
            
            hdd = this;
            hdd.w = wNew.*this.w;
            hdd.w = hdd.w/sum(hdd.w);
        end

        function s = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       one sample per column
            assert(isscalar(n));
            ids = discretesample(this.w,n);
            s = this.d(:,ids);
        end
        
        function result = entropy(this)
            % Calculates the entropy analytically 
            %
            % Returns:
            %   result (scalar)
            %       entropy of the distribution
            warning('ENTROPY:DISCRETE','entropy is not defined in a continous sense')
            % return discrete entropy
            % The entropy only depends on the weights!
            result = -sum(this.w.*log(this.w));
        end       
        
        function entropyNumerical(~)
            error('PDF:UNDEFINED', 'not supported');
        end   
        
        function wd = toWD(this)
            % Convert to a WD distribution (only in 2D case)
            %
            % Returns:
            %   wd (WDDistribution)
            %       WDDistribution with same parameters
            assert(this.dim == 2);
            wd = WDDistribution(atan2(this.d(2,:), this.d(1,:)), this.w);
        end   
        
        function sampleMetropolisHastings(~, ~)
            % Disable sampling algorithm relying on pdf
            error('PDF:UNDEFINED', 'not supported');
        end            

        function pdf(~, ~)
            % Placeholder, pdf does not exist for Dirac distributions
            error('PDF:UNDEFINED', 'pdf is not defined')
        end     
    end
end

