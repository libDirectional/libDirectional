classdef HypertoroidalWDDistribution < AbstractHypertoroidalDistribution & AbstractDiracDistribution
    % Wrapped Dirac distribution on the hypertorus with dirac positions d and
    % weights w.
    
    methods
        function this = HypertoroidalWDDistribution(d_, w_)
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (dim x L)
            %       Dirac locations in [0,2pi)^dim
            %   w_ (1 x L)
            %       weights for each Dirac
            arguments
                d_ (:,:) double
                w_ (1,:) double = ones(1,size(d_,2))/size(d_,2);
            end
            this@AbstractDiracDistribution(mod(d_, 2*pi), w_);
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
            arguments
                this (1,1) HypertoroidalWDDistribution
            end
            arguments (Repeating)
                varargin
            end
            switch this.dim
                case 1
                    p = stem(this.d, this.w, varargin{:});
                case 2
                    p = stem3(this.d(1,:),this.d(2,:), this.w, varargin{:});
                case 3
                    [X,Y,Z] = sphere(4);
                    clf 
                    hold on
                    color = jet;
                    nColors = size(color, 1);
                    scale = 10;
                    maxW = max(this.w);
                    arrayfun(@(x,y,z,currSize) surf(scale*currSize*X+x,scale*currSize*Y+y,scale*currSize*Z+z, 'facecolor', color(1+floor((nColors-1)*currSize/maxW),:)),...
                        this.d(1,:), this.d(2,:), this.d(3,:), this.w);
                    hold off
                    setupAxisCircular('x','y','z')
                    ylabel('x_2')
                    zlabel('x_3')
                    view(40,20)
                    grid                    
                otherwise
                    error('Plotting for this dimension is currently not supported');
            end                        
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
            arguments
                this (1,1) HypertoroidalWDDistribution
                n (1,1) double {mustBeInteger,mustBePositive}
            end
            s = sample@AbstractDiracDistribution(this, n);
        end
        
        function dist = applyFunction(this,f)
            % Regular handling + mod 2*pi
            arguments
                this (1,1) HypertoroidalWDDistribution
                f (1,1) function_handle
            end
            dist = applyFunction@AbstractDiracDistribution(this,f);
            dist.d = mod(dist.d, 2*pi);
        end
        
        function m = trigonometricMoment(this,n)
            arguments
                this (1,1) HypertoroidalWDDistribution
                n (1,1) double {mustBeInteger}
            end
            % Calculate n-th trigonometric moment, i.e., 
            % E([e^(inx_1); e^(inx_2); ...; e^(inx_dim)])
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (dim x 1)
            %       n-th trigonometric moment (complex vector)  
            assert(isscalar(n));
            m = sum(exp(1i*n*this.d).*repmat(this.w,this.dim,1),2);
        end
                       
        function wd = marginalizeTo1D(this,dimension)
            arguments
                this (1,1) HypertoroidalWDDistribution
                dimension (1,1) double {mustBeInteger,mustBePositive}
            end
            % Get marginal distribution in i-th dimension
            %
            % Parameters:
            %   dimension (scalar)
            %       the marginal in which dimension to calculate,
            %       the other dimensions are marginalized out
            % Returns:
            %   wd (WDDistribution)
            %       marginal distribution (marginals are WD-distributed)
            wd = WDDistribution(this.d(dimension,:), this.w);
        end
        
        function wd = marginalizeOut(this, dimensions)
            arguments
                this (1,1) HypertoroidalWDDistribution
                dimensions (1,:) double {mustBeInteger,mustBePositive}
            end
            remainingDims = 1:this.dim;
            remainingDims(dimensions) = [];
            wd = WDDistribution(this.d(remainingDims,:), this.w);
        end
        
        function hd = shift(this, shiftAngles)
            % Shift distribution by the given angles
            %
            % Parameters:
            %   shiftAngles (dim x 1 column vector) 
            %       angles to shift by
            % Returns:
            %   hd (HypertoroidalWDDistribution)
            %       shifted distribution
            arguments
                this (1,1) HypertoroidalWDDistribution
                shiftAngles (:,1) double
            end
            assert(size(shiftAngles,1) == this.dim);
            
            hd = this;
            hd.d = mod(this.d+shiftAngles,2*pi); % Use implicit expansion of +
        end        
        
        function result = entropy(this)
            % Calculates the entropy analytically 
            %
            % Returns:
            %   result (scalar)
            %       entropy of the distribution
            arguments
                this (1,1) HypertoroidalWDDistribution
            end
            warning('ENTROPY:DISCRETE','entropy is not defined in a continous sense')
            % return discrete entropy
            % The entropy only depends on the weights!
            result = -sum(this.w.*log(this.w));
        end        
        
        function wd = toWD(this)
            % Convert to a WD distribution (only in 1D case)
            %
            % Returns:
            %   wd (WDDistribution)
            %       WDDistribution with same parameters
            arguments
                this (1,1) HypertoroidalWDDistribution
            end
            assert(this.dim == 1);
            wd = WDDistribution(this.d, this.w);
        end
        
        function twd = toToroidalWD(this)
            % Convert to a toroidal WD distribution (only in 2D case)
            %
            % Returns:
            %   twd (ToroidalWDDistribution)
            %       ToroidalWDDistribution with same parameters
            arguments
                this (1,1) HypertoroidalWDDistribution
            end
            assert(this.dim == 2);
            twd = ToroidalWDDistribution(this.d, this.w);
        end        

        function result = integral(this)
            arguments
                this (1,1) HypertoroidalWDDistribution
            end
            result = integral@AbstractDiracDistribution(this);
        end
    end
    
    methods(Static)
        function f = fromDistribution(distribution, noOfSamples)
            arguments
                distribution (1,1) AbstractHypertoroidalDistribution
                noOfSamples (1,1) {mustBePositive, mustBeInteger}
            end
            f = HypertoroidalWDDistribution(...
                distribution.sample(noOfSamples),ones(1,noOfSamples)/noOfSamples);
        end
    end
    
end

