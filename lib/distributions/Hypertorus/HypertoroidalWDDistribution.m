classdef HypertoroidalWDDistribution < AbstractHypertoroidalDistribution
    % Wrapped Dirac distribution on the hypertorus with dirac positions d and
    % weights w.
    
    properties
        d
        w
    end
    
    methods
        function this = HypertoroidalWDDistribution(d_, w_)
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (dim x L)
            %       Dirac locations in [0,2pi)^dim
            %   w_ (1 x L)
            %       weights for each Dirac            
            this.dim = size(d_,1);
            this.d = mod(d_, 2*pi);
            if (nargin<2)
                %all diracs have equal weights
                this.w = ones(1,size(this.d,2))/size(this.d,2);
            else
                assert(size(w_,1) == 1 ,'Weights must be given as a 1 x L vector.');
                assert(size(d_,2) == size(w_,2),'Number of Dircas and weights must match.');
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
        
        function m = trigonometricMoment(this,n)
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
            m=arrayfun(@(i)sum(exp(1i*n*this.d(i,:)).*this.w),1:size(this.d,1)).';
        end
                
        function hwd = applyFunction(this,f)
            % Apply a function f(x) to each Dirac component and obtain its new position
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2*pi)^dim to [0,2*pi)^dim
            % Returns:
            %   twd (ToroidalWDDistribution)
            %       distribution with new Dirac locations (and same
            %       weights as before)
            assert(isa(f,'function_handle'));
            
            d_ = zeros(size(this.d));
            for i=1:size(this.d,2)
                d_(:,i) = f(this.d(:,i));
            end
            hwd=this;
            hwd.d=d_;
        end
        
        function hwd = reweigh(this, f)
            % Uses a function f(x) to calculate the weight of each Dirac
            % component. The new weight is given by the product of the old 
            % weight and the weight obtained with f. Restores normalization
            % afterwards.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2*pi)^dim to [0, infinity)
            %       (needs to support vectorized inputs, i.e., dim x n matrices)
            % Returns:
            %   twd (ToroidalWDDistribution)
            %       distribution with new weights and same Dirac locations
            assert(isa(f,'function_handle'));
            
            wNew = f(this.d);
            assert(all(wNew >= 0));
            assert(sum(wNew) > 0);
            
            hwd=this;
            hwd.w=wNew.*this.w;
            hwd.w=hwd.w/sum(hwd.w);
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
                       
        function wd = marginalizeTo1D(this,dimension)
            % Get marginal distribution in i-th dimension
            %
            % Parameters:
            %   dimension (scalar)
            %       the marginal in which dimension to calculate,
            %       the other dimensions are marginalized out
            % Returns:
            %   wd (WDDistribution)
            %       marginal distribution (marginals are WD-distributed)
            assert(dimension>= 1 || dimension <= this.dim);
            wd = WDDistribution(this.d(dimension,:), this.w);
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
            assert(all(size(shiftAngles) == [this.dim, 1]));
            
            hd = this;
            hd.d = mod(this.d+repmat(shiftAngles, 1, size(this.d,2)),2*pi);
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
        
        function wd = toWD(this)
            % Convert to a WD distribution (only in 1D case)
            %
            % Returns:
            %   wd (WDDistribution)
            %       WDDistribution with same parameters
            assert(this.dim == 1);
            wd = WDDistribution(this.d, this.w);
        end
        
        function twd = toToroidalWD(this)
            % Convert to a toroidal WD distribution (only in 2D case)
            %
            % Returns:
            %   twd (ToroidalWDDistribution)
            %       ToroidalWDDistribution with same parameters
            assert(this.dim == 2);
            twd = ToroidalWDDistribution(this.d, this.w);
        end        
    end
    
    methods (Sealed)
        function l = logLikelihood(this, samples)
            error('PDF:UNDEFINED', 'not supported');
        end     

        function p = pdf(this, xa)
            % Placeholder, pdf does not exist for wrapped Dirac distributions
            p = 0; %HypertoroidalWDDistribution does not have a proper pdf
            warning('PDF:UNDEFINED', 'pdf is not defined')
        end    
        
        function result = integralNumerical(this, varargin)
            error('PDF:UNDEFINED', 'not supported')
        end        
        
        function m = trigonometricMomentNumerical(this,n)
            % Disable numerical calculation of angular moments since it relies on the pdf
            error('PDF:UNDEFINED', 'not supported');
        end   
        
        function s = sampleMetropolisHastings(this, n)
            % Disable sampling algorithm relying on pdf
            error('PDF:UNDEFINED', 'not supported');
        end        

        function d = squaredDistanceNumerical(this, other)
            error('PDF:UNDEFINED', 'not supported');
        end
        
        function kld = kldNumerical(this, other)
            error('PDF:UNDEFINED', 'not supported');
        end        
    end
end

