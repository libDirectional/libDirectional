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
            this.dim=size(d_,1);
            this.d = mod(d_, 2*pi);
            if (nargin<2)
                %all diracs have equal weights
                this.w = ones(1,size(this.d,2))/size(this.d,2);
            else
                assert(size(w_,1) == 1 );
                assert(size(d_,2) == size(w_,2));
                this.w = w_/sum(w_);
            end
        end
        
        function p = pdf(this, xa)
            % Placeholder, pdf does not exist for wrapped Dirac distributions
            p = 0; %HypertoroidalWDDistribution does not have a proper pdf
            warning('PDF:UNDEFINED', 'pdf is not defined')
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
        
        function m = trigonometricMomentNumerical(this,n)
            % Disable numerical calculation of angular moments since it relies on the pdf
            error('PDF:UNDEFINED', 'not supported');
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
               
        function s = sampleMetropolisHastings(this, n)
            % Disable sampling algorithm relying on pdf
            error('PDF:UNDEFINED', 'not supported');
        end
        
        function wd = marginal(this,dimension)
            % Get marginal distribution in i-th dimension
            %
            % Parameters:
            %   dimension (scalar)
            %       the marginal in which dimension to calculate,
            %       the other dimensions are marginalized out
            % Returns:
            %   wd (WDDistribution)
            %       marginal distribution (marginals are WD-distributed)
            assert(dimension>= 0 || dimension <= this.dim);
            wd = WDDistribution(this.d(dimension,:), this.w);
        end
    end
    methods (Sealed)
        function l = logLikelihood(this, samples)
            error('PDF:UNDEFINED', 'not supported');
        end     
    end
end

