classdef FIGFilter < AbstractCircularFilter
    
    properties
        gd
    end
    
    methods
        function this = FIGFilter(noOfCoefficients, transformation, enforcePdfNonnegative)
            % Constructor
            %
            % Parameters:
            %   noOfCoefficients (integer > 0)
            %       number of grid values to use
            if nargin == 1, transformation = 'identity'; end
            if nargin < 3, enforcePdfNonnegative = true; end
            this.gd = FIGDistribution.fromDistribution(CircularUniformDistribution(), noOfCoefficients, transformation);
            this.gd.enforcePdfNonnegative = enforcePdfNonnegative;
        end
        
        function setState(this, gd_)
            % Sets the current system state
            %
            % Parameters:
            %   gd_ (AbstractCircularDistribution)
            %       new state
            assert(isa(gd_, 'AbstractCircularDistribution'));
            if ~(isa(gd_, 'FIGDistribution'))
                warning('setState:nonGrid', 'gd_ is not a FIGDistribution. Transforming with a number of coefficients that is equal to that of the filter.');
                gd_ = FIGDistribution.fromDistribution(gd_, numel(this.gd.gridValues), this.gd.transformation);
            else
                if ~strcmp(this.gd.transformation, gd_.transformation)
                    warning('setState:transDiffer', 'New density is transformed differently.');
                end
                if ~isequal(size(this.gd.gridValues), size(gd_.gridValues))
                    warning('setState:noOfGridValuesDiffer', 'New density has different number of coefficients.')
                end
            end
            this.gd = gd_;
        end
        
        function gd = getEstimate(this)
            % Return current estimate
            %
            % Returns:
            %   gd (FIGDistribution)
            %       current estimate
            gd = this.gd;
        end
        
        function mean = getEstimateMean(this)
            mean = this.gd.circularMean;
        end
        
        function predictIdentity(this, dSys)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by dSys.
            %
            % Parameters:
            %   dSys (AbstractCircularDistribution or samples on grid)
            %       distribution of additive noise or samples on a grid. If
            %       values on a grid are given, take care that they are
            %       not rooted yet. Take care that if
            %       linspace(0,2*pi,noCoeffs) is used, that the last grid
            %       point has to be discarded. For example
            %       xvals=0:2*pi/noCoeffs:2*pi-2*pi/noCoeffs;
            assert(isa(dSys, 'AbstractCircularDistribution'))
            % If a distribution on a circle is given, use Hadamard
            % product approach
            if ~isa(dSys, 'FIGDistribution')
                warning('Predict:automaticConversion', ...
                    'dSys is not a FIGDistribution. Transforming with a number of coefficients that is equal to that of the filter. For non-varying noises, transforming once is much more efficient and should be preferred.');
                dSys = FIGDistribution.fromDistribution(dSys, numel(this.gd.gridValues), this.gd.transformation);
            end
            this.gd = this.gd.convolve(dSys);
        end
        
        function updateIdentity(this, dMeas, z)
            % Updates assuming identity measurement model, i.e.,
            % z(k) = x(k) + v(k)    mod 2pi,
            % where v(k) is additive noise given by dMeas.
            %
            % Parameters:
            %   dMeas (AbstractCircularDistribution)
            %       distribution of additive noise
            %   z (scalar)
            %       measurement in [0, 2pi)
            assert(isa(dMeas, 'AbstractCircularDistribution'));
            if ~isa(dMeas, 'FIGDistribution') % Include shift in transformation
                dMeas = FIGDistribution.fromFunction(@(x)dMeas.pdf(x-z), numel(this.gd.gridValues), this.gd.transformation);
            elseif z > 1e-8 % .shift is rather expensive, avoid if possible
                dMeas = dMeas.shift(z);
            end
            this.gd = this.gd.multiply(dMeas);
        end
        
        function predictNonlinear(this, f, noiseDistribution, truncateJointSqrt)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k)    mod 2pi,
            % where w(k) is additive noise given by noiseDistribution.
            % Using predictNonlinearViaTransitionDensity with a
            % pretransformed fTrans is to be preferred as the function is
            % always performed anew when using predictNonlinear.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi) to [0,2pi).
            %       fuction MUST support row vectors of different angles
            %       and output a row vector!
            %   noiseDistribution (AbstractCircularDistribution)
            %       distribution of additive noise
            
            assert(isa(noiseDistribution, 'AbstractCircularDistribution'));
            assert(isa(f, 'function_handle'));
            if nargin == 3, truncateJointSqrt = true; end
            fTrans = @(xkk, xk)reshape(noiseDistribution.pdf(xkk(:)'-f(xk(:)')), size(xk));
            this.predictNonlinearViaTransitionDensity(fTrans, truncateJointSqrt);
            
        end
        
        function predictNonlinearViaTransitionDensity(this, fTrans, truncateJointSqrt) %#ok<INUSD>
            error('Currently unsupported');
        end
        
        function updateNonlinear(this, likelihood, z) %measurement z, likelihood(z,x)=P(Z|X)
            fdMeas = FIGDistribution.fromFunction( ...
                ... % Assume likelihood can use implicit expansion (for scalars also possible in older Matlab versions)
                @(x)reshape(likelihood(z, x(:)'), size(x)), ...
                numel(this.gd.gridValues), this.gd.transformation);
            this.updateIdentity(fdMeas, zeros(size(z)));
        end
        
        function likelihoodVal = associationLikelihood(this, likelihood)
            assert(numel(this.getEstimate.a) == numel(likelihood.a));
            assert(numel(this.getEstimate.transformation) == numel(likelihood.transformation));
            
            if strcmp(this.getEstimate.transformation, 'identity')
                likelihoodVal = 2 * pi * real(this.getEstimate.c*likelihood.c');
            elseif strcmp(this.getEstimate.transformation, 'sqrt')
                likelihoodVal = 2 * pi * norm(conv(this.getEstimate.c, likelihood.c)).^2;
            else
                error('Transformation not supported');
            end
        end
        
    end
    
end
