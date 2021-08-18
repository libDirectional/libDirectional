classdef HypertoroidalGridFilter < AbstractHypertoroidalFilter & AbstractGridFilter
    methods
        function this = HypertoroidalGridFilter(noOfCoefficients, dim, gridType)
            % Constructor
            %
            % Parameters:
            %   noOfCoefficients (integer > 0)
            %       number of grid values to use
            %   dim (integer > 0)
            %       dimension of the space in which the hypersphere is
            %       embedded
            %	gridType (char)
            %       grid type to use for the HypertoroidalGridDistribution
            arguments
                noOfCoefficients {mustBeInteger,mustBePositive}
                dim {mustBeInteger,mustBePositive} = numel(noOfCoefficients) % For compatibility with other grid filters
                gridType char = 'CartesianProd'
            end
            % Only allow eq_point_set because large errors are made in the
            % prediction step for the spherical harmonics-based grid.
            this.gd = HypertoroidalGridDistribution.fromDistribution(HypertoroidalUniformDistribution(dim),...
                noOfCoefficients, gridType);
        end
        
        function setState(this, gd_)
            % Sets the current system state
            %
            % Parameters:
            %   sgd_ (AbstractHypertoroidalDistribution)
            %       new state
            arguments
                this (1,1) HypertoroidalGridFilter
                gd_ (1,1) HypertoroidalGridDistribution % Do not transform automatically because size cannot be automatically restored
            end
            assert(this.dim==gd_.dim);
            if ~isequal(this.gd.getGrid(), gd_.getGrid())
                warning('setState:gridDiffers', 'New density is defined on different grid.')
            end
            this.gd = gd_;
        end
        function predictIdentity(this, dSys)
            arguments
                this (1,1) HypertoroidalGridFilter
                dSys (1,1) AbstractHypertoroidalDistribution
            end
            if ~isa(dSys,'HypertoroidalGridDistribution')
                warning('PredictIdentity:automaticConversion','Transforming system noise. Consider pretransformation');
                dSys = HypertoroidalGridDistribution.fromDistribution(dSys, this.gd.noOfGridPoints,'CartesianProd');
            end
            if strcmp(this.gd.gridType,'CartesianProd')
                this.gd = this.gd.convolve(dSys);
            else
                warning('PredictIdentity:Inefficient','Using inefficient prediction because the gridType is not Cartesian product');
            
                trans = @(xkk,xk)cell2mat(arrayfun(@(i)...
                    dSys.shift(xk(:,i)).pdf(xkk),1:size(xk,2),'UniformOutput', false));
                tdtd = TdCondTdGridDistribution.fromFunction(trans, size(this.gd.grid,2), true, 'CartesianProd', 2*this.dim);
                this.predictNonlinearViaTransitionDensity(tdtd);
            end
        end
        
        function updateIdentity(this, measNoise, z)
            % This function was added for interface compatibility with the
            % other filters. measNoise must be a VMF or WatsonDistribution.
            arguments
                this (1,1) HypertoroidalGridFilter
                measNoise (1,1) AbstractHypertoroidalDistribution
                z (:,1) double
            end
            if nargin==3 && norm(z)>0
                measNoise = measNoise.shift(z);
            end
            if ~isa(measNoise,'HypertoroidalGridDistribution')
                % Do not throw warning since it is pretty normal that the
                % likelihood changes. 
                measNoise = HypertoroidalGridDistribution(this.gd.getGrid(),measNoise.pdf(this.gd.getGrid())');
            end
            this.gd = this.gd.multiply(measNoise);
        end
        
        function predictNonlinear(this, a, noiseDistribution)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k)    mod 2pi,
            % where w(k) is additive noise given by noiseDistribution.
            % Using predictNonlinearViaTransitionDensity with a
            % pretransformed fTrans is to be preferred as the function is
            % always performed anew when using predictNonlinear.
            %
            % Parameters:
            %   a (function handle)
            %       function using d input arguments giving d output
            %       arguments. Must support arbitrary dimensional tensors
            %       as input arguments (for vectorized evaluation that is
            %       later used for the transition density fTrans that is build based on f)
            %   noiseDistribution (AbstractHypertoroidalDistribution)
            %       distribution of additive noise
            gdTrans = this.getfTransAsGd(a,noiseDistribution);
            this.predictNonlinearViaTransitionDensity(gdTrans);
        end
        
        function gdTrans=getfTransAsGd(this,a,noiseDistribution)
            % Can call separately to use hfd in multiple time steps with
            % only one transformation
            assert(isa (noiseDistribution, 'AbstractHypertoroidalDistribution'));
            assert(isa(a,'function_handle'));
            fTrans= @(xkk,xk) reshape(noiseDistribution.pdf(xkk-a(xk)),repmat(this.gd.noOfGridPoints,[1,2]));
            gdTrans = TdCondTdGridDistribution.fromFunction(fTrans, this.gd.noOfGridPoints, false, 'CartesianProd', 2*this.dim);
        end
        
        function predictNonlinearViaTransitionDensity(this, fTrans)
            arguments
                this (1,1) HypertoroidalGridFilter
                fTrans (1,1) TdCondTdGridDistribution
            end
            assert(isequal(this.gd.getGrid(),fTrans.getGrid()), 'predictNonlinearViaTransitionDensity:gridDiffers',...
                'fTrans is using an incompatible grid.');
            % Multiplication could be realized via this.getEstimate.gridValues'.*fTrans.gridValues
            % combined with marginalization it would be
            % ((2*pi)^dim / size(this.grid,2))*sum(this.getEstimate.gridValues'.*fTrans.gridValues,2).
            % Faster with matrix multiplication
            this.gd.gridValues = this.gd.getManifoldSize/size(this.gd.gridValues,1)*fTrans.gridValues*this.gd.gridValues;
            this.gd = this.gd.normalize; % This also enforces a normalization if it is violated
        end
    end
    
end
