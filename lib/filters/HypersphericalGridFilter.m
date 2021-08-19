classdef HypersphericalGridFilter < AbstractHypersphericalFilter & AbstractGridFilter
    
    methods
        function this = HypersphericalGridFilter(noOfCoefficients, dim, gridType)
            % Constructor
            %
            % Parameters:
            %   noOfCoefficients (integer > 0)
            %       number of grid values to use
            %   dim (integer > 0)
            %       dimension of the space in which the hypersphere is
            %       embedded
            %	gridType (char)
            %       grid type to use for the HypersphericalGridDistribution
            arguments
                noOfCoefficients {mustBeInteger,mustBePositive}
                dim {mustBeInteger,mustBePositive}
                gridType char = 'eq_point_set'
            end
            % Only allow eq_point_set because large errors are made in the
            % prediction step for the spherical harmonics-based grid.
            this.gd = HypersphericalGridDistribution.fromDistribution(HypersphericalUniformDistribution(dim),...
                noOfCoefficients, gridType);
        end
        
        function setState(this, gd_)
            % Sets the current system state
            %
            % Parameters:
            %   sgd_ (AbstractHypersphericalDistribution)
            %       new state
            arguments
                this HypersphericalGridFilter
                gd_ AbstractHypersphericalDistribution
            end
            assert(this.dim==gd_.dim);
            if ~(isa(gd_, 'HypersphericalGridDistribution'))
                warning('setState:nonGrid', 'sgd_ is not a GridDistribution. Transforming with a number of coefficients that is equal to that of the filter.');
                gd_ = this.gd.fromDistribution(gd_, size(this.gd.gridValues,1), 'eq_point_set'); % Use this to preserve the class (e.g., SphericalGridFilter)
                assert(isequal(this.gd.getGrid(), gd_.getGrid()));
            else
                if ~isequal(this.gd.getGrid(), gd_.getGrid())
                    warning('setState:gridDiffers', 'New density is defined on different grid.')
                end
            end
            this.gd = gd_;
        end
        function predictIdentity(this, dSys)
            % Could transfer to spherical harmonics. Currently only
            % supporting VMFDistribution (predictIdentity only makes sense
            % for zonal densities).
            arguments
                this HypersphericalGridFilter
                dSys VMFDistribution
            end
            warning('PredictIdentity:Inefficient','Using inefficient prediction. Consider precalculating the SdCondSdGridDistribution and using predictNonlinearViaTransitionDensity.');
            trans = @(xkk,xk)cell2mat(arrayfun(@(i)...
                VMFDistribution(xk(:,i),dSys.kappa).pdf(xkk),1:size(xk,2),'UniformOutput', false));
            sdsd = SdCondSdGridDistribution.fromFunction(trans, size(this.gd.gridValues,1), true, 'eq_point_set', 2*this.dim);
            this.predictNonlinearViaTransitionDensity(sdsd);
        end
        
        function updateIdentity(this, measNoise, z)
            % This function was added for interface compatibility with the
            % other filters. measNoise must be a VMF or WatsonDistribution.
            arguments
                this HypersphericalGridFilter
                measNoise AbstractHypersphericalDistribution
                z (:,1) double
            end
            assert(isprop(measNoise,'mu'));
            assert(size(z,1)==this.dim);
            if nargin==3
                if norm(measNoise.mu-[zeros(this.dim-1,1); 1]) > 1E-6
                    error('UpdateIdentity:UnexpectedMeas', 'z needs to be [0;...; 0; 1] to use updateIdentity.');
                end
                measNoise.mu = z;
            end
            currGrid = this.gd.getGrid();
            this.gd = this.gd.multiply(HypersphericalGridDistribution(currGrid,measNoise.pdf(currGrid)'));
        end
        
        function predictNonlinearViaTransitionDensity(this, fTrans)
            arguments
                this HypersphericalGridFilter
                fTrans SdCondSdGridDistribution
            end
            assert(isequal(this.gd.getGrid(),fTrans.getGrid()), 'predictNonlinearViaTransitionDensity:gridDiffers',...
                'fTrans is using an incompatible grid.');
            % Multiplication could be realized via this.getEstimate.gridValues'.*fTrans.gridValues
            % combined with marginalization it would be
            % (domainSize/size(this.grid,2))*sum(this.getEstimate.gridValues'.*fTrans.gridValues,2).
            % Faster with matrix multiplication
            this.gd.gridValues = this.gd.getManifoldSize/size(this.gd.gridValues,1)*fTrans.gridValues*this.gd.gridValues;
            this.gd = this.gd.normalize; % This also enforces a normalization if it is violated
        end
    end
    
end
