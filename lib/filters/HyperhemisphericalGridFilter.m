classdef HyperhemisphericalGridFilter < AbstractGridFilter & AbstractHyperhemisphericalFilter % Could inherit from AbstractAxialFilter, but would be of no use and might be confusing.
    
    methods
        function this = HyperhemisphericalGridFilter(noOfCoefficients, dim, gridType)
            % Constructor
            %
            % Parameters:
            %   noOfCoefficients (integer > 0)
            %       number of grid values to use
            %   dim (integer > 0)
            %       dimension of the space in which the hypersphere is
            %       embedded
            arguments
                noOfCoefficients (1,1) double {mustBeInteger,mustBePositive}
                dim (1,1) double {mustBeInteger,mustBePositive}
                gridType char = 'eq_point_set_symm'
            end
            % Only allow eq_point_set because large errors are made in the
            % prediction step for the spherical harmonics-based grid.
            this.gd = HyperhemisphericalGridDistribution.fromDistribution(HyperhemisphericalUniformDistribution(dim),...
                noOfCoefficients, gridType);
        end
        
        function setState(this, gd_)
            % Sets the current system state
            %
            % Parameters:
            %   gd_ (AbstractDistribution)
            %       new state
            arguments
                this (1,1) HyperhemisphericalGridFilter
                gd_ (1,1) AbstractDistribution
            end
            assert(this.dim==gd_.dim);
            if isa(gd_, 'HyperhemisphericalGridDistribution')
                if ~isequal(this.gd.getGrid(), gd_.getGrid())
                    warning('setState:gridDiffers', 'New density is defined on different grid.')
                end
            elseif isa(gd_, 'HypersphericalGridDistribution')
                warning('Called setState with GridDistribution on the entire hypersphere. Please ensure it is at least symmetric.')
                HyperhemisphericalGridDistribution(this.gd.grid(1:size(this.gd.grid)/2), this.gd.gridValues(1:size(this.gd.grid)/2));
                if ~isequal(this.gd.getGrid(), gd_.getGrid())
                    warning('setState:gridDiffers', 'New density is defined on different grid.')
                end
            else
                warning('setState:nonGrid', 'sgd_ is not a HyperhemisphericalGridDistribution. Transforming with a number of coefficients that is equal to that of the filter.');
                gd_ = HyperhemisphericalGridDistribution.fromDistribution(gd_, size(this.gd.gridValues,1), 'eq_point_set');
                assert(isequal(this.gd.getGrid(), gd_.getGrid()));
            end
            this.gd = gd_;
        end
        
        function predictIdentity(this, dSys)
            % Could transfer to spherical harmonics. Currently only
            % supporting VMFDistribution (predictIdentity only makes sense
            % for zonal densities).
            arguments
                this (1,1) HyperhemisphericalGridFilter
                dSys (1,1) AbstractDistribution
            end
            assert(dSys.dim==this.dim);
            warning('PredictIdentity:Inefficient','Using inefficient prediction. Consider precalculating the SdHalfCondSdHalfGridDistribution and using predictNonlinearViaTransitionDensity.');
            
            sdHalfCondSdHalf = HyperhemisphericalGridFilter.sysNoiseToTransitionDensity(dSys, size(this.gd.gridValues,1));
            this.predictNonlinearViaTransitionDensity(sdHalfCondSdHalf);
        end
        
        function updateIdentity(this, measNoise, z)
            % This function was added for interface compatibility with the
            % other filters. measNoise must be a VMF or WatsonDistribution.
            arguments
                this (1,1) HyperhemisphericalGridFilter
                measNoise (1,1) AbstractDistribution
                z (:,1) double
            end
            assert(size(z,1)==this.dim);
            assert( isa(measNoise,'AbstractHyperhemisphericalDistribution')||...
                isa(measNoise,'WatsonDistribution') ||... % Watson is symmetric
                isa(measNoise,'HypersphericalMixtureDistribution')... % Allow for mixtures of two von mises
                && numel(measNoise.dists)==2 && all(measNoise.w==0.5) && isequal(measNoise{1}.mu,-measNoise{2}.mu)... 
                || isa(measNoise,'VMFDistribution') && z(end)==0); % Allow VMF with means that cause it to be symmetric
            if nargin==3
                if isa(measNoise,'HyperhemisphericalWatsonDistribution')
                    measNoise = measNoise.setMode(z);
                elseif isa(measNoise,'WatsonDistribution') || isa(measNoise,'VMFDistribution') && z(end)==0
                    if norm(measNoise.mu-[zeros(this.dim-1,1); 1]) > 1E-6
                        error('UpdateIdentity:UnexpectedMeas', 'z needs to be [0;...; 0; 1] to use updateIdentity.');
                    end
                    measNoise.mu = z;
                elseif isa(measNoise,'HypersphericalMixtureDistribution')
                    measNoise.dists{1} = z;
                    measNoise.dists{2} = -z;
                else
                    error('Unsupported distribution');
                end
            end
            currGrid = this.gd.getGrid();
            this.gd = this.gd.multiply(HyperhemisphericalGridDistribution(currGrid,2*measNoise.pdf(currGrid)'));
        end
        
        function updateNonlinear(this, likelihood, z) %measurement z, likelihood(z,x)=P(Z|X)
            this.gd.gridValues = this.gd.gridValues .* likelihood(z, this.gd.getGrid())'; % Evaluate on precisely same grid (.multiply would require GridDistribution)
            warnStruct = warning('off', 'Normalization:notNormalized');
            this.gd = this.gd.normalize;
            warning(warnStruct);
        end
        
        function predictNonlinearViaTransitionDensity(this, fTrans)
            arguments
                this (1,1) HyperhemisphericalGridFilter
                fTrans (1,1) SdHalfCondSdHalfGridDistribution
            end
            assert(isequal(this.gd.getGrid(),fTrans.getGrid()), 'predictNonlinearViaTransitionDensity:gridDiffers',...
                'fTrans is using an incompatible grid.');
            % Multiplication could be realized via this.getEstimate.gridValues'.*fTrans.gridValues
            % combined with marginalization it would be
            % (4*pi/size(this.grid,2)*sum(this.getEstimate.gridValues'.*fTrans.gridValues,2).
            % Faster with matrix multiplication
            this.gd = this.gd.normalize;
            gridValuesnew = this.gd.getManifoldSize/size(this.gd.gridValues,1)*fTrans.gridValues*this.gd.gridValues;
            this.gd = HyperhemisphericalGridDistribution(this.gd.getGrid(),gridValuesnew); % This also enforces a normalization if it is violated
        end
        
        function p = getPointEstimate(this)
            arguments
                this (1,1) HyperhemisphericalGridFilter
            end
            % For estimate: fit Bingham to it and obtain mode
            gdFullSphere = this.gd.toFullSphere();
            p = BinghamDistribution.fit(gdFullSphere.getGrid(), gdFullSphere.gridValues'/sum(gdFullSphere.gridValues)).mode();
            if p(end)<0
                p = -p;
            end
        end
    end
    methods (Static)
        function sdHalfCondSdHalf = sysNoiseToTransitionDensity(dSys, noGridPoints)
            arguments
                dSys (1,1) AbstractDistribution
                noGridPoints (1,1) double
            end
            if isa(dSys,'HyperhemisphericalWatsonDistribution')||isa(dSys,'WatsonDistribution')
                trans = @(xkk,xk)cell2mat(arrayfun(@(i)...
                    2*WatsonDistribution(xk(:,i),dSys.kappa).pdf(xkk),1:size(xk,2),'UniformOutput', false));
            elseif isa(dSys,'HypersphericalMixture') && numel(dSys.dists)==2 && all(dSys.w==0.5) && isequal(dSys.dists{1}.mu, -dSys.dists{2}.mu) && dSys.dists{1}.kappa==dSys.dists{2}.kappa
                trans = @(xkk,xk)cell2mat(arrayfun(@(i)... % * 0.5 and * 2 cancel out
                    VMFDistribution(xk(:,i),dSys.dists{1}.kappa).pdf(xkk)+VMFDistribution(xk(:,i),dSys.dists{1}.kappa).pdf(-xkk),1:size(xk,2),'UniformOutput', false));
            else
                error('Distribution not supported for predict identity. Must be zonal (rotationally symmetric around last dimension');
            end
            warning('PredictIdentity:Inefficient','Using inefficient prediction. Consider precalculating the SdHalfCondSdHalfGridDistribution and using predictNonlinearViaTransitionDensity.');
            sdHalfCondSdHalf = SdHalfCondSdHalfGridDistribution.fromFunction(trans, noGridPoints, true, 'eq_point_set_symm', 2*dSys.dim);
        end
    end
end
