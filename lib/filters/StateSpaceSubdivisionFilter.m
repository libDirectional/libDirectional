classdef StateSpaceSubdivisionFilter < AbstractHypercylindricalFilter
    properties
        apd StateSpaceSubdivisionGaussianDistribution = StateSpaceSubdivisionGaussianDistribution.empty
        % Initialize empty and let user overwrite it with setState
    end
    
    methods
        function setState(this, apd_)
            arguments
                this (1,1) StateSpaceSubdivisionFilter
                apd_ (1,1) StateSpaceSubdivisionGaussianDistribution
            end
            if ~isempty(this.apd)&&~isequal(numel(this.apd.gaussians),numel(apd_.gaussians))
                warning('LinPeriodic:NoComponentsDiffer','Number of components differ.');
            end
            this.apd = apd_;
        end
        
        function predictLinear(this, transitionDensity, covarianceMatrices, systemMatrices, linearInputVectors)
            arguments
                this (1,1) StateSpaceSubdivisionFilter
                transitionDensity % Transition density empty: Assume dirac (weights form identity matrix). (Do not check class allow it to stay empty.)
                covarianceMatrices (:,:,:,:) double = NaN(0,0,0,0) % If empty: Assume dirac, no covariance added
                systemMatrices (:,:,:,:) double = NaN(0,0,0,0) % If empty: Assume identity matrix
                linearInputVectors (:,1,:,:) double = NaN(0,0,0) % If empty: Assume zero input
            end
            assert(isempty(transitionDensity)...
                ||isa(this.apd.gd,'HyperhemisphericalGridDistribution')&&isa(transitionDensity,'SdHalfCondSdHalfGridDistribution')...
                ||isa(this.apd.gd,'HypersphericalGridDistribution')&&isa(transitionDensity,'SdCondSdGridDistribution')...
                ||isa(this.apd.gd,'HypertoroidalGridDistribution')&&isa(transitionDensity,'TdCondTdGridDistribution'),...
                'Conditional must match grid distribution.');
            noAreas = numel(this.apd.linearDistributions);
            % If matrices are non empty (i.e., zero/identity) the
            % dimensions must match those of the state
            assert(isempty(covarianceMatrices) || size(covarianceMatrices,1) == this.apd.linD && size(covarianceMatrices,2) == this.apd.linD);
            assert(isempty(systemMatrices) || size(systemMatrices,1) == this.apd.linD && size(systemMatrices,2) == this.apd.linD);
            % Must be either empty (then assume zero/identity matrices), a
            % matrix (assume same for all) or have the right dimension
            assert(size(covarianceMatrices,3)<=1 || size(covarianceMatrices,3)==noAreas);
            assert(size(covarianceMatrices,4)<=1 || size(covarianceMatrices,4)==noAreas);
            assert(size(systemMatrices,3)<=1 || size(systemMatrices,3)==noAreas);
            assert(size(systemMatrices,4)<=1 || size(systemMatrices,4)==noAreas);
            % Linear input can be either a vector (same for all) or have
            % noAreas along the 3rd or 4th dimension
            assert(isempty(linearInputVectors)||size(linearInputVectors,1)==this.apd.linD);
            assert(isempty(linearInputVectors)||size(linearInputVectors,2)==1);
            assert(size(linearInputVectors,3)<=1 || size(linearInputVectors,3)==noAreas);
            assert(size(linearInputVectors,4)<=1 || size(linearInputVectors,4)==noAreas);
            
            
            if isempty(transitionDensity) && isempty(covarianceMatrices) && isempty(systemMatrices) && isempty(linearInputVectors)
                %%%%%%%%%%%%%%%%%%%%%%%%% Case 1 %%%%%%%%%%%%%%%%%%%%%%%%%
                % All Diracs and no system is applied: Do nothing
                warning('StateSpaceSubdivisionFilter:NoParamsForPred','Nothing to do for this prediction step.');
                return
            elseif isempty(transitionDensity)
                %%%%%%%%%%%%%%%%%%%%%%%%% Case 2 %%%%%%%%%%%%%%%%%%%%%%%%%
                % Transition density is empty: Transition density is
                % certain mapping with Dirac! Only change the linear parts.
                % In this case, system matrices *must* be on third
                % dimension. It should not need to depend on both the
                % current and next orientation because the circular
                % transition "density" is a Dirac.
                assert(size(covarianceMatrices,4)<=1&&size(systemMatrices,4)<=1&&size(linearInputVectors,4)<=1,...
                    'When we work without any uncertainty in the bounded/periodic domain, the system matrices, covariance matrices, and input vectors for the linear part must be stacked along the third dimension.') 
                
                % Generate indices for which sys/cov mat or input to use. If there is
                % only one, index is always 1, otherwise it is l
                if size(systemMatrices,3)<=1,sysMatIndices = ones(1,noAreas);else,sysMatIndices = 1:noAreas;end
                if size(linearInputVectors,3)<=1,inputIndices = ones(1,noAreas);else,inputIndices = 1:noAreas;end
                if size(covarianceMatrices,3)<=1,covMatIndices = ones(1,noAreas);else,covMatIndices = 1:noAreas;end
                
                if ~isempty(systemMatrices)
                    for i=1:noAreas
                        this.apd.linearDistributions(i).mu = systemMatrices(:,:,sysMatIndices(i))*this.apd.linearDistributions(i).mu;
                        this.apd.linearDistributions(i).C =...
                            systemMatrices(:,:,sysMatIndices(i))*this.apd.linearDistributions(i).C*systemMatrices(:,:,sysMatIndices(i))';
                    end
                end
                if ~isempty(linearInputVectors)
                    for i=1:noAreas
                        this.apd.linearDistributions(i).mu = this.apd.linearDistributions(i).mu + linearInputVectors(:,:,inputIndices(i));
                    end
                end
                if ~isempty(covarianceMatrices)
                    for i=1:noAreas
                        this.apd.linearDistributions(i).C = this.apd.linearDistributions(i).C + covarianceMatrices(:,:,covMatIndices(i));
                    end
                end
                
            else % Now the case in which there is a transition density
                %%%%%%%%%%%%%%%%%%%%%%%%% Case 3 %%%%%%%%%%%%%%%%%%%%%%%%%
                % A transition density for the periodic part is given. This
                % will require updating everything.
                
                % Calculate weights for joint density to get both
                % periodic and linear part
                weightsJoint = transitionDensity.gridValues.*this.apd.gd.gridValues';
                % Determine new periodic part
                this.apd.gd.gridValues = this.apd.gd.getManifoldSize/size(this.apd.gd.gridValues,1)*sum(weightsJoint,2);
                
                % Calculate linear part. To parallelize this, we reshape 
                % xEst into a linD x 1    x 1 x noAreas tensor and
                % Cest into a linD x linD x 1 x noAreas tensor
                xEsts = cat(4,this.apd.linearDistributions.mu);
                CEsts = cat(4,this.apd.linearDistributions.C);
                
                % Apply system matrix, if any is given
                if ~isempty(systemMatrices)
                    xPreds = pagemtimes(systemMatrices,xEsts);
                    CPreds = pagemtimes(pagemtimes(systemMatrices,CEsts),'none',systemMatrices,'transpose');
                else
                    xPreds = xEsts;
                    CPreds = CEsts;
                end
                % Add system noise covariance, if any. Works for single
                % matrices due to implicit expansion
                if ~isempty(covarianceMatrices)
                    CPreds = CPreds + covarianceMatrices;
                end
                % Add input vector, if any. Implicit expansion works here
                % if only one vector is given
                if ~isempty(linearInputVectors)
                    xPreds = xPreds + linearInputVectors;
                end
                
                % Determine all linear distributions
                % The fourth dimension of xPreds and Cpreds should always be occoupied. If
                % the third dimension has not been occupied with the
                % previous operations, then we need to only use the first
                if size(xPreds,3)<=1,xPredsIndex = ones(1,noAreas);else,xPredsIndex = 1:noAreas;end
                if size(CPreds,3)<=1,CPredsIndex = ones(1,noAreas);else,CPredsIndex = 1:noAreas;end
                
                for i=1:noAreas
                    % Do not use reshape and avoid using constructor. We
                    % can overwrite the linearDistributions 
                    [this.apd.linearDistributions(i).mu,this.apd.linearDistributions(i).C]... % Do not use constructor to save computation
                        = GaussianMixtureDistribution.mixtureParametersToGaussianParamters(...
                            reshape(xPreds(:,:,xPredsIndex(i),:),size(xPreds,1),[]),... % Means
                            reshape(CPreds(:,:,CPredsIndex(i),:),size(CPreds,1),size(CPreds,2),[]),.... % Covs
                            weightsJoint(i,:)/sum(weightsJoint(i,:))); % Weights
                end
            end
        end
         
        function update(this, likelihoodPeriodicGrid, likelihoodsLinear)
            arguments
                this (1,1) StateSpaceSubdivisionFilter
                % If likelihoodPeriodicGrid is empty: assume uniform
                % distribution on the bounded domain (i.e., measurement
                % does not contain any information about it)
                likelihoodPeriodicGrid
                % If likelihoodsLinear is empty: assume uniform
                % distribution on the bounded domain (i.e., measurement
                % does not contain any information about it)
                likelihoodsLinear (:,1) GaussianDistribution = GaussianDistribution.empty
            end
            if isa(likelihoodPeriodicGrid,'AbstractDistribution')
                likelihoodPeriodicGrid = likelihoodPeriodicGrid.pdf(this.apd.gd.getGrid())';
            end
            % Either have no grid values or as many as we have grid points
            assert(isempty(likelihoodPeriodicGrid) || isequal(size(likelihoodPeriodicGrid),size(this.apd.gd.gridValues)));
            assert(isempty(likelihoodsLinear) || all([likelihoodsLinear.dim]==this.apd.linD));
            % Either have zero linear distributions (assume uniform), or
            % have 1 (assume same for all regions) or as many as we have
            % regions
            assert(numel(likelihoodsLinear)<=1 || numel(likelihoodsLinear) == size(this.apd.gd.gridValues,1));
            if isempty(likelihoodPeriodicGrid) && isempty(likelihoodsLinear)
                warning('StateSpaceSubdivisionFilter:NoParamsForUpdate','Nothing to do for this update step.');
                return
            end
            if ~isempty(likelihoodPeriodicGrid)
                this.apd.gd.gridValues = this.apd.gd.gridValues.*likelihoodPeriodicGrid;
            end
            if ~isempty(likelihoodsLinear)
                % If there are linear distributions given, also consider
                % their influence on the grid values
                muPreds = [this.apd.linearDistributions.mu];
                muLikelihoods = [likelihoodsLinear.mu];
                % If only one linear distribtion is available, this will use implicit expansion
                covs = cat(3,this.apd.linearDistributions.C)+cat(3,likelihoodsLinear.C);
                % If only one linear distribution is available, this mvnpdf
                % will return the same as if we repeated it for all
                this.apd.gd.gridValues = this.apd.gd.gridValues.*mvnpdf(muPreds',muLikelihoods',covs);
                % Update current linear distributions
                j = 1;
                for i=1:numel(this.apd.linearDistributions)
                    if numel(likelihoodsLinear)>1 % If the linear distribution is identical for all, we should use this for all
                        j = i;
                    end
                    CEstInvCurr = inv(this.apd.linearDistributions(i).C) + inv(likelihoodsLinear(j).C);
                    muEstCurr = CEstInvCurr\(this.apd.linearDistributions(i).C\this.apd.linearDistributions(i).mu...
                        + likelihoodsLinear(j).C\likelihoodsLinear(j).mu);
                    % Can overwrite quantities (using constructor would be
                    % more expensive)
                    this.apd.linearDistributions(i).mu = muEstCurr;
                    this.apd.linearDistributions(i).C = inv(CEstInvCurr);
                end
            end
            this.apd.gd = this.apd.gd.normalize(warnUnnorm=false);
        end
        
        function rgbd = getEstimate(this)
            arguments
                this (1,1) StateSpaceSubdivisionFilter
            end
            rgbd = this.apd;
        end
        
        function mu = getPointEstimate(this)
            arguments
                this (1,1) StateSpaceSubdivisionFilter
            end
            mu = this.apd.hybridMean();
        end
    end
    
end