classdef CustomHypertoroidalDistribution < AbstractHypertoroidalDistribution & CustomDistribution
    % Hypertoroidal distribution with custom pdf.
    
    methods
        function this = CustomHypertoroidalDistribution(f_,dim_)
            % Constructor, it is the user's responsibility to ensure that f is a valid
            % hypertoroidal density and takes arguments of the same form as
            % .pdf, i.e., it needs to be vectorized.
            % 
            % Parameters:
            %   f_ (function handle)
            %       pdf of the distribution
            %   dim_ (scalar)
            %       dimension of the hypertorus
            this@CustomDistribution(f_,dim_);
        end

        function ccd = toCustomCircular(this)
            % Convert to a custom circular distribution (only in 1D case)
            %
            % Returns:
            %   ccd (CustomCircularDistribution)
            %       CustomCircularDistribution with same parameters
            assert(this.dim == 1);
            ccd = CustomCircularDistribution(this.f);
        end
        
        function ctd = toCustomToroidal(this)
            % Convert to a custom toroidal distribution (only in 2D case)
            %
            % Returns:
            %   ctd (CustomToroidalDistribution)
            %       CustomToroidalDistribution with same parameters
            assert(this.dim == 2);
            ctd = CustomToroidalDistribution(this.f);
        end
        
        function hd = shift(this, shiftAngles)
            % Shifts the distribution by shiftAngles. If not overloaded,
            % this will return a CustomHypertoroidalDistribution.
            %
            % Parameters:
            %   shiftAngles (dim x 1 column vector)
            %       angles to shift by
            % Return:
            %   hd (CustomHypertoroidalDistribution)
            %       shifted distribution
            arguments
                this (1,1) CustomHypertoroidalDistribution
                shiftAngles (:,1) double
            end
            hd = shift@CustomDistribution(this, shiftAngles);
        end
        
        function chd = scale(this, scalingFactor)
            arguments
                this (1,1) CustomHypertoroidalDistribution
                scalingFactor (1,1) double {mustBeNonzero,mustBeNonNan,mustBeFinite}
            end
            if scalingFactor<0
                warning('CustomHypertoroidalDistribution:NegativeScaling','Scaling factor is smaller than 0, check if this is intended.');
            end
            chd = this;
            chd.scaleBy = chd.scaleBy*scalingFactor;
        end
        
        function chdMarginalized = marginalizeOut(this, dimensions)
            arguments
                this (1,1) CustomHypertoroidalDistribution
                dimensions (1,:) double {mustBeInteger,mustBePositive}
            end
            assert(all(dimensions<=this.dim),'Cannot perform this operation for a dimension that is higher than the dimensionality of the distribution.');
            remainingDims = 1:this.dim;
            remainingDims(dimensions) = [];
            chdMarginalized = CustomHypertoroidalDistribution(...
                @(xs)arrayfun(@(i)this.sliceAt(remainingDims,xs(:,i)).integral(),1:size(xs,2)),numel(remainingDims));
        end
        
        function chdCond = conditionOn(this, dims, val)
            arguments
                this (1,1) CustomHypertoroidalDistribution
                dims (1,:) double {mustBeInteger,mustBePositive}
                val (:,1) double
            end
            assert(all(dims<=this.dim),'Cannot perform this operation for a dimension that is higher than the dimensionality of the distribution.');
            chdSliced = this.sliceAt(dims, val);
            chdCond = chdSliced.normalize();
        end
        
        function chdLikelihood = likelihood(this, dims, val)
            arguments
                this (1,1) CustomHypertoroidalDistribution
                dims (1,:) double {mustBeInteger,mustBePositive}
                val (:,1) double
            end
            chdSliced = this.sliceAt(dims, val);
            chdMarginalized = this.marginalizeOut(dims);
            chdLikelihood = CustomHypertoroidalDistribution(@(x)...
                chdSliced.pdf(x)./chdMarginalized.pdf(x),this.dim-numel(dims));
        end
        
        function chdSliced = sliceAt(this, dims, val)
            arguments
                this (1,1) CustomHypertoroidalDistribution
                dims (1,:) double {mustBeInteger,mustBePositive}
                val (:,1) double
            end
            assert(all(dims<=this.dim),'Cannot perform this operation on a dimension that is higher than the dimensionality of the distribution.');
            assert(all(diff(dims)>0),'Dimensions must be given in ascending order.');
            dimsRemaining = setdiff(1:this.dim,dims);
            chdSliced = CustomHypertoroidalDistribution(@newPdf,this.dim-numel(dims));
            % We directly constructed it based on the function. Now we need
            % to transfer the shifting
            dimsRemaining = setdiff(1:this.dim,dims);
            chdSliced.shiftBy = this.shiftBy(dimsRemaining);
            chdSliced.scaleBy = this.scaleBy;
            function p = newPdf(xs)
                inputMat = NaN(this.dim,size(xs,2));
                % Have to incorporate the shifting because we directly take
                % the function
                inputMat(dims,:) = repmat(val-this.shiftBy(dims),[1,size(xs,2)]);
                inputMat(dimsRemaining,:) = xs;
                p = this.f(inputMat);
            end
        end
    end
    
    methods (Static)
        function chd = fromDistribution(dist)
            % Creates a CustomHypertoroidalDistribution from some other distribution
            %
            % Parameters:
            %   dist (AbstractHypertoroidalDistribution)
            %       distribution to convert
            % Returns:
            %   chd (CustomHypertoroidalDistribution)
            %       CustomHypertoroidalDistribution with identical pdf
            arguments
                dist (1,1) AbstractHypertoroidalDistribution
            end
            
            chd = CustomHypertoroidalDistribution(@(xa)dist.pdf(xa),dist.dim);
        end
    end
end
