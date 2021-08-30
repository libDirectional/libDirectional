classdef CustomDistribution < AbstractDistribution
    % Distribution with custom pdf.
    
    properties
        f function_handle
        % Addtional parameters to prevent having to generate functions
        % based on functions too often (can cause large call stack)
        shiftBy (:,1) double {mustBeNonNan,mustBeFinite}
        scaleBy (1,1) double {mustBeNonNan,mustBeFinite,mustBeNonzero} = 1
    end
    
    methods
        function this = CustomDistribution(f_,dim_)
            arguments
                f_ (1,1) function_handle
                dim_ (1,1) double {mustBeNonempty, mustBeInteger, mustBePositive}
            end
            % Constructor, it is the user's responsibility to ensure that f is a valid
            % density and takes arguments of the same form as
            % .pdf, i.e., it needs to be vectorized. Further, the user has
            % to ensure the density is normalized.
            % 
            % Parameters:
            %   f_ (function handle)
            %       pdf of the distribution
            %   dim_ (scalar)
            %       dimension of the real space in which the hypersphere is
            %       embedded
            this.dim = dim_;
            this.shiftBy = zeros(this.dim,1);
            this.f = f_;
        end
        
        function p = pdf(this, xa)
            arguments
                this (1,1) CustomDistribution
                xa (:,:) double {mustBeNonempty}
            end
            % Evaluate pdf at each column of xa.
            %
            % Parameters:
            %   xa (dim x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            assert(size(xa,1)==this.dim);
            p = this.scaleBy * this.f(xa - this.shiftBy);
            assert(isequal(size(p),[1,size(xa,2)])); % Validate output format of pdf is as expected
        end
        
        function cd = shift(this, shiftVector)
            arguments
                this (1,1) CustomDistribution
                shiftVector (:,1) double
            end
            assert(this.dim == size(shiftVector,1));
            cd = this;
            cd.shiftBy = this.shiftBy + shiftVector;
        end
        
        function cd = normalize(this, verify)
            arguments
                this (1,1) CustomDistribution
                verify (:,1) logical = [] % When not setting manually, it will be determined based on the warning status
            end
            cd = this;
            if isempty(verify)
                lastwarn('') % Rest warning in order to prevent considering previous warnings
            end
            cd.scaleBy = cd.scaleBy/this.integral();
            if isempty(verify) % Set it to check the normalization if and only if the integration halted due reaching the maximum function evaluation
                [~,msg] = lastwarn('');
                verify = contains(msg,'maxFunEvalsPass');
            end    
            if verify && abs(cd.integral()-1)>0.001
                warning('Custom:NotYetNormalized','Density is not yet properly normalized.');
            end
        end
        
    end

end
