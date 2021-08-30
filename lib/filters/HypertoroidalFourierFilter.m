classdef HypertoroidalFourierFilter < AbstractHypertoroidalFilter
    % Filter based on Fourier series on the hypertorus
    %
    % Florian Pfaff, Gerhard Kurz, Uwe D. Hanebeck,
    % Multivariate Angular Filtering Using Fourier Series
    % Journal of Advances in Information Fusion, December 2016.
    properties
        hfd HypertoroidalFourierDistribution
    end
    
    methods
        function this = HypertoroidalFourierFilter(noOfCoefficients,transformation)
            % Constructor. Pay attention that the length of
            % noOfCoefficients determines the dimensionality of the
            % underlying hypertoroidal Fourier distribution
            if nargin==1,transformation='sqrt';end
            this.hfd=HypertoroidalFourierDistribution.fromDistribution(...
                HypertoroidalUniformDistribution(numel(noOfCoefficients)),noOfCoefficients,transformation);
        end
        
        function setState(this, hfd_)
            % Sets the current system state
            %
            % Parameters:
            %   hfd_ (HypertoroidalFourierDistribution)
            %       new state
            assert(isa(hfd_,'AbstractHypertoroidalDistribution'));
            if ~(isa(hfd_,'HypertoroidalFourierDistribution'))
                warning('setState:nonFourier','hfd_ is not a HypertoroidalFourierDistribution. Transforming with a number of coefficients that is equal to that of the filter.');
                sizeHfdC=size(this.hfd.C); % Needed for workaround for 1D case
                hfd_=HypertoroidalFourierDistribution.fromDistribution(hfd_,sizeHfdC(sizeHfdC>1),this.hfd.transformation);
            elseif ~strcmp(this.hfd.transformation,hfd_.transformation)
                warning('setState:transDiffer','New density is transformed differently.');
            end
            if ndims(this.hfd.C)~=ndims(hfd_.C)
                warning('setState:noOfDimsDiffer','New desity has different dimensionality.');
            elseif ~isequal(size(this.hfd.C),size(hfd_.C))
                warning('setState:noOfCoeffsDiffer','New density has different number of coefficients.');
            end
            this.hfd = hfd_;
        end

        function hfd = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   hfd (HypertoroidalFourierDistribution)
            %       current estimate
            hfd=this.hfd;
        end
        
        function mean=getEstimateMean(this)
            mean=this.hfd.circularMean;
        end

        function predictIdentity(this, dSys)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by dSys.
            % The modulo operation is carried out componentwise.
            %
            % Parameters:
            %   dSys (HypertoroidalFourierDistribution)
            %       distribution of additive noise
            sizeHfdC=size(this.hfd.C); % Needed for workaround for 1D case
            if ~(isa(dSys,'HypertoroidalFourierDistribution'))
                warning('PredictIdentity:automaticConversion',...
                    'dSys is not a HypertoroidalFourierDistribution. Transforming with a number of coefficients that is equal to that of the filter. For non-varying noises, transforming once is much more efficient and should be preferred.');
                dSys=HypertoroidalFourierDistribution.fromDistribution(dSys,sizeHfdC(sizeHfdC>1),this.hfd.transformation);
            end
            this.hfd=this.hfd.convolve(dSys,sizeHfdC(sizeHfdC>1));
        end

        function updateIdentity(this, dMeas, z)
            % Updates assuming identity measurement model, i.e.,
            % z(k) = x(k) + v(k)    mod 2pi,
            % where v(k) is additive noise given by dMeas.
            % The modulo operation is carried out componentwise.
            %
            % Parameters:
            %   dMeas (AbstractHypertoroidalDistribution)
            %       distribution of additive noise
            %   z (dim x 1 vector)
            %       measurement in [0, 2pi)^dim
            assert(isa(dMeas,'AbstractHypertoroidalDistribution'));
            if ~(isa(dMeas,'HypertoroidalFourierDistribution'))
                warning('Update:automaticConversion',...
                    'dMeas is not a HypertoroidalFourierDistribution. Transforming with an amount of coefficients that is equal to that of the filter. For non-varying noises, transforming once is much more efficient and should be perferred.');
                sizeHfdC=size(this.hfd.C); % Needed for workaround for 1D case
                dMeas=HypertoroidalFourierDistribution.fromDistribution(dMeas,sizeHfdC(sizeHfdC>1),this.hfd.transformation);
            end
            assert(isequal(size(z),[this.hfd.dim,1]));
            
            dMeasShifted=dMeas.shift(z);
            this.hfd=this.hfd.multiply(dMeasShifted,size(this.hfd.C));
        end
        
        function predictNonlinear(this, a, noiseDistribution, truncateJointSqrt)
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
            if nargin==3,truncateJointSqrt=true;end
            % Input checks are done in getfTrans
            hfdTrans = this.getfTransAsHfd(a,noiseDistribution);
            
            this.predictNonlinearViaTransitionDensity(hfdTrans,truncateJointSqrt);
        end
        
        function hfdTrans=getfTransAsHfd(this,a,noiseDistribution)
            % Can call separately to use hfd in multiple time steps with
            % only one transformation
            assert(isa (noiseDistribution, 'AbstractHypertoroidalDistribution'));
            assert(isa(a,'function_handle'));
            dimC=size(this.hfd.C);
            hfdTrans=HypertoroidalFourierDistribution.fromFunction(@fTrans,dimC([1:this.hfd.dim,1:this.hfd.dim]),this.hfd.transformation);
            function p=fTrans(varargin) % Can be written as lambda function but is not faster and less legible
                assert(all(diff(cellfun(@numel,varargin))==0),...
                    'All input arguments need to be equally sized. Use ndgrid to generate appropriate input arguments.'); 
                fout=cell(1,this.hfd.dim);
                % Second d dimenisons are for x(k), propagate through a.
                [fout{:}]=a(varargin{this.hfd.dim+1:2*this.hfd.dim});
                % Concatenate on the last dimension to ensure that it is iterated
                % through LAST when using reshape. Make sure to make the
                % concatenation dependent on ndims(varargin{1}) because the dimensionality
                % of the input may be higher than the dimensionality of the
                % distribution!
                ws=cat(ndims(varargin{1}),varargin{1:this.hfd.dim})-cat(ndims(varargin{1}),fout{:});
                % Flatten tensor to dim x (numel(varargin{1}) tensor for
                % evaluation in pdf of the noise. Because the values for
                % the different dimensions are iterated through last, use
                % reshape into (numel(varargin{1}) x dim first and
                % transpose.
                wsReshaped=reshape(ws,[],noiseDistribution.dim)';
                pdfvals=noiseDistribution.pdf(wsReshaped);
                % Since this is conditional, integrating over all
                % dimensions will yield (2*pi)^this.dim instead of 1.
                % Also restore to input shape.
                p=reshape(pdfvals,size(varargin{1}))/((2*pi)^this.dim);
            end
        end
        
        function predictNonlinearViaTransitionDensity(this, fTrans,truncateJointSqrt)
            % Predicts assuming a nonlinear system model using a
            % probabilistic model. See also FourierFilter that
            % implements this for for the 1-D case in a real Fourier
            % series representation.
            %
            % Parameters:
            %   fTrans
            %       transition density f(x(k+1)|x(k))
            %       Either given as HypertoroidalFourierDistribution or as
            %       function.
            %       If fTrans is a HypertoroidalFourierDistribution: 
            %       First d dimension are for x(k+1), next d dimensions for
            %       x(k).
            %       If fTrans is a function: It is a function of 2d input
            %       arguments. The first d input arguments are for x(k+1)
            %       and the second d input arguments are for x(k).
            %
            % To reuse code, fTrans is treated as a distribution and is 
            % thus normalized to one although it should integrate to
            % (2*pi)^d. This is because fTrans is not a 2d-variate density but a
            % conditional density that, when integrated over x_k+1
            % always yields 1. By also integrating over x_k, we get
            % (2*pi)^d when d is the dimension of x_k.

            if nargin==2,truncateJointSqrt=true;end
            
            % As explained above, fTrans acutally integrates to (2*pi)^d.
            % This would result in a warning due to the lack of
            % normalization. Warning is therefore disabled.
            dimC=size(this.hfd.C);
            warnStruct=warning('off','Normalization:notNormalized');
            if isa(fTrans,'function_handle')
                assert(nargin(fTrans)==(2*this.hfd.dim)||nargin(fTrans)==-1); % -1 if varargin
                fTrans=HypertoroidalFourierDistribution.fromFunction(fTrans,dimC([1:this.hfd.dim,1:this.hfd.dim]),this.hfd.transformation);
            else
                assert(isa(fTrans,'HypertoroidalFourierDistribution')&&strcmp(this.hfd.transformation,fTrans.transformation)&&fTrans.dim==(2*this.hfd.dim),'predictNonlinear:fTransInvalid',...
                    'fTrans must be given as a function or a HypertoroidalFourierDistribution of appropriate dimension and transformation.');
            end
            
            if strcmp(this.hfd.transformation,'identity')
                % Multiply and marginalize in one step. See test case in 
                % FourierFilter test that illustrates the validity of this 
                % approach. On top of the (2*pi)^dim from the marginalization, 
                % we also have to respect that we have fTrans/((2*pi)^dim).
                % Thus, we multiply by (2*pi)^(2*dim).
                % Use reshape because it is more efficient than permute.
                % Using conv and valid, we could calculate this according
                % to
                % CPredictedId=(2*pi)^(2*this.hfd.dim)*convnc(fTrans.C,...
                %    reshape(this.hfd.C,[ones(1,this.hfd.dim),size(this.hfd.C)]),'valid');
                if this.hfd.dim~=1
                    CPredictedId=(2*pi)^(2*this.hfd.dim)*convnc(fTrans.C,...
                    reshape(this.hfd.C,[ones(1,this.hfd.dim),size(this.hfd.C)]),'valid');
                    % We could reformulate this as, e.g., in 2D as 
                    % sum(sum(a.*b(:,:,end:-1:1,end:-1:1),3),4). But in
                    % some testing, this was not deemed more efficient. In
                else
                    % In 1-D, we can further reformulate it into an
                    % significantly more efficient form.
                    CPredictedId=(2*pi)^(2*this.hfd.dim)*fTrans.C*flip(this.hfd.C);
                end
            elseif strcmp(this.hfd.transformation,'sqrt')
                % The rational is analogous to the one in the Fourier
                % filter. We can truncate and save one convolution if we
                % first calculate CJointSqrt
                if ~truncateJointSqrt
                    % Since we are still in square root representation, we can
                    % scale fTrans with sqrt((2*pi)^dim) to ensure normalization
                    CJointSqrt=convnc(...
                        (2*pi)^(this.hfd.dim/2)*fTrans.C,reshape(this.hfd.C,[ones(1,this.hfd.dim),size(this.hfd.C)]),'full');
                else
                    % See commments above. Due to the additional
                    % truncation, normalization is not ensured even
                    % though we scale it by sqrt((2*pi)^dim)
                    CJointSqrt=convnc(...
                        (2*pi)^(this.hfd.dim/2)*fTrans.C,reshape(this.hfd.C,[ones(1,this.hfd.dim),size(this.hfd.C)]),'same');
                end
                % Do NOT simply convolve with 'same' to prevent truncation
                % along the preserved dimension. We pad along these
                % dimensions to preserve the relevant entries. Any truncation along
                % these dimensions in the identity representation could
                % otherwise lead to negative function values. The entries
                % along the other dimensions are to be discarded as the
                % dimensions are marginalized out. An additional factor
				% (2*pi)^dim is required for preserving the normalization in the
                % the marginalization.
                % It is faster to use this padding and then call convn with
                % 'valid' than to calculate 'full' and then truncate.
                CJointSize=size(CJointSqrt);
                CJointPadded=zeros([3*CJointSize(1:this.hfd.dim)-2,CJointSize(this.hfd.dim+1:2*this.hfd.dim)]);
                indices=[arrayfun(@(currSize){currSize:2*currSize-1},CJointSize(1:this.hfd.dim)),repmat({':'},1,this.hfd.dim)];
                CJointPadded(indices{:})=CJointSqrt;
                CPredictedId=(2*pi)^(this.hfd.dim)*convnc(CJointPadded,CJointSqrt,'valid');
            end
            if strcmp(this.hfd.transformation,'identity')||~truncateJointSqrt
                % If identity transformation is used or no truncation has
                % been performed, no normalization issues are expected
                % (except numerical issues)
                warning(warnStruct); % Reenable warnings first as should be no problem with identity
                this.hfd=HypertoroidalFourierDistribution(CPredictedId,'identity');
            else % Normalization is not ensured. Disable warning afterward
                this.hfd=HypertoroidalFourierDistribution(CPredictedId,'identity');
                warning(warnStruct);
            end
            
            if strcmp(fTrans.transformation,'sqrt')
                this.hfd=this.hfd.transformViaFFT('sqrt',dimC(dimC>1)); % DimC(dimC>1) as workaround for 1-D
            end
        end

        function updateNonlinear(this, likelihood, z) %measurement z, likelihood(z,x)=P(Z|X)
            % Performs and update for an arbitrary likelihood function and
            % a measurement. If the measurement z is not given, assume that
            % likelihood (for varying x) is given as a hfd. Otherwise,
            % transform it.
            %
            % Parameters:
            %   likelihood f(z|x)
            %       Either given as HypertoroidalFourierDistribution or as
            %       function. If given as a function, we assume it takes 
            %       d x #points matrices (same convention as .pdf) as input 
            %       for both measurement and state.
            %   measurement z
            %       Used as input for likelihood. Is repmatted if
            %       likelihood is to be evaluated at multiple points.
            
            if nargin==2 % Can use likelihood without z, but then we assume it is given as a FourierDistribution
                assert(isa(likelihood,'HypertoroidalFourierDistribution'));
            else % If z is given, assume likelihood is a function
                likelihood=HypertoroidalFourierDistribution.fromFunction(...
                    ...% Could assume that likelihood can use implicit expansion
                    @(varargin)reshape(likelihood(repmat(z,[1,numel(varargin{1})]),cell2mat(cellfun(@(c){c(:)'},varargin)')),size(varargin{1})),...
                    size(this.hfd.C),this.hfd.transformation);
            end
            this.hfd=this.hfd.multiply(likelihood,size(this.hfd.C));
        end
        
    end
    
end

