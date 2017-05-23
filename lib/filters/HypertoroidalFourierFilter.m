classdef HypertoroidalFourierFilter < AbstractHypertoroidalFilter
    % Filter based on Fourier series on the hypertorus
    %
    % Florian Pfaff, Gerhard Kurz, Uwe D. Hanebeck,
    % Multivariate Angular Filtering Using Fourier Series
    % Journal of Advances in Information Fusion, December 2016.
    properties
        hfd
    end
    
    methods
        function this = HypertoroidalFourierFilter(noOfCoefficients,transformation)
            % Constructor. Pay attention that the length of
            % noOfCoefficients determines the dimensionality of the
            % underlying hypertoroidal Fourier distribution
            if nargin==1,transformation='sqrt';end
            dim=numel(noOfCoefficients);
            if dim==1
                C_=zeros(noOfCoefficients,1);
            else 
                C_=zeros(noOfCoefficients);
            end
            indexc00=num2cell((size(C_)+1)/2);
            if strcmp(transformation,'sqrt')
                C_(indexc00{:})=1/sqrt((2*pi)^dim);
            else
                C_(indexc00{:})=1/(2*pi)^dim;
            end
            this.hfd=HypertoroidalFourierDistribution(C_,transformation);
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
        
        function predictIdentity(this, hfdSys)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by hfdSys.
            % The modulo operation is carried out componentwise.
            %
            % Parameters:
            %   hfdSys (HypertoroidalFourierDistribution)
            %       distribution of additive noise
            if ~(isa(hfdSys,'HypertoroidalFourierDistribution'))
                warning('Predict:automaticConversion',...
                    'tfdSys is not a HypertoroidalFourierDistribution. Transforming with a number of coefficients that is equal to that of the filter. For non-varying noises, transforming once is much more efficient and should be perferred.');
                sizeHfdC=size(this.hfd.C); % Needed for workaround for 1D case
                hfdSys=HypertoroidalFourierDistribution.fromDistribution(hfdSys,sizeHfdC(sizeHfdC>1),this.hfd.transformation);
            end
            this.hfd=this.hfd.convolve(hfdSys,size(this.hfd.C));
        end
        
        function predictNonlinear(this, f, noiseDistribution,truncateJointSqrt)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k)    mod 2pi,
            % where w(k) is additive noise given by noiseDistribution.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi)^d to [0,2pi)^d
            %   noiseDistribution (AbstractCircularDistribution)
            %       distribution of additive noise
            assert(isa (noiseDistribution, 'AbstractHypertoroidalDistribution'));
            assert(isa(f,'function_handle'));
            if nargin==3,truncateJointSqrt=true;end
            function p=fTrans(varargin) % Can be written as lambda function but is not faster and less legible
                % Concatenate on dim+1 to ensure that tensor can be
                % correctly flattened
                ws=cat(2*this.hfd.dim+1,varargin{this.hfd.dim+1:2*this.hfd.dim})-f(cat(2*this.hfd.dim+1,varargin{1:this.hfd.dim}));
                % Flatten tensor to matrix as input for the pdf
                wsReshaped=reshape(ws,[],noiseDistribution.dim)';
                pdfvals=noiseDistribution.pdf(wsReshaped);
                p=reshape(pdfvals,size(varargin{1})); % Restore to input size
            end
            this.predictNonlinearViaTransitionDensity(@fTrans,truncateJointSqrt);
        end
        
        function predictNonlinearViaTransitionDensity(this, fTrans,truncateJointSqrt)
            % Predicts assuming a nonlinear system model using a
            % probabilistic model.
            %
            % Parameters:
            %   fTrans as HypertoroidalFourierDistribution
            %       transition density f(x(k+1)|x(k))
            if nargin==2,truncateJointSqrt=true;end
            
            warnStruct=warning('off','Normalization:notNormalized'); % fTrans might be unnormalized and the multiplication result below will not be normalized either
            if isa(fTrans,'function_handle')
                dimC=size(this.hfd.C);    % As as intermediate result for indexing as done below (necessary for dim 1)
                fTrans=HypertoroidalFourierDistribution.fromFunction(fTrans,dimC([1:this.hfd.dim,1:this.hfd.dim]),this.hfd.transformation);
            else
                assert(isa(fTrans,'HypertoroidalFourierDistribution')&&strcmp(this.hfd.transformation,fTrans.transformation)&&fTrans.dim==(2*this.hfd.dim),'predictNonlinear:fTransInvalid',...
                    'fTrans must be given as a function or a HypertoroidalFourierDistribution of appropriate dimension and transformation.');
            end
            
            if strcmp(this.hfd.transformation,'identity')
                CPredictedId=reshape(convn(fTrans.C,this.hfd.C,'valid'),size(this.hfd.C)); % Not using squeeze to support 1D
            elseif strcmp(this.hfd.transformation,'sqrt')
                % The rational is analogous to the one in the Fourier
                % filter. We can truncate and save one convolution if we
                % first calculate CJointSqrt
                if truncateJointSqrt
                    CJointSqrt=convnc(fTrans.C,this.hfd.C,'same');
                else
                    CJointSqrt=convnc(fTrans.C,this.hfd.C,'full'); 
                end
                % Pad sufficiently to prevent any truncation along non
                % marginalized dimensions
                CJointSize=size(CJointSqrt);
                CJointPadded=zeros([CJointSize(1:this.hfd.dim),3*CJointSize(this.hfd.dim+1:2*this.hfd.dim)-2]);
                indices=[repmat({':'},1,this.hfd.dim),arrayfun(@(currSize){currSize:2*currSize-1},CJointSize(this.hfd.dim+1:end))];
                CJointPadded(indices{:})=CJointSqrt;
                CPredictedId=shiftdim(convnc(CJointPadded,CJointSqrt,'valid'),this.hfd.dim);
            end
            hfdPredictedId=HypertoroidalFourierDistribution(CPredictedId,'identity');
            warning(warnStruct); % Reenable warnings
            if strcmp(this.hfd.transformation,'sqrt')
                this.hfd=hfdPredictedId.transformViaFFT('sqrt',size(this.hfd.C));
            else
                this.hfd=hfdPredictedId;
            end
        end
        
        function updateIdentity(this, hfdMeas, z)
            % Updates assuming identity measurement model, i.e.,
            % z(k) = x(k) + v(k)    mod 2pi,
            % where v(k) is additive noise given by twnMeas.
            % The modulo operation is carried out componentwise.
            %
            % Parameters:
            %   hfdMeas (HypertoroidalFourierDistribution)
            %       distribution of additive noise
            %   z (dim x 1 vector)
            %       measurement in [0, 2pi)^dim
            assert(isa(hfdMeas,'AbstractHypertoroidalDistribution'));
            if ~(isa(hfdMeas,'HypertoroidalFourierDistribution'))
                warning('Update:automaticConversion',...
                    'hfdMeas is not a HypertoroidalFourierDistribution. Transforming with an amount of coefficients that is equal to that of the filter. For non-varying noises, transforming once is much more efficient and should be perferred.');
                sizeHfdC=size(this.hfd.C); % Needed for workaround for 1D case
                hfdMeas=HypertoroidalFourierDistribution.fromDistribution(hfdMeas,sizeHfdC(sizeHfdC>1),this.hfd.transformation);
            end
            assert(isequal(size(z),[ndims(hfdMeas.C),1]));
            
            hfdShifted=hfdMeas.shift(z);
            this.hfd=this.hfd.multiply(hfdShifted,size(this.hfd.C));
        end

        function updateNonlinear(this, likelihood, z) %measurement z, likelihood(z,x)=P(Z|X)
            hfdMeas=HypertoroidalFourierDistribution.fromFunction(...
                @(varargin)reshape(likelihood(repmat(z,[1,numel(varargin{1})]),cell2mat(cellfun(@(c){c(:)'},varargin)')),size(varargin{1})),...
                size(this.hfd.C),this.hfd.transformation);
            this.updateIdentity(hfdMeas,zeros(size(z)));
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
        
    end
    
end

