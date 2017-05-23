classdef FourierFilter < AbstractCircularFilter
    % Fourier filters.  Depending on the parameter in the initialization, 
    % either Fourier identity or Fourier square root filter is used.
    %
    % see Florian Pfaff, Gerhard Kurz, and Uwe D. Hanebeck,
    % Multimodal Circular Filtering Using Fourier Series
    % Proceedings of the 18th International Conference on Information Fusion (Fusion 2015), 
    % Washington, D.C., USA, July 2015.
    %
    % Florian Pfaff, Gerhard Kurz, Uwe D. Hanebeck,
    % Nonlinear Prediction for Circular Filtering Using Fourier Series
    % Proceedings of the 19th International Conference on Information Fusion (Fusion 2016), 
    % Heidelberg, Germany, July 2016.
    
    properties
        fd
    end
    
    methods
        function this = FourierFilter(noOfCoefficients,transformation)
            % Constructor
            %
            % Parameters:
            %   noOfCoefficients (integer > 0)
            %       number of Fourier coefficients to use
            if nargin==1,transformation='sqrt';end
            assert(strcmp(transformation,'sqrt')||strcmp(transformation,'identity'));
            this.fd = FourierDistribution.fromDistribution(CircularUniformDistribution(),noOfCoefficients,transformation);
        end
        
        function setState(this, fd_)
            % Sets the current system state
            %
            % Parameters:
            %   fd_ (FourierDistribution)
            %       new state
            assert(isa(fd_,'AbstractCircularDistribution'));
            if ~(isa(fd_,'FourierDistribution'))
                warning('setState:nonFourier','fd_ is not a FourierDistribution. Transforming with a number of coefficients that is equal to that of the filter.');
                fd_=FourierDistribution.fromDistribution(fd_,2*length(this.fd.a)-1,this.fd.transformation);
            else
                if ~strcmp(this.fd.transformation,fd_.transformation)
                    warning('setState:transDiffer','New density is transformed differently.');
                end
                if ~(length(fd_.a)==length(this.fd.a))
                    warning('setState:noOfCoeffsDiffer','New density has different number of coefficients.')
                end
            end
            this.fd = fd_;
        end
        
        function predictIdentity(this,dSys)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by dSys.
            %
            % Parameters:
            %   dSys (AbstractCircularDistribution)
            %       distribution of additive noise
            assert(isa(dSys,'AbstractCircularDistribution'));
            if ~isa(dSys,'FourierDistribution')
                warning('Predict:automaticConversion',...
                    'dSys is not a FourierDistribution. Transforming with a number of coefficients that is equal to that of the filter. For non-varying noises, transforming once is much more efficient and should be perferred.');
                dSys=FourierDistribution.fromDistribution(dSys,2*length(this.fd.a)-1,this.fd.transformation);
            end
            this.fd=this.fd.convolve(dSys);  
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
            assert(isa(dMeas,'AbstractCircularDistribution'));
            if ~isa(dMeas,'FourierDistribution')
                warning('Update:automaticConversion',...
                    'dMeas is not a FourierDistribution. Transforming with a number of coefficients that is equal to that of the filter. For non-varying noises, transforming once is much more efficient and should be perferred.');
                dMeas=FourierDistribution.fromDistribution(dMeas,2*length(this.fd.a)-1,this.fd.transformation);
            end
            dMeasNew=dMeas.shift(z);
            fdtmp=this.fd.multiply(dMeasNew);
            this.fd=fdtmp.truncate(2*length(this.fd.a)-1);
        end
        
        function fd = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   fd (FourierDistribution)
            %       current estimate
            fd = this.fd;
        end
        
        function predictNonlinear(this, f, noiseDistribution)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k)    mod 2pi,
            % where w(k) is additive noise given by noiseDistribution.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi) to [0,2pi)
            %   noiseDistribution (AbstractCircularDistribution)
            %       distribution of additive noise
            
            assert(isa (noiseDistribution, 'AbstractCircularDistribution'));
            assert(isa(f,'function_handle'));
            
            fTrans=@(xk,xkk)reshape(noiseDistribution.pdf(xkk(:)'-f(xk(:)')),size(xk));
            this.predictNonlinearViaTransitionDensity(fTrans);

        end
        
        function predictNonlinearViaTransitionDensity(this, fTrans,truncateJointSqrt)
            % Predicts assuming a nonlinear system model using a
            % probabilistic model.
            %
            % Parameters:
            %   fTrans
            %       transition density f(x(k+1)|x(k))
            if nargin==2,truncateJointSqrt=true;end
            
            warnStruct=warning('off','Normalization:notNormalized'); % fTrans might be unnormalized and the multiplication result below will not be normalized either
            if ~isa(fTrans,'ToroidalFourierDistribution')
                assert(nargin(fTrans)==2);
                fTrans=ToroidalFourierDistribution.fromFunction(fTrans,(2*length(this.fd.a)-1)*[1,1],this.fd.transformation);
            else
                assert(strcmp(this.fd.transformation,fTrans.transformation));
            end
            
            if strcmp(this.fd.transformation,'identity')
                cPredictedId=conv2(fTrans.C,this.fd.c.','valid');
            elseif strcmp(this.fd.transformation,'sqrt')
                % There are two ways to ensure nonnegativity of the result.
                % 1) Calculating the full coefficient vectors/tensors for 
                % the identity transformation for both ('full' is needed as
                % coefficint vectors must not be truncated) and then by
                % convolving them with 'valid' is sufficient to obtain the
                % result of the marginalization. This can be implemented as
                % conv2(conv2(fTrans.C,fTrans.C,'full'),conv2(this.fd.c.',this.fd.c.','full'),'valid')
                % 2) Calculating the (possibly truncated) result for the
                % coefficient tensor for the square root of the joint and
                % then convolve the result with itself. We use this
                % approach as it is more flexible and featured a better run
                % time performance (even if we chose to not truncate the
                % joint density).
                if truncateJointSqrt
                    cJointSqrt=conv2(fTrans.C,this.fd.c.','same');
                else
                    cJointSqrt=conv2(fTrans.C,this.fd.c.','full'); 
                end
                additionalColumns=2*length(this.fd.b); % If we set this value to length(this.fd.b) instead of 2*length(this.fd.b), we would only get the desired number of coefficients but we could no longer guarantee that the function is always nonnegative.
                cPredictedId=conv2(... % It is faster to use this padding and then call conv2 with 'valid' than to calculate 'full' and then truncate
                        [zeros(size(cJointSqrt,1),additionalColumns),cJointSqrt,zeros(size(cJointSqrt,1),additionalColumns)],...
                        cJointSqrt,'valid');
            end
            fdPredictedId=FourierDistribution.fromComplex(cPredictedId,'identity');
            warning(warnStruct); % Reenable warnings
            if strcmp(this.fd.transformation,'sqrt')
                this.fd=fdPredictedId.transformViaFFT('sqrt',2*length(this.fd.a)-1);
            else
                this.fd=fdPredictedId;
            end
        end
        
        function mean=getEstimateMean(this)
            mean=this.fd.circularMean;
        end
        
    end
    
end

