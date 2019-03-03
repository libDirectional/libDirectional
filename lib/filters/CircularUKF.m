classdef CircularUKF < AbstractCircularFilter
    % A modified unscented Kalman filter for circular
    % distributions, interprets circle as 1D interval [0, 2pi)
    %
    % see for example
    % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
    % Recursive Bayesian Filtering in Circular State Spaces
    % arXiv preprint: Systems and Control (cs.SY), January 2015.
    
    properties
       state 
       ukf
    end
    
    methods
        function this = CircularUKF()
            % Constructor
            state_ = GaussianDistribution(0, 1);
            this.setState(state_)
            this.ukf = UKF();
        end
        
        function setState(this, state_)
            % Sets the current system state
            %
            % Parameters:
            %   state_ (GaussianDistribution)
            %       new state (1D Gaussian)
            if ~isa(state_,'GaussianDistribution')
                state_ = state_.toGaussian();
            end
            assert(length(state_.mu)==1);
            this.state = state_;
        end
        
        function predictIdentity(this, gaussSys)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by gaussSys.
            %
            % Parameters:
            %   gaussSys (GaussianDistribution)
            %       distribution of additive noise
            if ~isa(gaussSys,'GaussianDistribution')
                gaussSys = gaussSys.toGaussian();
            end
            mu = mod(this.state.mu + gaussSys.mu,2*pi);
            C = this.state.C + gaussSys.C;
            this.state = GaussianDistribution(mu, C);
        end
        
        function predictNonlinear(this, f, gaussSys)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k)    mod 2pi,
            % where w(k) is additive noise given by gaussSys.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi) to [0,2pi)
            %   gaussSys (GaussianDistribution)
            %       distribution of additive noise
            assert(isa(f,'function_handle'));
            if ~isa(gaussSys,'GaussianDistribution')
                gaussSys = gaussSys.toGaussian();
            end

         	model = SysModelWrapper(@(x) arrayfun(f, x));
            model.setNoise(Gaussian(gaussSys.mu, gaussSys.C));
            
            this.ukf.setState(Gaussian(this.state.mu, this.state.C));
            this.ukf.predict(model);
            [mu_, C_] = this.ukf.getStateMeanAndCov();
            
            this.state = GaussianDistribution(mod(mu_, 2*pi), C_);
        end

        function updateIdentity(this, gaussMeas, z) 
            % Updates assuming identity measurement model, i.e.,
            % z(k) = x(k) + v(k)    mod 2pi,
            % where v(k) is additive noise given by gaussMeas.
            %
            % Parameters:
            %   gaussMeas (GaussianDistribution)
            %       distribution of additive noise
            %   z (scalar)
            %       measurement in [0, 2pi)
            assert(isscalar(z));
            if ~isa(gaussMeas,'GaussianDistribution')
                gaussMeas = gaussMeas.toGaussian();
            end
            z = mod(z - gaussMeas.mu, 2*pi);
            
            % move measurement if necessary
            if abs(this.state.mu - z) > pi
                z = z + 2*pi*sign(this.state.mu-z);
            end
            
            % Kalman update
            K = this.state.C/(this.state.C + gaussMeas.C);
            mu = this.state.mu + K*(z-this.state.mu);
            C = (1-K)*this.state.C;
            
            mu = mod(mu,2*pi);
            this.state = GaussianDistribution(mu,C);
        end
        
        function updateNonlinear(this, f, gaussMeas, z, measurementPeriodic)
            % Updates assuming nonlinear measurement model
            % z(k) = f(x(k)) + v_k
            % or 
            % z(k) = f(x(k)) + v_k mod 2 pi (if  measurementPeriodic is true)
            % where z is the measurement and v_k is the measurement noise
            % given by gaussMeas
            % 
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi) to R^d, where R^d is
            %       the measurement space containing z
            %   z (d x 1 vector)
            %       measurement
            %   measurementPeriodic (boolean)
            %       is the measurement a periodic quantity?
            assert(size(z,2)==1);
            assert(isa(f,'function_handle'));
            
            if ~isa(gaussMeas,'GaussianDistribution')
                gaussMeas = gaussMeas.toGaussian();
            end
            
            if nargin < 5
                measurementPeriodic = false;
            end
            
            if measurementPeriodic
                %move measurement if necessary
                for dim=1:length(z)
                    if abs(this.state.mu(dim) - z(dim)) > pi
                        z(dim) = z(dim) + 2*pi*sign(this.state.mu(dim)-z(dim));
                    end      
                end
            end   
            
            %UKF Update
            this.ukf.setState(Gaussian(this.state.mu, this.state.C));
            model = MeasModelWrapper(f);
            model.setNoise(Gaussian(gaussMeas.mu, gaussMeas.C));
            this.ukf.update(model, z);
            [this.state.mu, this.state.C] = this.ukf.getStateMeanAndCov();
            
            this.state.mu = mod(this.state.mu,2*pi);
        end
        
        function state = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   state (GaussianDistribution)
            %       current estimate
            state = this.state;
        end
        function mean=getEstimateMean(this)
            mean=this.state.mu;
        end
        
    end
    
end
