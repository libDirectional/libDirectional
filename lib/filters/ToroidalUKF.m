classdef ToroidalUKF < AbstractToroidalFilter
    % A modified unscented Kalman filter for torus distributions, 
    % interprets the torus as Cartesian product [0, 2pi)^2 of the interval
    % [0,2pi).
    %
    % Gerhard Kurz, Igor Gilitschenski, Maxim Dolgov, Uwe D. Hanebeck,
    % Bivariate Angular Estimation Under Consideration of Dependencies Using Directional Statistics
    % Proceedings of the 53rd IEEE Conference on Decision and Control (CDC 2014),
    % Los Angeles, California, USA, December 2014.
    
    properties
        state
        ukf
    end
    
    methods
        function this = ToroidalUKF()
            % Constructor
            state_ = GaussianDistribution([0;0], eye(2,2));
            this.setState(state_)
            this.ukf = UKF();
        end   
        
        function setState(this, state_)
            % Sets the current system state
            %
            % Parameters:
            %   state_ (GaussianDistribution)
            %       new state (2D Gaussian)
            if ~isa(state_,'GaussianDistribution')
                state_ = state_.toGaussian();
            end
            assert(size(state_.mu,1)==2);
            assert(size(state_.mu,2)==1);
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
            %       function from [0,2pi)^2 to [0,2pi)^2
            %   gaussSys (GaussianDistribution)
            %       distribution of additive noise
            assert(isa(f,'function_handle'));
            if ~isa(gaussSys,'GaussianDistribution')
                gaussSys = gaussSys.toGaussian();
            end
                        
         	model = SysModelWrapper(@(x) f(x));
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
            %   z (2 x 1 vector)
            %       measurement in [0, 2pi)^2
            assert(size(z,1) == 2);
            assert(size(z,2) == 1);
            if ~isa(gaussMeas,'GaussianDistribution')
                gaussMeas = gaussMeas.toGaussian();
            end
            z = mod(z - gaussMeas.mu, 2*pi);
            
            %move measurement if necessary
            for dim=1:2
                if abs(this.state.mu(dim) - z(dim)) > pi
                    z(dim) = z(dim) + 2*pi*sign(this.state.mu(dim)-z(dim));
                end
            end
            
            %Kalman update
            K = this.state.C/(this.state.C + gaussMeas.C);
            mu = this.state.mu + K*(z-this.state.mu);
            C = (eye(2,2)-K)*this.state.C;
            
            mu = mod(mu,2*pi);
            this.state = GaussianDistribution(mu,C);
        end
        
        function updateNonlinear(this, f, gaussMeas, z, measurementPeriodic)
            % Updates assuming nonlinear measurement model
            % z(k) = f(x(k)) + v_k
            % or 
            % z(k) = f(x(k)) + v_k mod 2 pi (in case measurementPeriodic is
            % true)
            % where z is the measurement and v_k is the measurement noise
            % given by gaussMeas
            % 
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi)^2 to R^d, where R^d is
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

