classdef ConstrainedUKF < AbstractCircularFilter
    % A modified unscented Kalman filter for circular
    % distributions, interprets circle as constraint on R^2.
    %
    % The interface does not expose the 2D view, so it is compatible with the other
    % circular filters.
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
        function this = ConstrainedUKF()
            % Constructor
            state_ = GaussianDistribution(1, 1);
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
           this.state = ConstrainedUKF.gauss1Dto2D(state_);
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
            gaussSys = ConstrainedUKF.gauss1Dto2D(gaussSys);
            
            function y=g(x)
                x=x/norm(x);
                tmp = f(atan2(x(2),x(1)));
                y = [cos(tmp);sin(tmp)];
            end
            model = SysModelWrapper(@g);
            model.setNoise(Gaussian(gaussSys.mu, gaussSys.C));
            
            this.ukf.setState(Gaussian(this.state.mu, this.state.C));
            this.ukf.predict(model);
            [this.state.mu, this.state.C] = this.ukf.getStateMeanAndCov();
            
            if norm(this.state.mu)~=0
                this.state.mu = this.state.mu/norm(this.state.mu);
            end
        end
        
        function updateNonlinear(this, f, gaussMeas, z)
           % Updates assuming nonlinear measurement model
            %todo: with or without mod 2 pi?
            % z(k) = f(x(k)) + v_k
            % where z is the measurement and v_k is the measurement noise
            % given by gaussMeas
            % 
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi) to R^d, where R^d is
            %       the measurement space containing z
            %   z (d x 1 vector)
            %       measurement            
            if ~isa(gaussMeas,'GaussianDistribution')
                gaussMeas = gaussMeas.toGaussian();
            end
            %gaussMeas = ConstrainedUKF.gauss1Dto2D(gaussMeas);
            
            %move measurement if necessary
            %if abs(this.state.mu - z) > pi
            %    z = z + 2*pi*sign(this.state.mu-z);
            %end            
            
            %UKF Update
            function y=g(x)
                x=x/norm(x);
                y = f(atan2(x(2),x(1))); %todo: is this smart? would it not make more sense to convert the measurement to 2D?
            end
            model = MeasModelWrapper(@g);
            model.setNoise(Gaussian(gaussMeas.mu, gaussMeas.C));
            
            this.ukf.setState(Gaussian(this.state.mu, this.state.C));
            this.ukf.update(model, z);
            [this.state.mu, this.state.C] = this.ukf.getStateMeanAndCov();
            
            if norm(this.state.mu)~=0 
                this.state.mu = this.state.mu/norm(this.state.mu);
            else
                error('mu was 0');
            end
        end
        
        function state = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   state (GaussianDistribution)
            %       current estimate
            state = ConstrainedUKF.gauss2Dto1D(this.state);
        end
    end
    
    methods (Static)
       
        function g2d = gauss1Dto2D(g1d)
            % Convert 1D Gaussian to 2D Gaussian.
            % The 1D Gaussian is interpreted as a Gaussian on the unit
            % circle, the 2D gaussian is a constrained Gaussien with
            % ||mu||=1 in the plane.
            %
            % Parameters:
            %   g1d (GaussianDistribution)
            %       one-dimensional Gaussian to approximate
            %   g2d (GaussianDistribution)
            %       two-dimensional approximation
            assert(isa(g1d, 'GaussianDistribution'));
            assert(all(size(g1d.mu) == [1 1]));
            mu = [cos(g1d.mu);sin(g1d.mu)];
            
            fx  = @(phi) g1d.pdf(phi).*(cos(phi)-mu(1)).^2;
            fxy = @(phi) g1d.pdf(phi).*(cos(phi)-mu(1)).*(sin(phi)-mu(2));
            fy  = @(phi) g1d.pdf(phi).*(sin(phi)-mu(2)).^2;
            C(1,1) = integral(fx, 0, 2*pi);
            C(1,2) = integral(fxy, 0, 2*pi);
            C(2,1) = C(1,2);
            C(2,2) = integral(fy, 0, 2*pi);
            
            g2d = GaussianDistribution(mu,C);
        end
        
        function g1d = gauss2Dto1D(g2d)
            % Convert 2D Gaussian to 1D Gaussian
            % The 1D Gaussian is interpreted as a Gaussian on the unit
            % circle, the 2D gaussian is a constrained Gaussien with
            % ||mu||=1 in the plane.
            %
            % Parameters:
            %   g2d (GaussianDistribution)
            %       two-dimensional Gaussian to approximate
            %   g1d (GaussianDistribution)
            %       one-dimensional approximation
            assert(isa(g2d, 'GaussianDistribution'));
            assert(all(size(g2d.mu) == [2 1]));
            mu = mod(atan2(g2d.mu(2),g2d.mu(1)),2*pi);
                                   
            f = @(x,y) reshape(g2d.pdf([x(:)';y(:)']).* angularError(mu,atan2(y(:)',x(:)')).^2, size(x,1),size(x,2));
            C = integral2(f, -3, 3, -3, 3);
            
            g1d = GaussianDistribution(mu,C);
        end
    end
end

