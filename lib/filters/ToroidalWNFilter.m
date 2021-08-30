classdef ToroidalWNFilter < AbstractToroidalFilter
    % A filter based on the ToroidalWNDistribution.
    %  
    % Gerhard Kurz, Igor Gilitschenski, Maxim Dolgov, Uwe D. Hanebeck,
    % Bivariate Angular Estimation Under Consideration of Dependencies Using Directional Statistics
    % Proceedings of the 53rd IEEE Conference on Decision and Control (CDC 2014),
    % Los Angeles, California, USA, December 2014.
    % 
    % Gerhard Kurz, Florian Pfaff, Uwe D. Hanebeck,
    % Nonlinear Toroidal Filtering Based on Bivariate Wrapped Normal Distributions (to appear)
    % Proceedings of the 20th International Conference on Information Fusion (Fusion 2017), 
    % Xi'an, China, July 2017.
    
    properties
        twn ToroidalWNDistribution
    end
    
    methods
        function this = ToroidalWNFilter()
            % Constructor
            twn_ = ToroidalWNDistribution([0;0],eye(2,2));
            this.setState(twn_);
        end
        
        function setState(this, twn_)
            % Sets the current system state
            %
            % Parameters:
            %   twn_ (ToroidalWNDistribution)
            %       new state
            assert(isa(twn_, 'ToroidalWNDistribution'));
            this.twn = twn_;
        end
        
        function predictIdentity(this, twnSys)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by twnSys.
            % The modulo operation is carried out componentwise.
            %
            % Parameters:
            %   twnSys (ToroidalWNDistribution)
            %       distribution of additive noise
            assert(isa(twnSys, 'ToroidalWNDistribution'));
            this.twn = this.twn.convolve(twnSys);
        end
            
        function predictNonlinear(this, f, twnSys) %system function f
            twd = this.twn.toToroidalWD5();
            twdF = twd.applyFunction(f);
            try
                twnF = twdF.toToroidalWN();
            catch
                twnF = twdF.toToroidalWNmixedMLE();
            end
            this.twn = twnF.convolve(twnSys);
        end
        
        function updateIdentity(this, twnMeas, z)
            % Updates assuming identity measurement model, i.e.,
            % z(k) = x(k) + v(k)    mod 2pi,
            % where v(k) is additive noise given by twnMeas.
            % The modulo operation is carried out componentwise.
            %
            % Parameters:
            %   twnMeas (ToroidalWNDistribution)
            %       distribution of additive noise
            %   z (2 x 1 vector)
            %       measurement in [0, 2pi)^2
            assert(isa(twnMeas, 'ToroidalWNDistribution'));
            assert(size(z,1) == 2);
            assert(size(z,2) == 1);
            
            muWnew = mod(z - twnMeas.mu,2*pi);
            twnMeasShifted = ToroidalWNDistribution(muWnew, twnMeas.C);
            try
                this.twn = this.twn.multiplyMomentBased(twnMeasShifted);
            catch 
                warning('multiplication failed, using previous estimate');
            end
        end
        
        function updateNonlinear(this, likelihood, z) %measurement z, likelihood(z,x)=P(Z|X)
            assert(size(z,2) == 1);
            twd = this.twn.toToroidalWD5();
            twdNew = twd.reweigh(@(x) likelihood(z,x));
            this.twn = twdNew.toToroidalWN();
        end

        function updateNonlinearParticle(this, likelihood, z)
            % Updates assuming nonlinear measurement model given by a
            % likelihood function likelihood(z,x) = f(z|x), where z is the
            % measurement. The function can be created using the
            % LikelihoodFactory.
            %
            % This method uses a nondeterministic Monte Carlo algorithm.
            % 
            % Parameters:
            %   likelihood (function handle)
            %       function from Z x [0,2pi)^2 to [0, infinity), where Z is
            %       the measurement space containing z
            %   z (arbitrary)
            %       measurement
            assert(isa(likelihood,'function_handle'));
            
            n = 200;
            samples = this.twn.sample(n);
            twd = ToroidalWDDistribution(samples);
            twdNew = twd.reweigh(@(x) likelihood(z,x));
            this.twn = twdNew.toToroidalWN();
        end
        
        function updateNonlinearProgressive(this, likelihood, z, tau) 
            % Updates assuming nonlinear measurement model given by a
            % likelihood function likelihood(z,x) = f(z|x), where z is the
            % measurement. The function can be created using the
            % LikelihoodFactory.
            %
            % This method uses a progressive algorithm.
            % 
            % Parameters:
            %   likelihood (function handle)
            %       function from Z x [0,2pi)^2 to [0, infinity), where Z is
            %       the measurement space containing z
            %   z (arbitrary)
            %       measurement
            %   tau (scalar)
            %       threshold for progression
            assert(isa(likelihood,'function_handle'));
            
            if nargin<4
                tau = 0.02; 
            else
                assert(isscalar(tau));
            end           
            lambda = 1;
            steps = 0;
            while lambda> 0
                twd = this.twn.toToroidalWD5();
                l = zeros(1, size(twd.d,2));
                for i=1:size(twd.d,2)
                    l(i) = likelihood(z,twd.d(:,i));
                end
                % avoid numerical issues by replacing zero with the smallest
                % possible double value
                l(l==0) = realmin; 
                lmin = min(l);
                lmax = max(l);
                if lmax <= realmin
                    warning('progressive update failed because likelihood is 0 everwhere')
                    return
                end                   
                currentLambda = min(log(tau)/log(lmin/lmax), lambda);
                if currentLambda <= 0
                    warning('progressive update with given threshold impossible')
                    currentLambda = 0.001;
                end                
                twdNew = twd.reweigh(@(x) likelihood(z,x).^currentLambda);
                try
                    this.twn = twdNew.toToroidalWN();
                catch
                    this.twn = twdNew.toToroidalWNmixedMLE();
                end                
                lambda = lambda - currentLambda;
                steps = steps + 1;
            end
        end

        function twn = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   twn (ToroidalWNDistribution)
            %       current estimate
            twn = this.twn;
        end
        
    end
    
end

