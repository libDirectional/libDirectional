classdef ToroidalWNFilter < AbstractToroidalFilter
    % A filter based on the ToroidalWNDistribution.
    %  
    % Gerhard Kurz, Igor Gilitschenski, Maxim Dolgov, Uwe D. Hanebeck,
    % Bivariate Angular Estimation Under Consideration of Dependencies Using Directional Statistics
    % Proceedings of the 53rd IEEE Conference on Decision and Control (CDC 2014),
    % Los Angeles, California, USA, December 2014.
    
    properties
        twn
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

