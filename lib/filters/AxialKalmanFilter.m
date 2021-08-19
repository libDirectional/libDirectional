classdef AxialKalmanFilter < AbstractAxialFilter
    % Kalman Filter with modifications to handle quaternions, 
    % works for antipodally symmetric complex numbers as well.
    %
    % Gerhard Kurz, Igor Gilitschenski, Simon Julier, Uwe D. Hanebeck,
    % Recursive Bingham Filter for Directional Estimation Involving 180 Degree Symmetry
    % Journal of Advances in Information Fusion, 9(2):90 - 105, December 2014.    
        
    properties
        gauss
    end
    
    methods
        function this = AxialKalmanFilter()
            % Constructor
            this.setState(GaussianDistribution([1;0;0;0],eye(4,4)));
        end
        
        function setState(this, g)
            % Sets the current system state
            %
            % Parameters:
            %   g (GaussianDistribution)
            %       new state
            assert(isa(g, 'GaussianDistribution'));
            assert(size(g.mu,2)==1, 'mu must be column vector');
            assert(size(g.mu,1)==2 || size(g.mu,1)==4, 'only 2d and 4d supported');
            assert(abs(norm(g.mu)-1)<1E-5, 'mean must be a unit vector');
            this.gauss = g;
            this.setCompositionOperator();
        end
        
        function predictIdentity(this, gaussW)           
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) (+) w(k)    
            % where w(k) is noise given by gaussW.
            % The composition operator (+) refers to a complex or quaternion
            % multiplication.
            %
            % Parameters:
            %   gaussW (GaussianDistribution)
            %       distribution of noise
            assert(isa(gaussW, 'GaussianDistribution'));            
            assert(abs(norm(gaussW.mu)-1)<1E-5, 'mean must be a unit vector');
            mu_ = this.compositionOperator(this.gauss.mu, gaussW.mu);
            C_ = this.gauss.C + gaussW.C;
            this.gauss = GaussianDistribution(mu_, C_);
        end
        
        function updateIdentity(this, gaussV, z)           
            % Updates assuming identity measurement model, i.e.,
            % z(k) = x(k) (+) v(k)    mod 2pi,
            % where v(k) is additive noise given by vmMeas.
            % The composition operator (+) refers to a complex or quaternion
            % multiplication.
            %
            % Parameters:
            %   gaussV (GaussianDistribution)
            %       distribution of additive noise
            %   z (dim x 1 vector)
            %       measurement on the unit hypersphere            
            assert(isa(gaussV, 'GaussianDistribution'));            
            assert(abs(norm(gaussV.mu)-1)<1E-5, 'mean must be a unit vector');
            assert(size(gaussV.mu,1) == size(this.gauss.mu,1));
            assert(all(size(z) == size(this.gauss.mu)));
            assert(abs(norm(z) - 1) < 1E-5, 'measurement must be a unit vector');
            
            muVconj = [gaussV.mu(1); -gaussV.mu(2:end)];
            z = this.compositionOperator(muVconj, z);
            
            if dot(z,this.gauss.mu)<0 % mirror z if necessary
                z= -z;
            end
            
            d = this.dim;
            H = eye(d,d); % measurement matrix
            IS = H*this.gauss.C*H' + gaussV.C; %innovation covariance
            K = this.gauss.C*H'/IS; % Kalman gain
            IM = z - H*this.gauss.mu; % measurement residual
            mu_ = this.gauss.mu + K*IM; %updated mean
            C_ = (eye(d,d) - K*H)*this.gauss.C; %updated covariance
            
            mu_ =mu_/norm(mu_); % enforce unit vector
            
            this.gauss = GaussianDistribution(mu_, C_);
        end
        
        function g = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   g(GaussianDistribution)
            %       current estimate
            g = this.gauss;
        end
            
    end
    
end
