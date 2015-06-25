classdef AxialKalmanFilter < AbstractAxialFilter
    % Kalman Filter with modifications to handle quaternions, 
    % works for antipodally symmetric complex numbers as well.
        
    properties
        gauss
    end
    
    methods
        function this = AxialKalmanFilter()
            % Constructor
            this.setState(GaussianDistribution([1;0;0;0],eye(4,4)));
        end
        
        function setState(this, g)
            % Sets the state
            assert(isa(g, 'GaussianDistribution'));
            assert(size(g.mu,2)==1, 'mu must be column vector');
            assert(size(g.mu,1)==2 || size(g.mu,1)==4, 'only 2d and 4d supported');
            assert(abs(norm(g.mu)-1)<1E-5, 'mean must be a unit vector');
            this.gauss = g;
            this.d = size(g.mu,1);
            this.setCompositionOperator();
        end
        
        function predictIdentity(this, gaussW)           
            % Perform prediction with system noise gaussW
            assert(abs(norm(gaussW.mu)-1)<1E-5, 'mean must be a unit vector');
            mu_ = this.compositionOperator(this.gauss.mu, gaussW.mu);
            C_ = this.gauss.C + gaussW.C;
            this.gauss = GaussianDistribution(mu_, C_);
        end
        
        function updateIdentity(this, gaussV, z)
            % Perform update with measurement z and measurement noise (muV, Cv)
            assert(abs(norm(gaussV.mu)-1)<1E-5, 'mean must be a unit vector');
            muVconj = [gaussV.mu(1); -gaussV.mu(2:end)];
            z = this.compositionOperator(muVconj, z);
            
            if dot(z,this.gauss.mu)<0 %mirror z if necessary
                z= -z;
            end
            
            d = this.d;
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
            g = this.gauss;
        end
            
    end
    
end
