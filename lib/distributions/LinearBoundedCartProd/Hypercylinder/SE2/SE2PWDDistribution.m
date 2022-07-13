classdef SE2PWDDistribution < HypercylindricalDiracDistribution & AbstractSE2Distribution
    % Partially wrapped Dirac distribution with Dirac positions d and
    % weights w
    % 
    % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
    % The Partially Wrapped Normal Distribution for SE(2) Estimation
    % Proceedings of the 2014 IEEE International Conference on Multisensor
    % Fusion and Information Integration (MFI 2014), Beijing, China, September 2014.
    %
    % Gerhard Kurz
    % Directional Estimation for Robotic Beating Heart Surgery
    % Karlsruhe Institute of Technology, 
    % Intelligent Sensor-Actuator-Systems Laboratory, 2015     
    
    methods
        function this = SE2PWDDistribution(d_, w_)
            arguments
                d_ (3,:) double
                w_ (1,:) double = ones(1,size(d_,2))/size(d_,2);
            end
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (3 x L)
            %       Dirac locations (in [0,2pi) x R^2)
            %   w_ (1 x L)
            %       weights for each Dirac
           this@HypercylindricalDiracDistribution(1, d_, w_);
        end
        
        function mu = mean4D(this)
            % Calculates 4D mean of [cos(x1), sin(x1), x2, x3]
            %
            % Returns:
            %   mu (4 x 1)
            %       expectation value of [cos(x1), sin(x1), x2, x3]
            s = this.d;
            S = [cos(s(1,:)); sin(s(1,:)); s(2:3,:)];
            mu = sum(repmat(this.w,4,1).*S,2);
        end
        
        function C = covariance4D(this)
            % Calculates 4D covariance of [cos(x1), sin(x1), x2, x3]
            %
            % Returns:
            %   C (4 x 4)
            %       covariance matrix of [cos(x1), sin(x1), x2, x3]
            s = this.d;
            S = [cos(s(1,:)); sin(s(1,:)); s(2:3,:)];
            C = S*diag(this.w)*S';
        end
    end
    
    methods (Static)
        function ddist = fromDistribution(dist, nParticles)
            arguments
                dist (1,1) AbstractHypercylindricalDistribution
                nParticles (1,1) {mustBeInteger,mustBePositive}
            end
            assert(dist.boundD==1,dist.linD==2);
            ddist = SE2PWDDistribution(dist.sample(nParticles), 1/nParticles * ones(1,nParticles));
        end
    end
    
end