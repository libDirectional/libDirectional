classdef WatsonDistribution < AbstractHypersphericalDistribution
    % Represents a Watson distribution.
    %
    % Watson, G. S. 
    % Equatorial Distributions on a Sphere
    % Biometrika, 1965, 52, 193-201
    
    properties
        kappa (1,1) double  % concentration (scalar)
        mu (:,1) double      % mean as a unit vector
        C {mustBeNonzero}    % normalization constant
    end
    
    methods
        function W = WatsonDistribution(mu_, kappa_)
            % Constructor
            %
            % Parameters:
            %   mu_ (d x 1)
            %       location parameter (unit vector)
            %   kappa_ (scalar)
            %       concentration parameter (>=0)
            arguments
                mu_ (:,1) double
                kappa_ (1,1) {mustBeReal} % Can be negative
            end
            epsilon = 1E-6;
            assert(size(mu_,2) == 1, 'mu must be a row vector');
            assert(abs(norm(mu_) - 1)<epsilon, 'mu must be normalized');
            
            W.mu = mu_;
            W.kappa = kappa_;
            
            W.dim = size(mu_,1);
            W.C = gamma(W.dim/2)/2/pi^(W.dim/2)/hypergeom(1/2, W.dim/2, W.kappa);
        end
        
        function p = pdf(this, xa)
            % Evaluates pdf at each column of xa
            % Parameters:
            %   xa (d x n matrix)
            %       each column represents one of the n points in R^d that the
            %       pdf is evaluated at; can be just a (d x 1) vector as well
            % Returns:
            %   p (1 x n row vector)
            %       values of the pdf at each column of xa
            assert(size(xa,1) == this.dim);
            p = this.C * exp( this.kappa * (this.mu' * xa).^2);
        end      
        
        function B = toBingham(this)
            % Converts this distribution to a Bingham distribution. The
            % conversion is exact because a Watson distribution is an
            % isotropic Bingham distribution.
            %
            % Returns:
            %   B (BinghamDstribution)
            %       the Bingham distribution with identical pdf
            if this.kappa<0
                error ('conversion to Bingham is not implemented yet for kappa<0');
            end
            
            M = repmat(this.mu, 1, this.dim);
            % make vectors in M linear independent
            E = eye(this.dim, this.dim); 
            E(1,1)=0;
            M = M + E;
            % make vectors in M orthogonal
            [Q,~] = qr(M);
            M = [Q(:,2:end)  Q(:,1)]; % put mu at in the last column
            Z = [repmat(-this.kappa,this.dim-1,1); 0];
            B = BinghamDistribution(Z,M);
        end
        
        function s = sample(this,n)
            % Generate samples from Watson distribution.
            %
            % Parameters:
            %   n (scalar)
            %       number of samples to generate
            % Returns:
            %   s (d x n matrix)
            %       generated samples (one sample per column)
            if this.dim~=3
                s = this.toBingham.sample(n); % use Bingham sampling
            else
                % Algorithm LW from the Paper "Random Sampling From the
                % Watson Distribution" by Kim-Hung Li and Carl Ka-Fai Wong
                s=NaN(this.dim,n);
                theta=NaN(1,n);
                phi=NaN(1,n);
                k=this.kappa;
                rho=4*k/(2*k+3+sqrt((2*k+3)^2-16*k));
                r=(3*rho/(2*k))^3*exp(-3+2*k/rho);
                i=0;
                while i<n
                    U=rand(n-i,3);
                    S=U(:,1).^2./(1-rho*(1-U(:,1).^2));
                    V=r*U(:,2).^2./(1-rho*S).^3;
                    valid=V<=exp(2*k*S);
                    if ~any(valid)
                        continue
                    end
                    U=U(valid,:);
                    S=S(valid);
                    
                    thetasNew=acos(sqrt(S));
                    U3ltHalf=U(:,3)<0.5;
                    thetasNew=pi*U3ltHalf+(-1).^(U3ltHalf).*thetasNew;
                    % The 2*pi*(2*U(:,3)-1) in the otherwise case in the
                    % paper is equivalent as this yields 4*pi*U(:,3)-2*pi,
                    % which can be neglected due to the periodicity.
                    phisNew=4*pi*U(:,3); 
                    theta(i+1:i+size(U,1))=thetasNew;
                    phi(i+1:i+size(U,1))=phisNew;
                    i=i+size(U,1);
                end
                Rz=@(alpha)[cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1];
                Ry=@(beta)[cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)];
                [muphi,mutheta]=cart2sph(this.mu(1),this.mu(2),this.mu(3));
                [s(1,:),s(2,:),s(3,:)]=sph2cart(phi,-theta+pi/2,1);
                s=Rz(muphi-pi)*Ry(mutheta+pi/2)*s;
            end
        end
        
        function m = mode(this)
            % Calculate the mode of a Watson distribution
            % Returns:
            %   m (d x 1 column vector)
            %       mode of the distribution (note that -m is the mode as well)
            if this.kappa>=0
                m = this.mu; %todo: this is only correct for kappa>=0
            else
                m = this.modeNumerical(); % fallback to numerical
            end
        end

        function dist = setMode(this, newMode)
            arguments
                this (1,1) WatsonDistribution
                newMode (:,1) double
            end
            assert(isequal(size(newMode),size(this.mu)));
            dist = this;
            dist.mu = newMode;
        end

        function distShifted = shift(this, offsets)
            % There is no true shifting for the hypersphere. This is a function for compatibility and only works when mu is [0,0,...,1].
            arguments
                this (1,1) WatsonDistribution
                offsets (:,1) double {mustBeNonempty}
            end
            assert(isequal(this.mu,[zeros(this.dim-1,1);1]), 'There is no true shifting for the hypersphere. This is a function for compatibility and only works when mu is [0,0,...,1].');
            distShifted = this.setMode(offsets);
        end
    end
end