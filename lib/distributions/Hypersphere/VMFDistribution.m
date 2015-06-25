classdef VMFDistribution < AbstractHypersphericalDistribution
    % Represents a von Mises-Fisher distribution.
    %
    % Notation:
    % In this class, d represents the dimension of the distribution.
    %
    % see R. Fisher, "Dispersion on a sphere," Proceedings of the Royal Society of
    % London. Series A, Mathematical and Physical Sciences, vol. 217, no. 1130,
    % pp. 295-305, 1953.
    
    properties
        kappa   % Concentration (scalar)
        mu      % mean as vector
        C       % normalization constant
    end
    
    methods
        function VFM = VMFDistribution(mu_, kappa_)
            % Constructor
            %
            % Parameters:
            %   mu_ (d x 1)
            %       location parameter (unit vector)
            %   kappa_ (scalar)
            %       concentration parameter (>=0)
            epsilon = 1E-6;
            assert(size(mu_,2) == 1, 'mu must be a row vector');
            assert(abs(norm(mu_) - 1)<epsilon, 'mu must be a normalized');
            assert(kappa_>=0, 'kappa must be postive');
            
            VFM.mu = mu_;
            VFM.kappa = kappa_;
            
            VFM.d = size(mu_,1);
            VFM.C = kappa_^(VFM.d/2-1)/ ( (2*pi)^(VFM.d/2) * besseli(VFM.d/2-1, kappa_));
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
            assert (size(xa,1) == size(this.mu,1));
            
            p = zeros(1, size(xa,2));
            for i=1:size(xa,2)
                x = xa(:,i);
                p(i) = this.C * exp( this.kappa * this.mu' * x);
            end
        end
        
        function vmf = multiply(this, other)
            assert(isa(other, 'VMFDistribution'));
            assert(length(this.mu) == length(other.mu));
            
            mu_ = this.kappa*this.mu + other.kappa*other.mu;
            kappa_ = norm(mu_);
            mu_ = mu_ /kappa_;
            vmf = VMFDistribution(mu_,kappa_);
        end

    end
    
    methods (Static)
        function V = fit(samples, weights)
            % Fits VMF parameters to a set of samples
            % Parameters:
            %   samples (d x n matrix)
            %       matrix that contains one sample per column
            % Returns:
            %   V (VMFDistribution)
            %       the MLE estimate for a VMF distribution given the
            %       samples
            %
            % Sra, S. 
            % A Short Note on Parameter Approximation for von Mises--Fisher
            % Distributions: And a Fast Implementation of Is (x)
            % Computational Statistics, Springer-Verlag, 2012, 27, 177-190

            d = size(samples,1);    
            n = size(samples,2);
            if nargin<2
                weights = ones(1,n)/n;
            end
            assert(size(weights,1)==1, 'weights needs to be a row vector');
            assert(size(weights,2)==n, 'number of weights and samples needs to match');
            assert(abs(sum(weights)-1) < 0.001, 'weights must sum to 1'); %check normalization
            
            vecSum = sum(samples*diag(weights),2);
            mu_ = vecSum/norm(vecSum);
            Rbar = norm(vecSum);
            kappa_= VMFDistribution.AdInverse(d, Rbar);
            
            V = VMFDistribution(mu_, kappa_);
        end
        
        function result = Ad(d, kappa)
            % attention Ad is defined differently then in
            % besselratio/besselratioinverse
            result = besseli(d/2, kappa)/besseli(d/2-1, kappa);
        end
        
        function result = AdInverse(d, x)
            % approximation by Banaerjee et al.
            kappa_ = x*(d-x^2)/(1-x^2);

            % newton approximation by Sra 2012
            newtonStep = @(kappa) kappa - (VMFDistribution.Ad(d, kappa)-x)/(1- VMFDistribution.Ad(d, kappa)^2 - (d-1)/kappa*VMFDistribution.Ad(d, kappa));
            
            maxSteps = 20;
            epsilon = 1E-7; % stopping condition
            
            for i=1:maxSteps
                kappaOld = kappa_;
                kappa_ = newtonStep(kappa_);
                if abs(kappa_ - kappaOld) < epsilon
                    break;
                end
            end
            
            result = kappa_;
        end
    end
end