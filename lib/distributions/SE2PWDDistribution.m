classdef SE2PWDDistribution
    % Partially wrapped dirac distribution with Dirac positions d and
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
    
    properties
        d
        w
    end
    
    methods
        function this = SE2PWDDistribution(d_, w_)
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (3 x L)
            %       Dirac locations (in [0,2pi) x R^2)
            %   w_ (1 x L)
            %       weights for each Dirac
            assert(size(d_,1) == 3 );
            d_(1,:) = mod(d_(1,:),2*pi);
            this.d = d_;
            if (nargin<2)
                %all Diracs have equal weights by default
                this.w = ones(1,size(this.d,2))/size(this.d,2);
            else
                assert(size(w_,1) == 1 );
                assert(size(d_,2) == size(w_,2));
                this.w = w_/sum(w_);
            end
        end
        
        function p = pdf(~, ~)
            % Placeholder, pdf does not exist for wrapped Dirac distributions
            p = 0; %PWD does not have a proper pdf
            warning('PDF:UNDEFINED', 'pdf is not defined')
        end

        function pwd = applyFunction(this,f)
            % Apply a function f(x) to each Dirac component and obtain its new position
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2*pi) x R^2 to [0,2*pi) x R^2
            % Returns:
            %   pwd (WDDistribution)
            %       distribution with new Dirac locations (and same
            %       weights as before)
            d_ = zeros(size(this.d));
            for i=1:size(this.d,2)
                d_(:,i) = f(this.d(:,i));
            end
            pwd = SE2PWDDistribution(d_, this.w);
        end
        
        function twd = reweigh(this, f)
            % Uses a function f(x) to calculate the weight of each Dirac
            % component. The new weight is given by the product of the old 
            % weight and the weight obtained with f. Restores normalization
            % afterwards.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2*pi) x R^2 to [0, infinity)
            % Returns:
            %   pwd (SE2PWDDistribution)
            %       distribution with new weights and same Dirac locations
            wNew = zeros(1, length(this.w));
            for i=1:length(this.d)
                wNew(i) = f(this.d(:,i));
            end
            twd = SE2PWDDistribution(this.d, wNew .* this.w);
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
                
        function wd = marginalizeLinear(this)
            % Marginalizes the linear dimensions, returns a WDDistribution
            % of the remaining circular dimension
            %
            % Returns:
            %   wd (WDDistribution)
            %       marginal distribution depending only on x_1
            wd = WDDistribution(this.d(1,:), this.w);
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

        function s = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (3 x n)
            %       n samples on [0,2pi) x R^2
            ids = discretesample(this.w,n);
            s = this.d(:,ids);
        end

    end
    
end

