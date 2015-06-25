classdef GeneralCircularMixture < AbstractCircularDistribution
    % Mixture of multiple circular distributions. The distributions may belong
    % to different classes.
    
    properties
        cds     % Array of circular distributions
        w       % Weights
    end
    
    methods
        function this = GeneralCircularMixture(cds, w)
            % Constructor
            assert(all(size(cds) == size(w)),'size of wns and w must be equal');
            assert(all(isa(cds,'AbstractCircularDistribution')),'cds must consist only of circular distributions');
            if all(isa(cds,'FourierDistribution'))
                warning('Mixtures of Fourier distributions can be built by combining the Fourier coefficients so using a mixture may not be necessary');
            end
            if all(isa(cds,'WDDistribution'))
                warning('WDDistributions are mixtures by themselves');
            end
            this.cds = cds;
            this.w = w/sum(w);
        end
        
        function p = pdf(this, xa)
            % Evaluate pdf at each column of xa
            %
            % Parameters:
            %   xa (1 x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            assert(size(xa,1)==1);
            p = zeros(1, size(xa,2));
            for i=1:length(this.cds);
                p = p + this.w(i)*this.cds(i).pdf(xa); % Calculate pdf using individual pdfs
            end
        end
        
        function m = trigonometricMoment(this,n)
            % Calculate n-th trigonometric moment analytically
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (scalar)
            %       n-th trigonometric moment (complex number)
            m = 0;
            for i=1:length(this.cds);
                m = m + this.w(i)*this.cds(i).trigonometricMoment(n); % Calculate moments using moments of each component
            end
        end
        
        function s = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (1 x n)
            %       n samples on the circle
            
            % Sample component first, then sample from the chosen component
            d = discretesample(this.w,n);
            s = zeros(1,n);
            for i=1:n
                s(i) = this.cds(d(i)).sample(1);
            end
        end
    end
end

