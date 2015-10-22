classdef GeneralCircularMixture < AbstractCircularDistribution
    % Mixture of multiple circular distributions. The distributions may belong
    % to different classes.
    
    properties
        cds     % Cell array of circular distributions
        w       % Weights
    end
    
    methods
        function this = GeneralCircularMixture(cds, w)
            % Constructor
            assert(isa(cds,'cell') && all(cellfun(@(cd)isa(cd,'AbstractCircularDistribution'),cds)),...
                'cds must be a cell array of circular distributions');
            assert(all(size(cds) == size(w)),'size of cds and w must be equal');
            
            if all(cellfun(@(cd)isa(cd,'FourierDistribution'),cds))
                warning('Mixtures of Fourier distributions can be built by combining the Fourier coefficients so using a mixture may not be necessary');
            end
            if all(cellfun(@(cd)isa(cd,'WDDistribution'),cds))
                warning('Mixtures of WDDistributions can usually be combined into one WDDistribution.');
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
                p = p + this.w(i)*this.cds{i}.pdf(xa); % Calculate pdf using individual pdfs
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
                m = m + this.w(i)*this.cds{i}.trigonometricMoment(n); % Calculate moments using moments of each component
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
            occurrences=histc(d,1:length(this.cds));
            count=1;
            for i=1:length(this.cds)
                s(count:count+occurrences(i)-1) = this.cds{d(i)}.sample(occurrences(i));
                count=count+occurrences(i);
            end
        end
    end
end

