classdef (Abstract) AbstractDistribution
    % Abstract base class for all distributions
    
    properties (SetAccess = protected) % Cannot be immutable becaues it needs to be set in a subclass constructor
        dim {mustBePositive, mustBeInteger, mustBeScalarOrEmpty}
    end
    
    methods (Abstract)
        % Evaluate pdf at positions stored in xa
        pdf(this, xa);
        mean(this);
    end
    
    methods
        function s = sample(this, n)
            % Obtain n samples from the distribution
            % use metropolics hastings by default
            % individual distributions can override this with more
            % sophisticated solutions
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       one sample per column
            arguments
                this (1,1) AbstractDistribution
                n (1,1) {mustBePositive,mustBeInteger}
            end
            s = this.sampleMetropolisHastings(n);
        end
        
        function s = sampleMetropolisHastings(this, n, proposal, startPoint, burnIn, skipping)
            % Metropolis Hastings sampling algorithm
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       one sample per column
            %
            % Hastings, W. K. 
            % Monte Carlo Sampling Methods Using Markov Chains and Their Applications 
            % Biometrika, 1970, 57, 97-109
            arguments
                this (1,1) AbstractDistribution
                n (1,1) {mustBePositive,mustBeInteger}
                proposal (1,1) function_handle % Default proposals are set in inherting classes
                startPoint (:,1) double % Default starting points are set in inherting classes
                burnIn (1,1) double = 10
                skipping (1,1) double = 5
            end
            
            totalSamples = burnIn+n*skipping;
            s = NaN(this.dim,totalSamples);
            % Start at mean
            % Start at mean or mode or some other reasonable point
            x = startPoint;
            % A better proposal distribution could be obtained by roughly estimating
            % the uncertainty of the true distribution.
            i=1;
            pdfx = this.pdf(x);
            while i<=totalSamples
                xNew = proposal(x); %generate new sample
                pdfxNew = this.pdf(xNew);
                a = pdfxNew/pdfx;
                if a>1 || a>rand(1) % Due to short circuiting of ||, rand is only evaluated if the first condition is not met
                    %keep sample
                    s(:,i)=xNew;
                    x = xNew;
                    pdfx = pdfxNew;
                    i=i+1;
                else
                    %reject sample
                end
            end
            s = s(:,burnIn+1:skipping:end);
        end
    end
    
end