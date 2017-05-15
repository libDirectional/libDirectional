classdef GvMDistribution < AbstractCircularDistribution
    % Generalized von Mises distribution of arbitrary order
    %
    % see Riccardo Gatto, Sreenivasa Rao Jammalamadaka ,
    % "The Generalized von Mises Distribution",
    % Statistical Methodology, 2007
    
    properties
        mu
        kappa
    end
    
    methods
        function this = GvMDistribution(mu_,kappa_)
            % Constructor
            %
            % Parameters
            %   mu_ (k x 1 vector)
            %       location parameters
            %   kappa_ (k x 1) vector
            %       concentration parameters
            assert(size(mu_,2)==1);
            assert(isequal(size(mu_),size(kappa_)));
            assert(all(kappa_>0));
            this.mu=mu_;
            this.kappa=kappa_;
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
            
            % Because this is a value class, we can only store values as attributes in the constructor
            normConst=integral(@(x)this.pdfUnnormalized(x),0,2*pi); 
            p=1/normConst*this.pdfUnnormalized(xa);
        end
        
        function p = pdfUnnormalized(this, xa)
            % Evaluate pdf without respective the normalization constant at each column of xa
            %
            % Parameters:
            %   xa (1 x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            assert(size(xa,1)==1);
            p=exp(sum(bsxfun(@times,cos(bsxfun(@times,bsxfun(@minus,xa,this.mu),(1:numel(this.mu))')),this.kappa),1));
        end
    end
end

