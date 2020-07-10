classdef CustomDistribution < AbstractDistribution
    % Distribution with custom pdf.
    
    properties
        f function_handle
    end
    
    methods
        function this = CustomDistribution(f_,dim_)
            arguments
                f_ function_handle
                dim_
            end
            % Constructor, it is the user's responsibility to ensure that f is a valid
            % density and takes arguments of the same form as
            % .pdf, i.e., it needs to be vectorized. Further, the user has
            % to ensure the density is normalized.
            % 
            % Parameters:
            %   f_ (function handle)
            %       pdf of the distribution
            %   dim_ (scalar)
            %       dimension of the real space in which the hypersphere is
            %       embedded
            this.dim = dim_;
            this.f = f_;
        end
        
        function p = pdf(this, xa)
            % Evaluate pdf at each column of xa.
            %
            % Parameters:
            %   xa (dim x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            assert(size(xa,1)==this.dim);
            p = this.f(xa);
            assert(isequal(size(p),[1,size(xa,2)])); % Validate output format of pdf is as expected
        end
        
    end

end
