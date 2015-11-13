classdef CustomHypertoroidalDistribution < AbstractHypertoroidalDistribution
    % Hypertoroidal distribution with custom pdf.
    
    properties
        f
    end
    
    methods
        function this = CustomHypertoroidalDistribution(f_,dim)
            % It is the user's responsibility to ensure that f is a valid
            % hypertoroidal density and takes arguments of the same form as
            % .pdf
            assert(isa(f_, 'function_handle'));
            this.dim=dim;
            this.f = f_;
        end
        
        function p = pdf(this, xa)
            p=this.f(xa);
        end
        function hd=shift(this,shiftAngles)
            hd=this;
            hd.f=@(xa)this.f(xa-repmat(shiftAngles,[1,size(xa,2)]));
        end
    end
    methods (Static)
        function chd=fromDistribution(dist)
            chd=CustomHypertoroidalDistribution(@(xa)dist.pdf(xa),dist.dim);
        end
    end
    
end
