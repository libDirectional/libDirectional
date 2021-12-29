classdef (Abstract) AbstractUniformDistribution < AbstractDistribution
    % Uniform distribution on the hypertorus.
    
    methods
        function p = pdf(this,xa)
            % Evaluate pdf at each column of xa.
            %
            % Parameters:
            %   xa (dim x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location            
            p = (1/this.getManifoldSize)*ones(1,size(xa,2));
        end 
    end
    methods (Sealed)
        function m = mode(~) %#ok<STOUT>
            error('Mode not available for uniform distribution');
        end
    end
end