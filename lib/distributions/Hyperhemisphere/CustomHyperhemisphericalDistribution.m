classdef CustomHyperhemisphericalDistribution < CustomDistribution & AbstractHyperhemisphericalDistribution
    % Hyperhemispherical distribution with custom pdf.
    
    methods (Static)
        function chhd = fromDistribution(dist)
            % Creates a CustomHypertoroidalDistribution from some other distribution
            %
            % Parameters:
            %   dist (AbstractHypertoroidalDistribution)
            %       distribution to convert
            % Returns:
            %   chd (CustomHypersphericalDistribution)
            %       CustomHypersphericalDistribution with identical pdf
            if isa(dist,'AbstractHyperhemisphericalDistribution')
                chhd = CustomHyperhemisphericalDistribution(@(xa)dist.pdf(xa),dist.dim);
            elseif isa(dist,'BinghamDistribution')
                % For antipodally symmetric densities on the hypersphere, can just cut in half
                % and double to preserve normalization
                chhd = CustomHyperhemisphericalDistribution(@(xa)2*dist.pdf(xa),dist.dim);
            elseif isa(dist,'AbstractHypersphericalDistribution')
                warning('FromDistribution:UsePdfHypersphere',...
                    ['You are creating a CustomHyperhemispherical distribution based on a distribution on the full hypersphere. ',...
                    'Using numerical integration to calculate the normalization constant.']);
                chhdUnnorm = CustomHyperhemisphericalDistribution(@(xa)dist.pdf(xa),dist.dim);
                normConstInv = chhdUnnorm.integral;
                chhd = CustomHyperhemisphericalDistribution(@(xa)dist.pdf(xa)./normConstInv,dist.dim);
            else
                error('Input variable dist is of wrong class');
            end 
        end
    end
end
