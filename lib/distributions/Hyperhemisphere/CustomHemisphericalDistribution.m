classdef CustomHemisphericalDistribution < CustomHyperhemisphericalDistribution
    % Hyperhemispherical distribution with custom pdf.
    methods
        function this = CustomHemisphericalDistribution(f_)
            arguments
                f_ (1,1) function_handle
            end
            this@CustomHyperhemisphericalDistribution(f_, 3);
        end
    end
    
    methods (Static)
        function chsd = fromDistribution(dist)
            % Creates a CustomHypertoroidalDistribution from some other distribution
            %
            % Parameters:
            %   dist (AbstractHypertoroidalDistribution)
            %       distribution to convert
            % Returns:
            %   chd (CustomHypersphericalDistribution)
            %       CustomHypersphericalDistribution with identical pdf
            assert(dist.dim==3);
            if isa(dist,'AbstractHyperhemisphericalDistribution')
                chsd = CustomHemisphericalDistribution(@(xa)dist.pdf(xa));
            elseif isa(dist,'BinghamDistribution')
                % For antipodally symmetric densities on the hypersphere, can just cut in half
                % and double to preserve normalization
                chsd = CustomHemisphericalDistribution(@(xa)dist.pdf(xa));
                chsd.scaleBy = 2;
            elseif isa(dist,'AbstractHypersphericalDistribution')
                warning('FromDistribution:UsePdfHypersphere',...
                    ['You are creating a CustomHyperhemispherical distribution based on a distribution on the full hypersphere. ',...
                    'Using numerical integration to calculate the normalization constant.']);
                chhdUnnorm = CustomHyperhemisphericalDistribution(@(xa)dist.pdf(xa),dist.dim);
                normConstInv = chhdUnnorm.integral();
                chsd = CustomHemisphericalDistribution(@(xa)dist.pdf(xa));
                chsd.scaleBy=1/normConstInv;
            else
                error('Input variable dist is of wrong class');
            end 
        end
    end
end
