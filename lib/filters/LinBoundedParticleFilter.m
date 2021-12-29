classdef LinBoundedParticleFilter < AbstractLinBoundedFilter & AbstractParticleFilter
    % SIR particle filter for the hypercylinder
    methods
        function this = LinBoundedParticleFilter(nParticles, boundD, linD)
            % Constructor
            %
            % Parameters:
            %   nParticles (integer > 0)
            %       number of particles to use 
            %   dim (integer >0)
            %       dimension
            arguments
                nParticles (1,1) double {mustBeInteger, mustBePositive}
                boundD (1,1) double {mustBeInteger, mustBeNonnegative}
                linD (1,1) double {mustBeInteger, mustBeNonnegative}
            end
            dim = boundD+linD;
            this.dist = LinBoundedDiracDistribution(boundD, zeros(dim, nParticles), 1/nParticles*ones(1,nParticles)); 
        end
        
        function setState(this, dist_)
            % Sets the current system state
            %
            % Parameters:
            %   distribution (AbstractHypertoroidalDistribution)
            %       new state
            arguments
                this (1,1) LinBoundedParticleFilter
                dist_ (1,1) AbstractLinBoundedDiracDistribution
            end
            assert(isequal(dist_size,this.dist.size));
            this.dist = dist_;
        end
    end
end

