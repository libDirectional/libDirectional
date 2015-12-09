classdef HypersphericalUniformDistribution < AbstractHypersphericalDistribution
    methods
        function this=HypersphericalUniformDistribution(dim_)
            this.dim=dim_;
        end
        
        function vals=pdf(this,xa)
            % pdf is always 1 divided by surface area (integrates to 1).
            vals = (1/AbstractHypersphericalDistribution.computeUnitSphereSurface(this.dim))*ones(1,size(xa,2));
        end
        
        function X=sample(this,n)
            % General algorithm based on "A note on a method for generating 
            % points uniformly on n-dimensional spheres" by Mervin E. Muller April 1959.
            % Algorithm for 2-sphere based on "Spherical sampling by archimedes' theorem"
            % by Min-Zhi Shao and Norman Badler, 1996
            if this.d==3
                X=NaN(this.dim,n);
                phi=2*pi*rand(1,n);
                X(3,:)=rand(n,1)*2-1;
                r=sqrt(1-X(3,:).^2);
                X(1,:)=r.*cos(phi);
                X(2,:)=r.*sin(phi);
            else
                X=cell2mat(cellfun(@(x){x/norm(x)},num2cell(randn(this.dim,n),1)));
            end
        end
    end 
end