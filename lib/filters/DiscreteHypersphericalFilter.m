classdef DiscreteHypersphericalFilter < AbstractHypersphericalFilter
    % A discrete filter on the hypersphere 
    
    properties
        d
        w
        DT
    end
    
    methods
        function this = DiscreteHypersphericalFilter(nParticles, dim, method)           
            % Constructor
            %
            % Parameters:
            %   nParticles (integer > 0)
            %       number of particles to use
            if nargin < 3
                method = 'equalarea';
            end            
            this.d = DiscreteHypersphericalFilter.computeDiscretization(nParticles, dim, method);
            nParticles = size(this.d,2); %the actual number of particles may be different from the requested number of particles
            this.w = ones(1,nParticles)/nParticles;
            %this.DT = delaunayTriangulation(this.d') ;
            this.DT = convhulln(this.d');
        end
        
        function setState(this, distribution)
            % Sets the current system state
            %
            % Parameters:
            %   distribution (AbstractHypersphericalDistribution)
            %       new state
            assert(isa(distribution, 'AbstractHypersphericalDistribution'));
            this.w = distribution.pdf(this.d);
            this.w = this.w/sum(this.w);
        end
        
        function d = dim(this)
            d = size(this.d,1);
        end
        
        function predictNonlinear(this, f, noiseDistribution)
            % where w(k) is noise given by noiseDistribution.
            % Predicts assuming a nonlinear system model, i.e.,
            % f(x(k+1)|x_k) = VMF(x(k+1) ; mu = f(x_k), kappa_k^w)
            % where w is noise given by sysNoise.
            %
            % Parameters:
            %   f (function handle)
            %       function from S^(dim-1) to S^(dim-1)
            %   noiseDistribution (VMFDistribution)
            %       distribution of noise
            assert(isa (noiseDistribution, 'VMFDistribution')); %todo generalize to Watson etc.
            assert(isa(f,'function_handle'));
                        
            %apply f
            nParticles = length(this.d);
            w_ = zeros(size(this.w)); % set all weights to zero
            for i=1:nParticles
                newLocation = f(this.d(:,i)); 
                noiseDistribution.mu = newLocation;
                w_ = w_ + this.w(i)*noiseDistribution.pdf(this.d);
            end
            
            this.w = w_;
        end
        
        function predictNonlinearNonAdditive(this, f, noiseSamples, noiseWeights, mode)       
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k), w(k))
            % where w(k) is non-additive noise given by samples and weights.
            %
            % Parameters:
            %   f (function handle)
            %       function from S^(dim-1) x W to S^(dim-1) (W is the space
            %       containing the noise samples)
            %   noiseSamples (d2 x n matrix)
            %       n samples of the noise as d2-dimensional vectors
            %   noiseWeights (1 x n vector)
            %       weight of each sample
            
            % (samples, weights) are discrete approximation of noise
            assert(size(noiseWeights,1) == 1, 'weights most be row vector')
            assert(size(noiseSamples,2) == size(noiseWeights,2), 'samples and weights must match in size');
            assert(isa(f,'function_handle'));
            
            noiseWeights = noiseWeights/sum(noiseWeights); % ensure normalization
            
            if nargin <= 4
                mode = 'nn';
            end

            % apply f
            nParticles = length(this.d);
            oldWeights = this.w;
            this.w = zeros(size(this.w)); % set all weights to zero
            for i=1:nParticles
                for j=1:length(noiseSamples)
                    newLocation = f(this.d(:,i),noiseSamples(:,j));
                    newWeight = oldWeights(i)*noiseWeights(j);
                    if strcmp(mode, 'nn')
                        %simple nearest neighbor assignment for now
                        idx = knnsearch(this.d', newLocation');
                        this.w(idx) = this.w(idx) + newWeight;
                    elseif strcmp(mode, 'delaunay')
                        %todo: only works in 3D so far?
                        if this.dim == 3
                            %find triangle
                            [t,u,v] = this.findTriangle(newLocation);
                            %interpolate in triangle
                            this.w(t(1)) = this.w(t(1)) + (1-u-v)*newWeight;
                            this.w(t(2)) = this.w(t(2)) + u*newWeight;
                            this.w(t(3)) = this.w(t(3)) + v*newWeight;
                        else
                            error('not supported')
                        end
                    elseif strcmp(mode, 'knn')
                        k=3;
                        [idx,D] = knnsearch(this.d', newLocation', 'K', k);
                        %todo test, verify that weights sum to 1
                        %todo this can lead to negative weights!
                        for l=1:k 
                            this.w(idx(l)) = this.w(idx(l)) + (sum(D)-(k-1)*D(l))/sum(D)*newWeight;
                        end
                    else
                        error('unsupported mode');
                    end
                end
            end            
        end
        
        function [triangle, u, v] = findTriangle(this, x)
            [int, ~, u, v] = TriangleRayIntersection([0,0,0], x, this.d(:,this.DT(:,1)), this.d(:,this.DT(:,2)), this.d(:,this.DT(:,3)));
            index = find(int,1);
            assert(~isempty(index)); %if this fails, no triangle was found due to numerical inaccuracy - how to handle that?
            triangle = this.DT(index,:);
            u = u(index);
            v = v(index);
        end
        
        function updateNonlinear(this, likelihood, z) 
            % Updates assuming nonlinear measurement model given by a
            % likelihood function likelihood(z,x) = f(z|x), where z is the
            % measurement. The function can be created using the
            % LikelihoodFactory.
            % 
            % Parameters:
            %   likelihood (function handle)
            %       function from Z x [0,2pi) to [0, infinity), where Z is
            %       the measurement space containing z
            %   z (arbitrary)
            %       measurement
            
            assert(isa(likelihood,'function_handle'));
            
            % You can either use a likelihood depending on z and x
            % and specify the measurement as z or use a likelihood that
            % depends only on x and omit z.
            if nargin==2
                this.w = this.w .* likelihood(this.d);
            else
                this.w = this.w .* likelihood(z, this.d);
            end
        end
        
        function [d,w] = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            d = this.d;
            w = this.w;
        end
        
        function plotEstimate(this)
            switch this.dim
                case 2
                    scatter3(this.d(1,:), this.d(2,:), this.w)
                case 3
                    scatter3(this.d(1,:), this.d(2,:), this.d(3,:), length(this.w)*10*this.w)
                otherwise
                    error('not implemented')
            end
        end
        
        function m = getEstimateMean(this)
            m = sum(this.d.*repmat(this.w, this.dim, 1),2);
            m = m/norm(m);
        end
    end
    
    methods (Static)
        function d = computeDiscretization(nParticles, dim, method)
            function r = goalFun (d)
                d = d./repmat(sqrt(sum(d.^2)),size(d,1),1);
                [~, dist] = knnsearch(d',d', 'K', 2);
                r = -sum((dist(:,2)));
            end            
            
            if strcmp(method, 'random')
                d = rand(dim, nParticles)-0.5;
                d = d./repmat(sqrt(sum(d.^2)),dim,1);
            elseif strcmp(method, 'stratified')
                assert(dim == 3);
                d = RandSampleSphere(nParticles, 'stratified')';
            elseif strcmp(method, 'optimize')
                if dim == 2
                    % use simple closed form solution in 2D
                    phi = linspace(0,2*pi,nParticles + 1);
                    phi = phi(1:end-1);
                    d = [cos(phi);sin(phi)];
                else
                    [pathstr,~,~] = fileparts(mfilename('fullpath'));
                    filename = sprintf('%s/../util/autogenerated/optimize-%i-%i.mat', pathstr, nParticles, dim);
                    if exist(filename, 'file')
                        load (filename);
                    else
                        %initialize with random values
                        %d = rand(dim, nParticles)-0.5;
                        %d = d./repmat(sqrt(sum(d.^2)),dim,1);
                        %initialize with equal area solution
                        d = eq_point_set(dim-1, nParticles);
                        %todo improve solution quality and performance
                        d = fminunc(@goalFun, d, optimoptions('fminunc', 'display', 'iter', 'MaxFunEvals', 1E6));
                        d = d./repmat(sqrt(sum(d.^2)),dim,1);
                        save(filename, 'd');
                    end
                end
            elseif strcmp(method, 'platonic')
                assert(dim == 3);
                %todo generate next larger polyhedron instead of insisting
                %on number of particles
                switch nParticles
                    case 4
                        d = DiscreteHypersphericalFilter.tetrahedron();
                    case 6
                        d = DiscreteHypersphericalFilter.octahedron();
                    case 8
                        d = DiscreteHypersphericalFilter.hexahedron();
                    case 12
                        d = DiscreteHypersphericalFilter.icosahedron();
                    case 20
                        d = DiscreteHypersphericalFilter.dodecahedron();
                    otherwise
                        error('unsupported number of particles');
                end
            elseif strcmp(method, 'spiral')
                assert(dim == 3);
                d = SpiralSampleSphere(nParticles)';
            elseif strcmp(method, 'reisz')
                assert(dim == 3);
                if nParticles<14 
                    nParticles = 14;
                end
                d = ParticleSampleSphere('N', nParticles)';
            elseif strcmp(method, 'cubesubdivision')
                %number of points: 6*4.^n+2
                assert(dim == 3);
                c = QuadCubeMesh();
                while size(c.vertices, 1) < nParticles
                    % subdivide untile we have at least nParticles
                    c = SubdivideSphericalMesh(c,1);
                end
                d = c.vertices';
            elseif strcmp(method, 'pentakissubdivision')
                %number of points 30*4.^n+2
                assert(dim == 3);
                c = DodecahedronMesh();
                while size(c.Points, 1) < nParticles
                    % subdivide untile we have at least nParticles
                    c = SubdivideSphericalMesh(c,1);
                end
                d = c.Points';  
            elseif strcmp(method, 'rhombicsubdivision')
                %number of points 12*4.^n+2
                assert(dim == 3);
                c = QuadRhombDodecMesh();
                while size(c.vertices, 1) < nParticles
                    % subdivide untile we have at least nParticles
                    c = SubdivideSphericalMesh(c,1);
                end
                d = c.vertices';  
            elseif strcmp(method, 'icosahedronsubdivision')
                %number of points: 10*4.^n+2
                %https://oeis.org/A122973
                assert(dim == 3);
                c = IcosahedronMesh();
                while size(c.Points, 1) < nParticles
                    % subdivide untile we have at least nParticles
                    c = SubdivideSphericalMesh(c,1);
                end
                d = c.Points';  
            elseif strcmp(method, 'equalarea')
                d = eq_point_set(dim-1, nParticles);
            elseif strcmp(method, 'gridsphere')
                assert(dim == 3);
                [lat,lon] = GridSphere(nParticles);
                lat = lat/180*pi;
                lon = lon/180*pi;

                x = cos(lat) .* cos(lon);
                y = cos(lat) .* sin(lon);
                z = sin(lat);
                d = [x'; y'; z'];
            elseif strcmp(method, 'schaefer')
                assert(dim == 4);
                %will return 2*8^n points for n=1,2,3,...
                d = tessellate_S3(nParticles);
            elseif strcmp(method, 'regular4polytope')
                assert(dim == 4);
                %todo implement remaining 4 polytopes?
                if nParticles <=24
                    d = DiscreteHypersphericalFilter.generate24cell();
                elseif nParticles <=120
                    d = DiscreteHypersphericalFilter.generate600cell();
                elseif nParticles<=600
                    d = DiscreteHypersphericalFilter.generate120cell();
                else
                    d = DiscreteHypersphericalFilter.generate120cell();
                    %todo warning/error?
                end
            elseif strcmp(method, 'tesseractsubdivision')
                assert(dim == 4);
                d = tesseractsubdivision(nParticles)';
            else
                error('unsupported method')
            end
        end
        
        function d = tetrahedron()
            % see https://en.wikipedia.org/wiki/Tetrahedron#Formulas_for_a_regular_tetrahedron
            d = [1 0 -1/sqrt(2);
                -1 0 -1/sqrt(2);
                0  1  1/sqrt(2);
                0 -1  1/sqrt(2)]';
            d = d./repmat(sqrt(sum(d.^2)),size(d,1),1);          
        end
        
        function d = hexahedron()
            % see https://en.wikipedia.org/wiki/Cube#Cartesian_coordinates
            d = [1 1 1;
                 1 1 -1;
                 1 -1 1;
                 1 -1 -1;
                 -1 1 1;
                 -1 1 -1;
                 -1 -1 1;
                 -1 -1 -1]'/sqrt(3);            
        end
        
        function d = octahedron()
            % see https://en.wikipedia.org/wiki/Octahedron#Cartesian_coordinates
            d = [1 0 0 ;
                 0 1 0;
                 0 0 1]';
            d = [d -d];
        end
        
        function d = dodecahedron()
            % see https://en.wikipedia.org/wiki/Dodecahedron#Cartesian_coordinates
            h = (sqrt(5)-1)/2; % inverse golden ratio
            a = 1+h;
            b = 1-h^2;
            d = [1 1 1;
                 1 1 -1;
                 1 -1 1;
                 1 -1 -1;
                 -1 1 1;
                 -1 1 -1;
                 -1 -1 1;
                 -1 -1 -1;
                 0 a b;
                 0 a -b;
                 0 -a b;
                 0 -a -b;
                 a, b, 0;
                 a, -b, 0;
                 -a, b, 0;
                 -a, -b, 0;
                 b, 0, a;
                 b, 0, -a;
                 -b, 0, a;
                 -b, 0, -a;
                 ]';     
            d = d./repmat(sqrt(sum(d.^2)),size(d,1),1);
        end
   
        function d = icosahedron()
            % see https://en.wikipedia.org/wiki/Regular_icosahedron#Cartesian_coordinates
            phi = (1+sqrt(5))/2; % golden ratio
            d = [0 1 phi;
                1 phi 0;
                phi 0 1;
                0 -1 phi;
                -1 phi 0;
                phi 0 -1]';
            d = [d -d];
            d = d./repmat(sqrt(sum(d.^2)),size(d,1),1);
        end
        
        function d = generate600cell()
            % Generates 600 cell in 4D
            % see https://en.wikipedia.org/wiki/600-cell
            % The 600 cell has 120 vertices!
            phi = (1+sqrt(5))/2; % golden ratio
            signs3 = 2*(dec2bin(0:7)-'0')-1;
            p = evenperms(4);
            
            %16 vertices on [+-0.5, +-0.5, +-0.5, +-0.5]
            d1 = dec2bin(0:2^4-1)-'0'-0.5;
            %8 vertices on permutations of [0 0 0 +-1]
            d2 = dec2bin(2.^(0:3))-'0';
            d3 = -d2;
            %96 vertices on even permutations of 0.5 [+-phi, +-1, +-1/phi,0]
            d4 = zeros(96,4);
            index = 1;
            for j=1:size(p,1) %iterate over permutations
                for i=1:size(signs3,1) %iterate over 8 different combinations of signs
                    currentPoint = 0.5* [ [phi , 1 , 1/phi].*signs3(i,:),  0];
                    d4(index,:) = currentPoint(p(j,:));
                    index = index + 1;
                end
            end
            d = [d1; d2; d3; d4]';
            d = d./repmat(sqrt(sum(d.^2)),size(d,1),1);
        end
        
        function d = generate120cell()
            % Generates 120 cell in 4D
            % see https://en.wikipedia.org/wiki/120-cell
            % and http://mathworld.wolfram.com/120-Cell.html
            % The 120 cell has 600 vertices!
            phi = (1+sqrt(5))/2; % golden ratio
            signs3 = 2*(dec2bin(0:7)-'0')-1;
            signs4 = 2*(dec2bin(0:15)-'0')-1;
            p = evenperms(4);
            shiftRight = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0];
            
            % permutations of [0, 0, +-2, +-2]
            %d1 = [zeros(4,2) (dec2bin(0:3)-'0')*4-2];
            d1 = [perms([0 0 -2 -2]); perms([0 0 -2 2]); perms([0 0 2 2])];
            d1 = unique(d1,'rows');
            
            % permutations of  [+-1, +-1, +-1,, +-sqrt(5)]
            d2 = signs4.*repmat([1 1 1 sqrt(5)], 16, 1);
            d2 = [d2; d2*shiftRight; d2*shiftRight^2; d2*shiftRight^3];
            
            % permutations of [+- phi^(-2), +-phi, +-phi, +-phi]
            d3 = signs4.*repmat([phi^(-2), phi, phi, phi], 16, 1);
            d3 = [d3; d3*shiftRight; d3*shiftRight^2; d3*shiftRight^3];
            
            % permutations of [+- phi^(-1), +- phi^(-1), +-phi^(-1), +-phi^2]
            d4 = signs4.*repmat([phi^(-1), phi^(-1), phi^(-1), phi^2], 16, 1);
            d4 = [d4; d4*shiftRight; d4*shiftRight^2; d4*shiftRight^3];
            
            %even permutations of [0, +-phi^(-2), +-1, +-phi^2]
            d5 = zeros(96,4);
            index = 1;
            for j=1:size(p,1) %iterate over permutations
                for i=1:size(signs3,1) %iterate over 8 different combinations of signs
                    currentPoint = [0,  [phi^(-2) , 1 , phi^2].*signs3(i,:)];
                    d5(index,:) = currentPoint(p(j,:));
                    index = index + 1;
                end
            end            
            
            %even permutations of [0, +-phi^(-1), +-phi, +-sqrt(5)]
            d6 = zeros(96,4);
            index = 1;
            for j=1:size(p,1) %iterate over permutations
                for i=1:size(signs3,1) %iterate over 8 different combinations of signs
                    currentPoint = [0,  [phi^(-1) , phi , sqrt(5)].*signs3(i,:)];
                    d6(index,:) = currentPoint(p(j,:));
                    index = index + 1;
                end
            end           
            
            %even permutations of [+-phi^(-1), +-1, +-phi, +-2]
            d7 = zeros(192,4);
            index = 1;
            for j=1:size(p,1) %iterate over permutations
                for i=1:size(signs4,1) %iterate over 16 different combinations of signs
                    currentPoint = [phi^(-1), 1, phi, 2].*signs4(i,:);
                    d7(index,:) = currentPoint(p(j,:));
                    index = index + 1;
                end
            end             
            
            d = [d1; d2; d3; d4; d5; d6; d7]';
            d = d./repmat(sqrt(sum(d.^2)),size(d,1),1);
        end
        
        function d = generate24cell()
            % Generates 24 cell in 4D
            % see https://en.wikipedia.org/wiki/24-cell
            % The 24 cell has 24 vertices!
            signs4 = 2*(dec2bin(0:15)-'0')-1;
            shiftRight = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0];
            
            %permutations of [+-1, 0, 0, 0]
            d1 = [1 0 0 0; -1 0 0 0];
            d1 = [d1; d1*shiftRight; d1*shiftRight^2; d1*shiftRight^3];
            
            %[+-0.5,+-0.5,+-0.5,+-0.5]
            d2 = signs4 * 0.5;
            
            d = [d1; d2;]';
            d = d./repmat(sqrt(sum(d.^2)),size(d,1),1);
        end
    end
    
end

function p = evenperms(n)
    %returns even permutations of 1,..,n
    p = perms(1:n);
    M = eye(n,n);
    e = false(1,size(p,1));
    for j=1:size(p,1) %iterate over permutations
        currentPerm = p(j,:);
        if det(M(currentPerm,:))==1 %check if permutation is even
            e(j)=1;
        end
    end
    p = p(e,:);
end
