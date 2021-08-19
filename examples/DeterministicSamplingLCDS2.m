%% DeterministicSamplingLCDS2

% Neat examples for on-sphere (S2) LCD sampling vor various spherical distributions. 
% Created by Daniel Frisch at ISAS for libDirectional, Nov. 2019 

% https://isas.iar.kit.edu/pdf/IFAC20_Frisch.pdf
% Daniel Frisch, Kailai Li, and Uwe D. Hanebeck
% Optimal Reduction of Dirac Mixture Densities on the 2-sphere
% Proceedings of the 1st Virtual IFAC World Congress (IFAC-V 2020), July, 2020. 


%% Dependencies
% [flist,plist] = matlab.codetools.requiredFilesAndProducts('DeterministicSamplingLCDS2.m'); [flist'; {plist.Name}']
%  - libdirectional-dev/lib/*
%  - Optimization Toolbox
%  - Deep Learning Toolbox
%  - Symbolic Math Toolbox
%  - Statistics and Machine Learning Toolbox


%% Definition & Sample Calculation

% number of deterministic samples 
nRed = 30; 

% select distribution
type = 'VonMisesFisher'; 
% VonMisesFisher, Bingham, FisherBingham, FisherBingham5b, Uniform, ProjectedGaussian, 
% Custom, Watson, SphericalHarmonics, SphericalHarmonicsComplex, AngularCentralGaussian, 
% HypersphericalMixture, ProjectedGaussianMixture 
% Note: Some of these densities are not in the public libDirectional library 
%       and may be shared on request. 

% initialize distribution object 
switch type
    
    case {'VonMisesFisher', 'VMF', 'VMFDistribution'}
        mu = [1; 0; 0];
        kappa = 7;        
        distr = VMFDistribution(mu, kappa);
        
    case {'Bingham', 'BinghamDistribution'}
        M = eye(3);
        Z = -[20; 2; 0];
        distr = BinghamDistribution(Z, M);
        
    case {'Uniform'}
        distr = HypersphericalUniformDistribution(3);
        
    case {'Custom', 'CustomHypersphericalDistribution'}
        % Cartesian to spherical coordinates  
        theta = @(c) atan2(vecnorm(c(1:2,:),2,1),c(3,:)); 
        phi   = @(c) atan2(c(2,:),c(1,:));
        % custom function
        fun = @(c) exp( -((mod(phi(c),2*pi)-mod(4*theta(c),2*pi))/1).^2 );        
        distr = CustomHypersphericalDistribution(fun, 3);
        % normalize
        const = 1/distr.integral;
        distr = CustomHypersphericalDistribution(@(c) const*fun(c), 3);
        
    case {'Watson', 'WatsonDistribution'}
        mu = [0; 1; 1];
        mu = mu/vecnorm(mu);
        kappa = 2;        
        distr = WatsonDistribution(mu, kappa);
        
    case {'SphericalHarmonics', 'SphericalHarmonicsDistributionReal'}
        unnormalizedCoeffs = rand(3, 5);
        distr = SphericalHarmonicsDistributionReal(unnormalizedCoeffs, 'sqrt');
        
    case {'SphericalHarmonicsComplex', 'SphericalHarmonicsDistributionComplex'}
        z1=.5; sigmaX=.1;
        z2=.5; sigmaY=.2;
        distr1 = SphericalHarmonicsDistributionComplex.fromFunctionFast(@(x)normpdf(x(1,:), z1, sigmaX), 11, 'sqrt'); 
        distr2 = SphericalHarmonicsDistributionComplex.fromFunctionFast(@(x)normpdf(x(2,:), z2, sigmaY), 11, 'sqrt'); 
        distr = distr1.multiply(distr2);
 
    case {'Mixture', 'HypersphericalMixture'}
        % Distribution 1: Bingham
        M = eye(3);
        Z = -[30; 3; 0];
        distr1 = BinghamDistribution(Z, M);
        % Distribution 2: Bingham
        M = eye(3);
        R = roty(90);
        Z = -[30; 3; 0];
        distr2 = BinghamDistribution(Z, R*M);
        % Mixture 
        distr = HypersphericalMixture({distr1,distr2}, [.5,.5]);
                
    otherwise
        error('Wrong density type: %s', type)
end


% select mode out of: Deterministic, Stochastic, PDFOnly
mode = 'DeterministicLCD';
switch mode
    case 'DeterministicLCD' % slow 
        % generate deterministic samples using LCD S2 Sample Reduction
        [samplesLCD,info] = distr.sampleDeterministicLCD(nRed);
        samplesRef = info.samplesRef.Cartesian;
    case 'Stochastic' % medium 
        samplesLCD = zeros(3,0);
        samplesRef = distr.sample(nRed*100);
    case 'PDFOnly' % fast, for quick view  
        samplesLCD = zeros(3,0);
        samplesRef = zeros(3,0);
    otherwise
        error('Wrong mode: %s',mode)
end

fname = sprintf('%s_nS=%u', class(distr), nRed);
fprintf('%s, %u samples \n', class(distr), nRed)


%% Plot Results on Sphere

% create & initialize figure
fig = figure(595943254); % use same figure when evaluated repeatedly 
clf(fig); 
set(fig, 'NumberTitle','off', 'Name',fname, 'Color','white')
ax = axes(fig);
set(ax, 'NextPlot','add')
set(ax, 'XGrid','on', 'YGrid','on', 'XMinorGrid','on', 'YMinorGrid','on')
axis(ax, 'image') % keep equal proportions along all axes 

% plot sphere
spherefun = @(x,y,z) x.^2 + y.^2 + z.^2 - 1;
hsph = fimplicit3(ax, spherefun, 'DisplayName','Unit Sphere');
set(hsph, 'EdgeColor','none')
colormap summer

% plot stochatic reference samples
hRef = scatter3( samplesRef(1,:), samplesRef(2,:), samplesRef(3,:), 'Parent',ax, 'DisplayName','Reference Samples' );
set(hRef, 'Marker','.', 'MarkerEdgeColor','blue', 'SizeData',30)
% hRef.Visible = false;

% plot deterministic samples
hReduce = scatter3( samplesLCD(1,:), samplesLCD(2,:), samplesLCD(3,:), 'Parent',ax, 'DisplayName','DetSamplingLCDS2' );
set(hReduce, 'Marker','.', 'MarkerEdgeColor','red', 'SizeData',1000)

% % legend? 
% lg = legend(ax, [hsph,hRef,hReduce]);
% lg.Interpreter = 'none';
% lg.FontName = 'Monospaced';
% lg.FontSize = 12;
% lg.Location = 'NorthEastOutside';

% axis? 
xlabel(ax, 'x')
ylabel(ax, 'y')
zlabel(ax, 'z')
axis off


%% Overlay Density Mesh

% get equally distributed grid points on sphere
grid_Spher = eq_point_set_polar(2, 2000); % (2 x nGrid) [theta; phi]
grid_Cart  = polar2cart(grid_Spher);     % (3 x nGrid) [x; y; z]

% get PDF at grid points 
dens = distr.pdf(grid_Cart); % (1 x nGrid)

% scaling for visualization 
height = 1.002 + 0.7*dens/max(dens); % (1 x nGrid)

% triangular mesh
DT = delaunayTriangulation(grid_Cart');
F = DT.freeBoundary(); % (nTriangles x 3)

PT = patch('Faces',F,'Vertices',grid_Cart'.*height', 'Parent',ax);
set(PT, 'FaceColor','none', 'EdgeAlpha',.8, 'LineWidth',1) % 
set(PT, 'EdgeAlpha',.6, 'LineWidth',.6)
% PT.Visible = 'on';

% view camera toolbar (better rotation of figure)
cameratoolbar(fig, 'Show');
cameratoolbar(fig, 'SetMode','orbit');
cameratoolbar(fig, 'SetCoordSys','none')



