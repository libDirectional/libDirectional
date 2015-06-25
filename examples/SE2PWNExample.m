%% Example for using the SE2 PWN Distributioon

%% Parameters
si1squared = 0.3;
si2squared = 1;
si3squared = 1;

rho12 = 0.9;
rho13 = 0.8;
rho23 = 0.5;
si12 = rho12*sqrt(si1squared*si2squared);
si13 = rho13*sqrt(si1squared*si3squared);
si23 = rho23*sqrt(si2squared*si3squared);
C = [si1squared si12 si13;
    si12 si2squared si23;
    si13 si23 si2squared];
mu = [1,20,10]';

%% Experiment with distribution
% Create distribution
p = SE2PWNDistribution(mu, C);

% check covariance
covarianceAnalytical = p.covariance4D
covarianceNumerical = p.covariance4DNumerical
covarianceError = covarianceAnalytical-covarianceNumerical

% stocahstic sampling
s = p.sample(50000);
p2 = SE2PWNDistribution.fromSamples(s);

%% create plot of samples
rotateFirst = true;

% apply rotation to translation
% this is a different interpretation (first rotate, then translate)
if rotateFirst
    for i=1:size(s,2)
        phi = s(1,i);
        R = [cos(phi), -sin(phi); sin(phi), cos(phi)];
        s(2:3, i) = R*s(2:3, i);
    end
end

%% plot samples
figure(1)
quiver(s(2,:),s(3,:),cos(s(1,:)),sin(s(1,:)));
hold on
if rotateFirst
    phi = mu(1);
    R = [cos(phi), -sin(phi); sin(phi), cos(phi)];
    muNew(1) = mu(1);
    muNew(2:3) = R*mu(2:3);
    quiver(muNew(2),muNew(3),cos(muNew(1)),sin(muNew(1)), 'color','red', 'linewidth', 2);
else
    quiver(mu(2),mu(3),cos(mu(1)),sin(mu(1)), 'color','red', 'linewidth', 2);
end
hold off
axis equal
xlabel('x')
ylabel('y')

%% plot 2D slice of density
figure(2)
step = 0.2;
plotSize = 3;
[x,y] = meshgrid([(mu(2)-plotSize):step:(mu(2)+plotSize)], [(mu(3)-plotSize):step:(mu(3)+plotSize)]);
z = zeros(size(x));
for i=1:size(x,1)
    for j=1:size(x,2)
        z(i,j) = p.pdf([mu(1), x(i,j), y(i,j)]');
    end
end
surf(x,y,z)
xlabel('x')
ylabel('y')
zlabel('f(\mu_1, x, y)');

