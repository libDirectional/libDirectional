% Example for using the discrete filter

% Initialize with 100 sampples
df = DiscreteFilter(100);

% Perform update
df.updateNonlinear(@(x) exp(-(x-3).^2));
clf
df.wd.plot2d('b');
hold on

% Perform prediction
df.predictNonlinear(@(x) x+0.1, WNDistribution(0,0.2))
df.wd.plot2d('gx');
hold off

% Some plot settigns
setupAxisCircular('x');
legend('updated density', 'predicted density');
xlabel('x'); ylabel('f(x)');
