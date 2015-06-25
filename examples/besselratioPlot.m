% This example illustrates different approximations for the ratio of Bessel
% functions

%% Calculate values
xBessel = 0:0.01:10;
yBesselMatlab = arrayfun(@(a) besseli(1,a,1)/besseli(0,a,1), xBessel);
yBesselAmos = arrayfun(@(a) besselratio(0,a), xBessel);
yBesselApproxStienne12 = arrayfun(@(a) besselratioApprox(0,a,'stienne12'), xBessel);
yBesselApproxStienne13 = arrayfun(@(a) besselratioApprox(0,a,'stienne13'), xBessel);

%% Plot results
plot(xBessel, yBesselAmos, 'b');
hold on
plot(xBessel, yBesselApproxStienne12, 'r');
plot(xBessel, yBesselApproxStienne13, 'k');
plot(xBessel(1:20:end), yBesselMatlab(1:20:end), 'kx');
hold off
xlabel('\kappa')
ylabel('I_1(\kappa)/I_0(\kappa)');
legend('Amos','Stienne12','Stienne13','Matlab','location','southeast')

