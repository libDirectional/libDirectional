% This example analyzes the results obtained from the  Metropolis-Hastings 
% algorithm. If the burn-in and/or skipping paramters are chosen too small, 
% a significant correlation can be introduced in the samples.

wn = WNDistribution(2,1.2);

%% one sample per call
n = 1000;
s = zeros(1,n);
for i=1:n
    s(i) = wn.sampleMetropolisHastings(1);
end

% plot
figure(1)
scatter(s(1:end-1), s(2:end));
setupAxisCircular('x','y');
xlabel('sample n')
ylabel('sample n+1')

% calculate circular correlation
twd = ToroidalWDDistribution([s(1:end-1); s(2:end)]);
twd.circularCorrelationJammalamadaka

%% all samples at once
s2 = wn.sampleMetropolisHastings(n);

% plot
figure(2)
scatter(s2(1:end-1), s2(2:end));
setupAxisCircular('x','y');
xlabel('sample n')
ylabel('sample n+1')

% calculate circular correlation
twd2 = ToroidalWDDistribution([s2(1:end-1); s2(2:end)]);
twd2.circularCorrelationJammalamadaka
