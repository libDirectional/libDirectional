% This function performs a benchmark of the mvnpdf and the mvnpdffast
% functions.

function mvnpdfbench
    figure(1)
    benchmark();
    figure(2)
    benchmarkSingle();

    %bench(1) ;
    %benchSingle(1);
end

function benchmark
    dims = 1:20;
    evals = 1000000;
    repeats = 10;
    t1 = zeros(1,repeats);
    t2 = zeros(1,repeats);
    times = zeros(2, length(dims));
    for i = dims
        for j=1:repeats
            [t1(j), t2(j)]= bench(i, evals);
        end
        times (1,i) = median(t1);
        times (2,i) = median(t2);
    end
    benchPlot(dims, times, sprintf('time for %i evaluations with one call', evals))
end

function benchmarkSingle
    dims = 1:20;
    evals = 1000;
    repeats = 10;
    t1 = zeros(1,repeats);
    t2 = zeros(1,repeats);    
    times = zeros(2, length(dims));
    for i = dims
        for j=1:repeats
            [t1(j), t2(j)]= benchSingle(i, evals);
        end
        times (1,i) = median(t1);
        times (2,i) = median(t2);
    end
    benchPlot(dims, times, sprintf('time for %i evaluations with individual calls', evals))
end

function benchPlot(dims, times, titletext)
    subplot(2,1,1);
    plot(dims, times(1,:));
    hold on
    plot(dims, times(2,:));
    hold off
    legend('mvnpdf', 'mvnpdffast', 'location', 'northwest')
    xlabel('dimension');
    ylabel('time (s)');
    title(titletext);
    
    subplot(2,1,2);
    plot(dims, times(1,:)./times(2,:));
    hold on
    plot(dims, times(1,:)./times(1,:), 'b--');
    hold off
    xlabel('dimension');    
    ylabel('speedup');
end

function [t1, t2]= bench(n, evals)
    % n-D
    rng default
    mu = rand(1,n);
    C = rand(n,n);
    C=C*C';
    x = rand(evals,n);

    tic
    r = mvnpdf(x, mu, C); %#ok<NASGU>
    t1 = toc;
    tic
    r = mvnpdffast(x, mu, C); %#ok<NASGU>
    t2 = toc;
    fprintf('%fs\n%fs\n', t1, t2);
end

function [t1, t2]= benchSingle(n,evals)
    % n-D
    rng default
    mu = rand(1,n);
    C = rand(n,n);
    C=C*C';
    x = rand(evals,n);

    tic
    for i=1:length(x)
        r = mvnpdf(x, mu, C); %#ok<NASGU>
    end
    t1 = toc;
    tic
    for i=1:length(x)
        r = mvnpdffast(x, mu, C); %#ok<NASGU>
    end
    t2 = toc;
    fprintf('%fs\n%fs\n', t1, t2);
end