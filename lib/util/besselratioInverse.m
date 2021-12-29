function kappa = besselratioInverse(v, x, type)
    % The inverse of the Ratio I_{v+1}(x)/I_v(x) is computed, where
    % I_v(x) is the Bessel-I function of order v evaluated at x.
    %
    % Parameters:
    %   v (scalar)
    %       order
    %   x (scalar)
    %       where to evaluate
    %   type (string)
    %       can be used to choose between several algorithms
    arguments
        v (1,1) double
        x (1,1) double
        type char = 'sraamos'
    end
    
    if x == 0
        kappa = 0;
        return
    end
    
    if strcmp(type,'sraamos')
        assert (v==0);
        kappa = besselInverseSraAmos(x);
    elseif strcmp(type,'fisher')
        assert(v==0);
        kappa = besselInverseFisher(x);
    elseif strcmp(type,'amosfsolve')
        kappa = besselInverseAmosFsolve(v,x);
    elseif strcmp(type,'matlabfsolve')
        kappa = besselInverseMatlabFsolve(v,x);
    elseif strcmp(type,'sra')
        kappa = besselInverseSra(v,x);
    end
end

function kappa = besselInverseFisher(x)
    % Approximation by Fisher for v=0
    %
    % Fisher, N. I. Statistical Analysis of Circular Data 
    % Cambridge University Press, 1995, eq (3.47)
    %
    % recommended if speed is more important than accuracy
    kappa = (x<0.53).*(2*x+x.^3+5*x.^5/6) + (x>=0.53).*(x<0.85).*(-0.4+1.39*x+0.43./(1-x)) + (x>=0.85)./(x.^3-4*x.^2+3*x);
end

function kappa = besselInverseAmosFsolve(v,x)
    % Use fsolve to calculate the inverse, use the approxmation by Amos for
    % the ratio of bessel functions, only works well until approx. kappa=1000
    f = @(t) besselratio(v,t) - x;
    start = 1;
    kappa = fsolve(f, start, optimset('Display', 'off', 'TolFun', 1e-40));
end

function kappa = besselInverseMatlabFsolve(v,x)
    % Use fsolve to calculate the inverse, use the Matlab for the ratio of
    % bessel functions, only works well until approx. kappa=1000
    f = @(t)  besseli(v+1,t,1)/besseli(v,t,1) - x;
    start = 1;
    kappa = fsolve(f, start, optimset('Display', 'off', 'TolFun', 1e-40));
end

function kappa = besselInverseSra(v,x)
    % Sra, S. 
    % A Short Note on Parameter Approximation for von Mises--Fisher Distributions: And a Fast Implementation of Is (x) 
    % Computational Statistics, Springer-Verlag, 2012, 27, 177-190
    d=2*v+2;
    % approximation by Banaerjee et al.
    kappa = x*(d-x^2)/(1-x^2);

    % newton approximation by Sra 2012
    % relies on matlab bessel functions
    % Attention: Sra defines A_d differently than we usually do!
    Ad = @(kappa) besseli(d/2, kappa, 1)/besseli(d/2-1, kappa, 1);
    newtonStep = @(kappa) kappa - (Ad(kappa)-x)/(1- Ad(kappa)^2 - (d-1)/kappa*Ad(kappa));
    for i=2:25
        kappaNew = newtonStep(kappa);
        if kappaNew == kappa
            break
        end
        kappa = kappaNew;
    end
end

function kappa = besselInverseSraAmos(x)
    % Combination of the methods by Sra and Amos
    % recommended for good if high accuracy is desired, also works for
    % large kappa
    
    % approximation by Banaerjee et al.
    kappa = x*(2-x^2)/(1-x^2);

    % newton approximation by Sra 2012
    Ad = @(kappa) besselratio(0,kappa); %use approximation by Amos to allow large values of kappa
    newtonStep = @(kappa) kappa - (Ad(kappa)-x)/(1- Ad(kappa)^2 - 1/kappa*Ad(kappa));
    for i=2:25
        kappaNew = newtonStep(kappa);
        if kappaNew == kappa
            break
        end
        kappa = kappaNew;
    end
end
