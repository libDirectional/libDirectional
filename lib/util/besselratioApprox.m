function res = besselratioApprox(v, x, type)
    % Computes the Ratio of Bessel-I Functions
    % The Ratio I_{v+1}(x)/I_v(x) is computed, where
    % I_v(x) is the Bessel-I function of order v evaluated at x.
    %
    % It uses approximations from
    % S.R. Jammalamadaka, A. SenGupta, Topics in Circular Statistics, World
    % Scientific Publ., New Jersey, 2001
    
    if nargin<3
        type='stienne13';
    end
    
    assert (v==0);

    if strcmp(type,'stienne13')
        %   as given in 
        %   Stienne, G., Reboul, S., Azmani, M., Choquel, J. and Benjelloun, M.
        %   A multi-temporal multi-sensor circular fusion filter
        %   Information Fusion, 2013
        % be aware that the resulting function is not continuous at x=0.6
        if x>=0.6
            %see Jammalamadaka 2001, page 290
            %see also Mardia 2000, page 350, eq. A13 with P=2
            res = 1-1/(2*x); 
        else
            %I'm not sure where this is stated, but it seems to correspond
            %to the approximation kappa = 1/sigma^2 when mapping VM and WN
            %however, this is commonly used for high concentrations
            res = exp(-1/(2*x));
        end
    elseif strcmp(type,'stienne12')    
        %   as given in 
        %   Stienne 12
        % be aware that the resulting function is not continuous at x=0.6
        if x>=0.6
            %see Jammalamadaka 2001, page 289, second order Taylor
            %expansion for I_p(kapap)
            %see also Mardia 2000, page 349, eq. A4
            res = (1-3/8/x-15/128/x^2)/(1+1/8/x + 9/128/x^2);
        else
            %see Jammalamadaka 2001, page 290
            %see also Mardia 2000, page 350, eq. A12 (with p=2)
            res = x/2; 
        end
    else
        error('invalid type')
    end
end
