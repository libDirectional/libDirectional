function res = circVMcdf(T,VK)
%circVMcdf cumulative Von-Mises distribution VM(0,k)
%   res = circVMcdf(T,VK)
%       T - angles at which to compute CDF
%       VK - kappa value for distribution
%
%  Directly converted from Fortran code published in
%
%       Algorithm 518: Incomplete Bessel Function I0.
%                      The Von Mises Distribution [S14]
%       ACM Transactions on Mathematical Software (TOMS)
%       Volume 3 ,  Issue 3  (September 1977)
%       Pages: 279 - 284
%       Author: Geoffrey W. Hill
%
%  (A BibTeX citation is in a comment at the end of this file)
%
%  By Shai Revzen, Berkeley 2006
%  Modified by Gerhard Kurz, KIT 2015


    % 8 digit accuracy
    % CK = 10.5;
    % 12 digit accuracy
    CK = 50;

    if length(VK) ~= 1
        error('circ:mustBeScalar','kurtosis must be a scalar')
    end
    %T = T-mu;
    Z = VK;
    U = mod(T+pi,2*pi);
    Y = U-pi;
    if Z>CK
        res = largeVK(Y, Z);
    elseif Z<=0
        res = (U*0.5)/pi;
    else
        V = smallVK(Y, Z);
        res = (U*0.5+V)/pi;
    end
    res(res<0)=0;
    res(res>1)=1;
end

function V = smallVK( Y, Z )
    % 8 digit accuracy
    %A1 = 12; A2 = 0.8; A3 = 8.0; A4 = 1.0;   
    % 12 digit accuracy    
    A1 = 28; A2 = 0.5; A3 = 100.0; A4 = 5.0;   
    
    IP = Z*A2 - A3/(Z+A4) + A1;
    P = round(IP);
    S = sin(Y);
    C = cos(Y);
    Y = P*Y;
    SN = sin(Y);
    CN = cos(Y);
    R = 0.0;
    V = 0.0;
    Z = 2.0/Z;
    for N=2:round(IP)
        P = P - 1.0;
        Y = SN;
        SN = SN.*C - CN.*S;
        CN = CN.*C + Y.*S;
        R = 1.0./(P*Z+R);
        V = (SN./P+V)*R;
    end
end

function res = largeVK( Y, Z )
    % 8 digit accuracy
    % C1 = 56;
    % 12 digit accuracy
    C1 = 50.1;

    C = 24.0 * Z;
    V = C - C1;
    R=sqrt((54.0/(347.0/V+26.0-C)-6.0+C)/12.0);
    Z=sin(Y*0.5)*R;
    S=Z.*Z*2.0;
    V = V - S + 3.0;
    Y = (C-S-S-16.0)/3.0;
    Y = ((S+1.75)*S+83.5)/V - Y;
    res=erf(Z-S/(Y.*Y).*Z)*0.5+0.5;
end

% BibTeX:
%  @article{355753,
%   author = {Geoffrey W. Hill},
%   title = {Algorithm 518: Incomplete Bessel Function I0.
%            The Von Mises Distribution [S14]},
%   journal = {ACM Trans. Math. Softw.},
%   volume = {3},
%   number = {3},
%   year = {1977},
%   issn = {0098-3500},
%   pages = {279--284},
%   doi = {http://doi.acm.org/10.1145/355744.355753},
%   publisher = {ACM Press},
%   address = {New York, NY, USA},
%  }