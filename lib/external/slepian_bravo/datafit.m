function [x,merr,mcov,chi2,L2err,rnk,Dw]=datafit(A,y,dcov,method,cnd)
% [x,merr,mcov,chi2,L2err,rnk,Dw]=DATAFIT(A,y,dcov,method,cnd)
%
% Solve linear system A*x=y for x, given a cov(y); y can be a matrix.
%
% Given m (orthogonal) functions evaluated at n points, where n>m, find the
% m best fit coefficients.
% 
% INPUT:
%
% A         matrix of size [M,N], N functions evaluated at M points.
% y         vector of M observations
% dcov      data coVARIANCE matrix (if vector, then diagonal)
% method    'svd'  By singular-value decomposition [default]
%           'gen'  By generalized inverse
%           'mp'   Moore-Penrose pseudo-inverse
% cnd       minimum acceptable condition number (ratio of smallest to
%           largest singular value) for SVD [default: 1e-6]
%
% OUTPUT:
%
% x        estimates of M best fit unknown coefficients
% merr     error in estimates of x
% mcov     covariance matrix of estimates x
% chi2     chi-square per degree of freedom
% err      L2-error: (y-A*x)'*(y-A*x)
% rnk      The rank of the matrix
% Dw       The eigenvalue spectrum
%
% Previously modified by fjsimons-at-alum.mit.edu, 07/29/2010
% Last modified by Florian Pfaff for libDirectional, 23/11/2017
%
% Using Numerical Recipes Chapters 2 and 15
% Using Aki and Richards, 1st Edition, Chapter 12.
printInfo=false;
if numel(y)==length(y)
  y=y(:);
end

% Define tolerance below which singular values will be considered equal
% to zero
switch nargin
    case 2
        dcov=NaN;
        method='svd';
        cnd=1e-6;
    case 3
        method='svd';
        cnd=1e-6;
    case 4
        cnd=1e-6;
end

[M,N]=size(A);
if N>M
  warning('More unknowns than equations. Add information');
end

flag=0;
if length(dcov)==1 
  if ~isnan(dcov)
    dcov=repmat(dcov,N,1);
  else
    flag=1;
  end
end
if numel(dcov)==length(dcov)
  dcov=diag(dcov);
end

% Calculate generalized inverse
switch method
 case 'svd' % Later substitute with PINV
  % Perform economy-size SVD 
  [U,w,V]=svd(A,0);
  % Threshold diagonal
  Dw=diag(w);
  if cnd~=1e-6
    fprintf('Condition number %3.3e\n',cnd)
  end
  invalid=Dw<cnd*max(Dw);
  Dw(invalid)=Inf;
  if M>N
    Gpi=V*diag(1./Dw)*U';
  else
    Gpi=V(:,1:M)*diag(1./Dw)*U';
  end
  % Get rank estimate
  rnk=sum(~invalid);
 case 'gen'
  Gpi=(A'*A)\A';
  rnk=[]; Dw=[];
 case 'mp'
  Gpi=pinv(A);
  rnk=[]; Dw=[];
end

% Calculate solution (y can be a matrix)
x=Gpi*y;

% L2 Error of the fit
err=(y-A*x);
L2err=err(:)'*err(:);
if nargout~=1 && printInfo
  fprintf('DATAFIT: L2 error norm %8.3e\n',L2err)
end

if numel(y)==length(y)
  % A posteriori statistics
  % Inverse data error matrix from data covariances
  if ~flag
    W=sqrt(inv(dcov));
    % Chi-square per degrees of freedom is the weighted error divided by the
    % number of parameters. This better be close to one.
    chi2=err'*W*err/N;
    % Model covariance matrix & standard error
    mcov=Gpi*dcov*Gpi';
  else
    mcov=Gpi*Gpi';
    chi2=err'*err/N;
  end
  merr=sqrt(diag(mcov));
  if nargout>3 && printInfo
    fprintf('Reduced chi-2 %8.3e\n',chi2)
    fprintf('Rank / Size %i / %i\n',rnk,length(Dw))
  end
else
  [merr,mcov,chi2]=deal([]);
end


