function [ res, normObjFun ] = gaussNewtonMLE( om, start )
% Computes maximum likelihood estimate for Bingham 
% distribution.
%   om - Eigenvalues of Covariance matrix of samples.
%   start - initial value

    m = length(om);
    normConstDeriv = zeros(m,1);
    normConstHessian = zeros(m,m);
    normObjFun=0;
    eps=1e-10;
  
    if nargin > 1
        res = -start;
    else
        %res = zeros(m,1)-eps;
        res = 1./om;
    end   
    
    for k=0:100
        res = (res-min(res));
        normConst = numericalSaddlepointWithDerivatives(res);
        normConst = normConst(3);

        for i=1:m
            [tmpNC, tmpDeriv] = ...
                numericalSaddlepointWithDerivatives(res([1:i i i:m]));

            normConstDeriv(i) = -tmpNC(3)/(2*pi);

            normConstHessian(i,1:i-1) = -tmpDeriv(3,1:i-1)/(2*pi);
            normConstHessian(i,i) = -3*tmpDeriv(3,i)/(2*pi); % Total Derivative.
            normConstHessian(i,i+1:m) = -tmpDeriv(3,(i+3:m+2))/(2*pi);
        end       
        
        objFun = normConstDeriv/normConst - om; 
        objFunJacobian =  ...
            (normConstHessian * normConst - normConstDeriv*normConstDeriv') ...
            / normConst^2;       
        
        % Gauss Newton step.
        %initval
        %normConst
        %normConstDeriv
        %normConstHessian
        %goalFun
        %goalFunJacobian
        oldNorm = normObjFun;
        normObjFun=norm(objFun);
        res = res - pinv(objFunJacobian'*objFunJacobian)*objFunJacobian' * objFun;
        
        if norm(oldNorm - normObjFun) < eps
            %k
            break;
        end
    end
   
	res = -(res-min(res));
end