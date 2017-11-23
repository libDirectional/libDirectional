function [w3j,j]=wigner3jm(L,l2,l3,m1,m2,m3)
% [w3j,j]=WIGNER3JM(L,l2,l3,m1,m2,m3)
%
% Calculates Wigner 3j symbols by recursion, for all values of j<=L
% allowed in the expression (L  l2 l3)
%                           (m1 m2 m3)
% There is no truncation at any bandwidth - they are all returned
% Note the selection rules:
% jmin = max(|l2-l3|, |m1|)
% jmax = l2 + l3
% m1 + m2 + m3 = 0. 
%
% INPUT:
%
% L           Maximum degree, bandwidth [default: that what's allowed]
% l2,l3       Other degrees in the symbol above, don't have to be integers
% m1,m2,m3    Orders in the symbol above, don't have to be integers, and
%             note that the defaults here are not zero!
%
% OUTPUT:
%
% w3j      The Wigner3j function
% j        The first degrees from 0 to L
%          With normalization check
%
% EXAMPLES:
% 
% wigner3jm('demo1') % Should return nothing if it all works
% wigner3jm('demo2') % Reproduces Table I of Schulten and Gordon
%
% See also: WIGNER0J, GUSEINOV, THREEJ, ZEROJ
%
% Last modified by fjsimons-at-alum.mit.edu, 06/15/2010
% Last modified by Florian Pfaff for libDirectional, 11/04/2016

% Note, if you truncate, like here, at the degree L, you're actually
% doing too much work; you could improve this.

% After Fortran (could you tell?) by Mark Wieczorek, who further notes:
% Returned values have a relative error less than ~1.d-8 when l2 and l3 are
% less than 103 (see below). In practice, this routine is probably usable up
% to 165.  This routine is based upon the stable non-linear recurrence
% relations of Luscombe and Luban (1998) for the "non classical" regions near
% jmin and jmax. For the classical region, the standard three term recursion
% relationship is used (Schulten and Gordon 1975). Note that this three term
% recursion can be unstable and can also lead to overflows. Thus the values
% are rescaled by a factor "scalef" whenever the absolute value of the 3j
% coefficient becomes greater than unity. Also, the direction of the iteration
% starts from low values of j to high values, but when abs(w3j(j+2)/w3j(j)) is
% less than one, the iteration will restart from high to low values. More
% efficient algorithms might be found for specific cases (for instance, when
% all m's are zero).

% Verification: 

% The results have been verified against this routine run in quadruple
% precision.  For 1.e7 acceptable random values of l2, l3, m2, and m3 between
% -200 and 200, the relative error was calculated only for those 3j
% coefficients that had an absolute value greater than 1.d-17 (values smaller
% than this are for all practical purposed zero, and can be heavily affected
% by machine roundoff errors or underflow). 853 combinations of parameters
% were found to have relative errors greater than 1.d-8. Here I list the
% minimum value of max(l2,l3) for different ranges of error, as well as the
% number of times this occured 1.d-7 < error <=1.d-8 = 103 # = 483 1.d-6 <
% error <= 1.d-7 = 116 # = 240 1.d-5 < error <= 1.d-6 = 165 # = 93 1.d-4 <
% error <= 1.d-5 = 167 # = 36

% Many times (maybe always), the large relative errors occur when the 3j
% coefficient changes sign and is close to zero. (I.e., adjacent values are
% about 10.e7 times greater in magnitude.) Thus, if one does not need to know
% highly accurate values of the 3j coefficients when they are almost zero
% (i.e., ~1.d-10) then this routine is probably usable up to about 160.

% These results have also been verified for parameter values less than 100
% using a code based on the algorithm of de Blanc (1987), which was originally
% coded by Olav van Genabeek, and modified by M. Fang (note that this code was
% run in quadruple precision, and only calculates one coefficient for each
% call. I also have no idea if this code was verified.) Maximum relative
% errors in this case were less than 1.d-8 for a large number of values
% (again, only 3j coefficients greater than 1.d-17 were considered here).  The
% biggest improvement that could be made in this routine is to determine when
% one should stop iterating in the forward direction, and start iterating from
% high to low values.

if ~ischar(L)
  switch nargin
      case 0
        l2=6;l3=5;
        L=l2+l3;
        m1=3;m2=2;m3=-5;
      case 1
        l2=6;l3=5;
        m1=3;m2=2;m3=-5;
      case 2
        l3=5;
        m1=3;m2=2;m3=-5;  
      case 3
        m1=3;m2=2;m3=-5;  
      case 4
        m2=2;m3=-5;
      case 5
        m3=-5; 
  end
  flag1=0;
  flag2=0;

  % Didn't realize this factor was wrong until 10/02/2006
  % But - Luscombe and Luban imply that the three-term recursion in the
  % classical, oscillatory region rarely suffers from the overflows - all
  % should be probably well
  scalef=1000;

  jmin=max(abs(l2-l3),abs(m1));
  jmax=l2+l3;
  jnum=jmax-jmin+1;

  % Initialize
  w3j=zeros(1,jnum);
  if abs(m2)>l2 || abs(m3)>l3
    w3j=zeros(L+1,1)';
    j=0:L; return
  elseif m1+m2+m3~=0 
    w3j=zeros(L+1,1)';
    j=0:L; return
  elseif jmax<jmin
    w3j=zeros(L+1,1)';
    j=0:L; return
  end

  % Only one value is allowed
  if jnum==1 
    w3j=1/sqrt(2*jmin+1);
    if (w3j<0 && (-1)^(l2-l3+m2+m3)>0) || ...
	  (w3j>0 && (-1)^(l2-l3+m2+m3)<0)
      w3j=-w3j;
    end
    [w3j,j]=output(jmin,l2,l3,L,w3j); return
  end

  % Calculate lower non-classical values for [jmin, jn]. If the second
  % term can not be calculated because the recursion relationsips give
  % rise to a 1/0, then set flag1 to 1. If all m's are zero, then this
  % is not a problem as all  odd terms must be zero.

 
  [rs,wu,wl]=deal(0);

  warning off
  rs(1)=-x(jmin,l2,l3,m1)/y(jmin,l2,l3,m1,m2,m3);
  warning on

  if m1==0 && m2==0 && m3==0
    wl(1)=1;
    wl(2)=0;
    jn=jmin+1;
  elseif y(jmin,l2,l3,m1,m2,m3)==0
    if x(jmin,l2,l3,m1)==0
      flag1=1;
      jn=jmin;
    else
      wl(1)=1;
      wl(2)=0;
      jn=jmin+1;
    end
  elseif rs(1)<=0
    wl(1)=1;
    wl(2)=-y(jmin,l2,l3,m1,m2,m3)/x(jmin,l2,l3,m1);
    jn=jmin+1;
  else
    jn=jmax;
    for j=jmin+1:jmax-1
      denom=y(j,l2,l3,m1,m2,m3)+z(j,l2,l3,m1)*rs(j-jmin);
      warning off
      rs(j-jmin+1)=-x(j,l2,l3,m1)/denom;
      warning on
      if (rs(j-jmin+1)>1 || rs(j-jmin+1) <= 0 || denom==0) 
	jn=j-1;
	break
      end
    end
    
    wl(jn-jmin+1)=1;
    
    for k=1:jn-jmin
      wl(jn-k-jmin+1)=wl(jn-k-jmin+2)*rs(jn-k-jmin+1);
    end
    if (jn==jmin) 					
      wl(2)=-y(jmin,l2,l3,m1,m2,m3)/x(jmin,l2,l3,m1);
      jn=jmin+1;
    end
  end

  if jn==jmax
    w3j(1:jnum)=wl(1:jnum);
    w3j=normw3j(w3j,jmin,jmax,jnum);
    w3j=fixsign(w3j,jnum,l2,l3,m2,m3);
    [w3j,j]=output(jmin,l2,l3,L,w3j); return
  end

  % Calculate upper non-classical values for  [jp, jmax].
  % If the second last term can not be calculated because the
  % recursion relations give a 1/0, then set flag2 to 1. 

  warning off
  rs(jnum)=-z(jmax,l2,l3,m1)/y(jmax,l2,l3,m1,m2,m3);
  warning on

  if m1==0 && m2==0 && m3==0
    wu(jnum)=1;
    wu(jmax-jmin)=0;
    jp=jmax-1;
  elseif y(jmax,l2,l3,m1,m2,m3)==0
    if z(jmax,l2,l3,m1)==0
      flag2=1;
      jp=jmax;
    else
      wu(jnum)=1;
      wu(jmax-jmin)=-y(jmax,l2,l3,m1,m2,m3)/z(jmax,l2,l3,m1);
      jp=jmax-1;
    end
  elseif rs(jnum)<=0 
    wu(jnum)=1;
    wu(jmax-jmin)=-y(jmax,l2,l3,m1,m2,m3)/z(jmax,l2,l3,m1);
    jp=jmax-1;
  else
    jp=jmin;
    for j=jmax-1:-1:jn
      % This appears to be Luscombe and Luban's Eq. (2)
      denom=y(j,l2,l3,m1,m2,m3)+x(j,l2,l3,m1)*rs(j-jmin+2);
      warning off
      rs(j-jmin+1)=-z(j,l2,l3,m1)/denom;
      warning on
      if (rs(j-jmin+1)>1 || rs(j-jmin+1) <= 0 || denom==0)
	jp=j+1;
	break
      end
    end	
    wu(jp-jmin+1)=1;
    for k=1:jmax-jp
      wu(jp+k-jmin+1)=wu(jp+k-jmin)*rs(jp+k-jmin+1);
    end
    
    if jp==jmax
      wu(jmax-jmin)=-y(jmax,l2,l3,m1,m2,m3)/z(jmax,l2,l3,m1);
      jp=jmax-1;
    end
  end

  % Calculate classical terms for [jn+1, jp-1] using standard three term
  % rercusion relationship. Start from both jn and jp and stop at the
  % midpoint. If flag1 is set, then perform the recursion solely from
  % high to low values. If flag2 is set, then perform the recursion
  % solely from low to high.

  if flag1==0
    % I think Fortran rounds like this
    jmid=round((jn+jp)/2);
    
    for j=jn:jmid-1
      wl(j-jmin+2)=-(z(j,l2,l3,m1)*wl(j-jmin)+...
			     y(j,l2,l3,m1,m2,m3)* ...
			     wl(j-jmin+1))/x(j,l2,l3,m1);  
      if abs(wl(j-jmin+2))>1
	wl(1:j-jmin+2)=...
	    wl(1:j-jmin+2)/scalef;
      end
      
      if abs(wl(j-jmin+2)/wl(j-jmin))<1 && ...
	    wl(j-jmin+2)~=0 
	jmid=j+1;				
	break
      end
    end
    wnmid=wl(jmid-jmin+1);

    warning off
    if abs(wnmid/wl(jmid-jmin))<1e-6 && wl(jmid-jmin)~=0
      wnmid=wl(jmid-jmin);
      jmid=jmid-1;
    end
    warning on
    
    for j=jp:-1:jmid+1
      wu(j-jmin)=-(x(j,l2,l3,m1)*wu(j-jmin+2)+...
			     y(j,l2,l3,m1,m2,m3)* ...
			     wu(j-jmin+1))/z(j,l2,l3,m1);  
      if abs(wu(j-jmin))>1
	wu(j-jmin:jnum)=...
	    wu(j-jmin:jnum)/scalef;
      end
    end

    wpmid=wu(jmid-jmin+1);
    
    if jmid==jmax
      w3j(1:jnum)=wl(1:jnum);
    elseif jmid==jmin
      w3j(1:jnum)=wu(1:jnum);
    else
      w3j(1:jmid-jmin+1)=wl(1:jmid-jmin+1)*wpmid/wnmid; 
      w3j(jmid-jmin+2:jnum)=...
	  wu(jmid-jmin+2:jnum);
    end
    
  elseif flag1==1 && flag2==0
    for j=jp:-1:jmin+1
      wu(j-jmin)=-(x(j,l2,l3,m1)*wu(j-jmin+2)+...
			     y(j,l2,l3,m1,m2,m3)* ...
			     wu(j-jmin+1))/z(j,l2,l3,m1);    
      if abs(wu(j-jmin))>1
	wu(j-jmin:jnum)=...
	    wu(j-jmin:jnum)/scalef;
      end
    end
    
    w3j(1:jnum)=wu(1:jnum);
    
  elseif flag2==1 && flag1==0
    
    for j=jn:jp-1
      wl(j-jmin+2)=-(z(j,l2,l3,m1)*wl(j-jmin)+...
			     y(j,l2,l3,m1,m2,m3)* ...
			     wl(j-jmin+1))/x(j,l2,l3,m1);  
      if abs(wl(j-jmin+2))>1
	wl(1:j-jmin+2)=...
	    wl(1:j-jmin+2)/scalef;
      end
    end
    
    w3j(1:jnum)=wl(1:jnum);
    
  elseif flag1==1 && flag2==1
    error('Can not calculate function for input values')
  end

  w3j=normw3j(w3j,jmin,jmax,jnum);
  w3j=fixsign(w3j,jnum,l2,l3,m2,m3);

  % Output: give output in all degrees from zero to the bandwidth
  if L<=l2+l3
    % Truncate
    w3j=[repmat(0,1,jmin) w3j(1:L-jmin+1)];
  else
    % Append
    w3j=[repmat(0,1,jmin) w3j repmat(0,1,L-l2-l3)];
  end
  j=0:L;

  % Perform normalization check
  if L==l2+l3
    norma=sum((2*[0:L]+1).*w3j.^2);
    difer(norma-1,[],[],NaN)
  end

elseif strcmp(L,'demo1')
  difer(wigner3jm(20,10,10,0,0,0)-wigner0j(20,10,10))
  difer(wigner3jm(400,200,200,0,0,0)-wigner0j(400,200,200))
  difer(wigner3jm(40,10,13,0,0,0)-wigner0j(40,10,13))
  difer(indeks(wigner3jm(8,6,5,3,2,-5),'end')-sqrt(42/4199))
  difer(indeks(wigner3jm(8,6,5,3,2,-5),'end-1')-35*sqrt(2/138567))
  difer(indeks(wigner3jm(8,6,5,3,2,-5),'end-2')-7*sqrt(1/2431))
  difer(indeks(wigner3jm(8,6,5,3,2,-5),'end-3')-sqrt(35/2431))
  difer(indeks(wigner3jm(8,6,5,3,2,-5),'end-4')-sqrt(5/858))
  difer(indeks(wigner3jm(3,2,2,1,0,-1),'end')--sqrt(1/35))
  difer(indeks(wigner3jm(0,1,1,0,1,-1),'end')-sqrt(1/3))
  difer(indeks(wigner3jm(1,1,1,0,1,-1),'end')-sqrt(1/6))
  difer(indeks(wigner3jm(2,1,1,0,1,-1),'end')-sqrt(1/30))
  difer(indeks(wigner3jm(3,2,3,-2,0,2),'end')-0)
  difer(indeks(wigner3jm(6,2,4,-4,1,3),'end')-4*sqrt(1/429))
  difer(indeks(wigner3jm(0,1,1,0,0,0),'end')--sqrt(1/3))
  difer(indeks(wigner3jm(2,1,1,0,0,0),'end')-sqrt(2/15))
  difer(indeks(wigner3jm(3,2,3,3,0,-3),'end')-(1/2)*sqrt(5/21))
  difer(indeks(wigner3jm(2,2,3,-2,0,2),'end')--sqrt(1/14))
  difer(indeks(wigner3jm(8,6,5,3,2,-5),'end-5')-sqrt(1/1001))
  difer(indeks(wigner3jm(20,15,9,-3,2,1),'end')+(311/2)*sqrt(115/1231496049))
  difer(indeks(wigner3jm(10,10,12,9,3,-12),'end')--(1/15)*sqrt(4199/9889))
elseif strcmp(L,'demo2')
  [w,L]=wigner3jm([],9/2,7/2,1,-7/2,5/2);
  disp(sprintf('------------------------------'))
  disp(sprintf('L1  Values of 3j coefficients'))
  disp(sprintf('------------------------------'))
  disp(sprintf('%2.2i %23.16i\n',[L(2:9)' w(2:9)']'))
  disp(sprintf('------------------------------'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [o,p]=output(jmin,l2,l3,L,w3j)
% Output: give output in all degrees from zero to the bandwidth
% Perform normalization check before possible completion or truncation
norma=sum((2*[jmin:l2+l3]+1).*w3j.^2);
difer(norma-1,[],[],NaN)
if L<=l2+l3
  % Truncate...
  o=[repmat(0,1,jmin) w3j(1:L-jmin+1)];
else
  % Append
  o=[repmat(0,1,jmin) w3j repmat(0,1,L-l2-l3)];
end
p=0:L;

% The following functions are straight from Table I of Luscombe and Luban 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o=x(j,l2,l3,m1)	
o=j*a(j+1,l2,l3,m1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o=y(j,l2,l3,m1,m2,m3)
% This is the Y-function in Table 1 of Luscombe and Luban
o=-(2*j+1)*(m1*(l2*(l2+1)-l3*(l3+1))-(m3-m2)*j*(j+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o=z(j,l2,l3,m1)
o=(j+1)*a(j,l2,l3,m1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o=a(j,l2,l3,m1)
o=sqrt((j^2-(l2-l3)^2)*((l2+l3+1)^2-j^2)*(j^2-m1^2));
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o=normw3j(w3j,jmin,jmax,jnum)
normo=0;
for j=jmin:jmax
  normo=normo+(2*j+1)*w3j(j-jmin+1)^2;
end
o(1:jnum)=w3j(1:jnum)/sqrt(normo);
			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o=fixsign(w3j,jnum,l2,l3,m2,m3)
if ((w3j(jnum)<0 & (-1)^(l2-l3+m2+m3)>0) |...
    (w3j(jnum)>0 & (-1)^(l2-l3+m2+m3)<0)) 
  o(1:jnum)=-w3j(1:jnum);
else
  o(1:jnum)=w3j(1:jnum);
end
		
