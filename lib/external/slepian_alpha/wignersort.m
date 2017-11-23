function [c,oddperm,phasefix]=wignersort(l1,l2,l3,m1,m2,m3)
% [c,oddperm,phasefix]=WIGNERSORT(l1,l2,l3,m1,m2,m3)
%
% Computes index value for a given 3j symbol for storage, whether it
% legally exists or not! So care must be taken not to ask too much. 
%
% INPUT:
%
% l1,l2,l3     Top row of the Wigner 3j symbol [may be same-length vectors]
% m1,m2,m3     Bottom row of the Wigner 3j symbol [may be same-length vectors]
%
% OUTPUT:
%
% c            Index value of the Wigner 3j symbol
% oddperm      True if there were an odd number of permutations involved...
% phasefix     ...in which case the phase needs to be fixed by this
%  
% This indexing scheme for Wigner 3j symbols takes account of the Regge
% symmetries and is described in detail by Rasch and Yu (2003). For any
% valid set of l's and m's, this function calculates the corresponding index
% for storage and retrieval of precomputed Wigner 3j values. The oddperm
% output value keeps track of whether an odd number of permutations were
% involved in sorting the Regge square. If oddperm=1, the value of the 3j
% symbol should be multiplied by (-1)^(l1+l2+l3) before storage, or after
% retrieval.
%
% References: Rasch, J.; Yu, A.C.H. (2003) Efficient storage scheme for 
% precalculated Wigner 3j, 6j and Gaunt coefficients. SIAM Journal of 
% Scientific Computing. 25 (4), 1416-1428.
%
% Last modified by Kevin W. Lewis, 6/14/10
% Last modified by fjsimons-at-alum.mit.edu, 05/24/2011
% Last modified by Florian Pfaff for libDirectional, 10/04/2016
switch nargin
    case 0
        l1=0;l2=0;l3=0;
        m1=0;m2=0;m3=0;
    case 1
        l2=0;l3=0;
        m1=0;m2=0;m3=0;
    case 2
        l3=0;
        m1=0;m2=0;m3=0;
    case 3
        m1=0;m2=0;m3=0;
    case 4
        m2=0;m3=0;
    case 5
        m3=0;
end
phasefix=1;

oddperm=zeros(length(l1),1);
regge=zeros(length(l1),5);

l1vec=0;

% Loop over the l1's
for i=1:length(l1)
    % Make Regge square
    if length(l2)==1
      % Multiple l1's, constant l2,l3,m's
      R=[-l1(i)+l2+l3 l1(i)-l2+l3 l1(i)+l2-l3;...
          l1(i)-m1    l2-m2       l3-m3;...
          l1(i)+m1    l2+m2       l3+m3];
      l1vec=1;
    elseif length(l2)==length(l1) && length(l2)==length(l3) && ...
          length(m1)==1 && length(m2)==1  && length(m3)==1
      % There are a list of Wigner symbols to be retrieved
      R=[-l1(i)+l2(i)+l3(i) l1(i)-l2(i)+l3(i) l1(i)+l2(i)-l3(i);...
          l1(i)-m1          l2(i)-m2          l3(i)-m3;...
          l1(i)+m1          l2(i)+m2          l3(i)+m3];
    elseif length(l2)==length(l1) && length(l2)==length(l3) && ...
          length(m1)==length(m2) && length(m2)==length(m3)
      % There are a list of Wigner symbols to be retrieved
      R=[-l1(i)+l2(i)+l3(i) l1(i)-l2(i)+l3(i) l1(i)+l2(i)-l3(i);...
          l1(i)-m1(i)       l2(i)-m2(i)       l3(i)-m3(i);...
          l1(i)+m1(i)       l2(i)+m2(i)       l3(i)+m3(i)];
    else
      error('Either a vector l1 and all others scalar or all vectors')
    end
    temp=zeros(1,2);
    [~,minIndex]=min(R(:));
    [temp(1),temp(2)]=ind2sub([3 3],minIndex);
    
    R=circshift(R,-(temp-1));
    maxinds=find(R==max(R(:)));
    [temp(1),temp(2)]=ind2sub([3 3],...
        maxinds(find(maxinds==2|maxinds==3|maxinds==4|maxinds==7,1)));
                              
    % Reorder Regge square
    if temp(2)==1
      R=R';
      if temp(1)==3
        R(:,2:3)=fliplr(R(:,2:3));
        oddperm(i)=true;
      end
    else
      if temp(2)==3
        R(:,2:3)=fliplr(R(:,2:3));
        oddperm(i)=true;
      end
    end
    
    if R(3,2)<R(2,2)
      oddperm(i)=1-oddperm(i);
      R(2:3,:)=flipud(R(2:3,:));
    elseif R(3,2)==R(2,2) && R(3,3)<R(2,3);
      R(2:3,:)=flipud(R(2:3,:));
      oddperm(i)=1-oddperm(i);
    end
    regge(i,1:5)=[R(1,2) R(2,1) R(3,3) R(2,2) R(1,1)]; 
    if ~issorted(fliplr(regge(i,:)))
      disp('WHOA! something''s wrong here...')
    end
    L=regge(i,1);
    X=regge(i,2);
    T=regge(i,3);
    B=regge(i,4);
    S=regge(i,5);
    % Rasch and Yu eq. (2.13)
    c(i)=L*(24+L*(50+L*(35+L*(10+L))))/120+...
         X*(6+X*(11+X*(6+X)))/24+T*(2+T*(3+T))/6+...
         B*(B+1)/2+S+1;
end

% Make oddperm logical
oddperm=logical(oddperm);

% Here is the phase fix
if l1vec==0
  phasefix=(-1).^(l1(oddperm)+l2(oddperm)+l3(oddperm));
else
  phasefix=(-1).^(l1(oddperm)+l2+l3);
end
