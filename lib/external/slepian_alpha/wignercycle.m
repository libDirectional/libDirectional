function wignercycle(L,emnz,simi)
% WIGNERCYCLE(L,emnz,simi)
%
% If WIGNER0J is the wheel, this is the bicycle, and ZEROJ the rider.
% If WIGNER3JM is the wheel, this is the bicycle, and THREEJ the rider.
% If WIGNER6J is the wheel, this is the bicycle, and SIXJ the rider.
%
% INPUT:
%
% L      Bandwidth of the data base that will be produced
% emnz   1 Do this for all possible orders [default]
%        0 Only do this for zero orders
%        6 Rather calculate the wigner 6j coefficients
% simi   0 Standard, very cumbersome storage scheme
%        1 Smart storage scheme a la Rasch and Yu (2003) for 3j only
% xver   1 Perform excessive verification (a wise choice)
%        0 Don't (an unwise choice)
%
% OUTPUT:
%
% Only to file, suitable to be read in by ZEROJ, THREEJ, SIXJ
%
% Creates a set of huge, but sparse matrices of Wigner 3j/6j 
% coefficients... up to a maximum degree L, in blocky ADDMOUT/ADDMABOUT
% ordering (for the 3j's). Saves them in a vector with column numbers and
% with data. For 3j symbols, when simi=1, uses the indexing scheme
% described by Rasch and Yu (2003), to take account of natural
% symmetries.  In this case, the values are stored as a sparse array in
% .mat format.
%
% All degrees/orders INTEGERS... not required by WIGNER3JM/WIGNER6J.
%
% See ZEROJ, THREEJ, SIXJ, WIGNERCYCLESYM, WIGNERSORT
% 
% Last modified by klewis-at-princeton.edu, 06/14/2010
% Last modified by fjsimons-at-alum.mit.edu, 05/24/2011
% Last modified by Florian Pfaff for libDirectional, 11/04/2016

t=cputime;

timewaster=0;
switch nargin
    case 0
        L=30;
        emnz=1;
        simi=1;
    case 1
        emnz=1;
        simi=1;
    case 2
        simi=1;
end
xver=0;

% Files to save the results in: one float, one integer
switch emnz
 case 1
  if simi==0
    mfn=mfilename('fullpath');
    fnpl1=fullfile(mfn(1:end-12),'WIGNER',sprintf('WIGNER3JCS-%i-C',L));
    fnpl2=fullfile(mfn(1:end-12),'WIGNER',sprintf('WIGNER3JCS-%i-S',L));
%     fnpl1=sprintf('%s/WIGNER3JCS-%i-C',...
%                   fullfile(getenv('IFILES'),'WIGNER'),L);
%     fnpl2=sprintf('%s/WIGNER3JCS-%i-S',...
%                   fullfile(getenv('IFILES'),'WIGNER'),L);
  elseif simi==1
    mfn=mfilename('fullpath');
    fnpl=fullfile(mfn(1:end-12),'WIGNER',sprintf('WIGNER3J-%i-sym',L));
%     fnpl=sprintf('%s/WIGNER3J-%i-sym',...
% 		 fullfile(getenv('IFILES'),'WIGNER'),L);
  end
  % Get orders and degrees in block form
  [EM,EL,mz,blkm,dblk]=addmout(L);
  EM=EM(blkm); EL=EL(blkm);
 case 0
   mfn=mfilename('fullpath');
   fnpl1=fullfile(mfn(1:end-12),'WIGNER',sprintf('WIGNER0JCS-%i-C',L));
   fnpl2=fullfile(mfn(1:end-12),'WIGNER',sprintf('WIGNER0JCS-%i-S',L));
%    fnpl1=sprintf('%s/WIGNER0JCS-%i-C',...
% 		 fullfile(getenv('IFILES'),'WIGNER'),L);
%    fnpl2=sprintf('%s/WIGNER0JCS-%i-S',...
% 		 fullfile(getenv('IFILES'),'WIGNER'),L);
   % Restrict to the zero orders only
   EM=repmat(0,1,L+1)';
   EL=[0:L]';
 case 6
   mfn=mfilename('fullpath');
   fnpl1=fullfile(mfn(1:end-12),'WIGNER',sprintf('WIGNER6JCS-%i-C',L));
   fnpl2=fullfile(mfn(1:end-12),'WIGNER',sprintf('WIGNER6JCS-%i-S',L));  
%    fnpl1=sprintf('%s/WIGNER6JCS-%i-C',...
% 		 fullfile(getenv('IFILES'),'WIGNER'),L);
%    fnpl2=sprintf('%s/WIGNER6JCS-%i-S',...
% 		 fullfile(getenv('IFILES'),'WIGNER'),L);
 otherwise
  error('Specify valid case')
end

% For format, see INTMAX and REALMAX
if simi==0
  fid1=fopen(fnpl1,'w'); fmt1='uint64';
  fid2=fopen(fnpl2,'w'); fmt2='float64';
end
  
if emnz~=6
  % Loop over m2, then over m3, then pick m1
  for ind2=1:length(EM)
    l2=EL(ind2); m2=EM(ind2);
    % Only pick an l3 if its m3 implies an abs(m1) that is
    % not bigger than l1 
    EL3AL=EL(abs(-m2-EM)<=L);
    EM3AL=EM(abs(-m2-EM)<=L);
    for ind3=1:length(EL3AL)
      l3=EL3AL(ind3); m3=EM3AL(ind3);
      m1=-m2-m3;
      
      % if timewaster==1
      %  % Note that display time is expensive!
      %  clc
      %  disp(sprintf('(%+3i %+3i %+3i)',L,l2,l3))
      %  disp(sprintf('(%+3i %+3i %+3i)',m1,m2,m3))
      % end
      
      % Do the calculation
      switch emnz
       case 1
	[w3j,l1]=wigner3jm(L,l2,l3,m1,m2,m3);
        if simi==1
          % So here is a modification by Kevin W. Lewis
          [c,oddperms,phasefix]=wignersort(l1,l2,l3,m1,m2,m3);
          % See also THREEJ
          w3j(oddperms)=w3j(oddperms).*phasefix;
          jmin=max(abs(l2-l3),abs(m1))+1;
          jmax=min(l2+l3+1,length(c));
          if exist('w3js','var') && length(w3js)>=max(c)
            if sum(w3js(c(jmin:jmax))~=0 & ...
                   w3js(c(jmin:jmax))-w3j(jmin:jmax)'>1e-10)
              disp(['Differences = ' num2str(w3js(c)-w3j)])
            end
          end
          if ~exist('w3js','var')
            % Rasch and Yu (2003) eq. (2.14)
            w3js=sparse(...
                2*L*(274+2*L*(225+2*L*(85+2*L*(15+2*L))))/120+1,1);
          end
          w3js(c(jmin:jmax),1)=w3j(jmin:jmax);
          if ~exist('ind','var')
            ind=0; 
          end
          totalindex(c(jmin:jmax))=true;
          ind=ind+jmax-jmin;
          % disp(num2str([ind]))
          if ind2==length(EM) && ind3==length(EL3AL)
            save(fnpl,'w3js');
          end
        elseif simi==0
          % Find the possible and non-trivial ones in this
          % If there are any, save them in a big array
          nontriv=w3j~=0;
          if sum(nontriv)>0
            % Now construct running indices for the non-zeroes only
            % This is the 3-D equivalent for SUB2IND which you may use
            C=            addmabout(L,l1(nontriv),m1)...
                +(L+1)^2*(addmabout(L,l2         ,m2)-1)...
                +(L+1)^4*(addmabout(L,l3         ,m3)-1);
            % Write this out on the fly
            % The column number, "C"
            fwrite(fid1,C,fmt1); 
            % The actual element, "S"
            fwrite(fid2,w3j(nontriv),fmt2); 
          end
        end
       case 0
	% Perform the orthogonality check or not
	[w3j,l1]=wigner0j(L,l2,l3,xver);
	if xver==1
	  % Check, and check again with the different routine
	  difer(w3j-wigner3jm(L,l2,l3,0,0,0));
	end
        % Find the possible and non-trivial ones in this
        % If there are any, save them in a big array
        nontriv=w3j~=0;
        if sum(nontriv)>0
          % Now construct running indices for the non-zeroes only
          % This is the 3-D equivalent for SUB2IND which you may use
          C=          addmabout(L,l1(nontriv),m1)...
              +(L+1)^2*(addmabout(L,l2,m2)-1)...
              +(L+1)^4*(addmabout(L,l3,m3)-1);
          % Write this out on the fly
          fwrite(fid1,C,fmt1); % The column number, "C"
          fwrite(fid2,w3j(nontriv),fmt2); % The actual element, "S"
        end
      end
    end
  end
elseif emnz==6
  % Here we can avoid loops since the arrays are much smaller, and then
  % we can use the selection rules directly; may use NDGRID for this
  l2=       gamini([0:L],(L+1)^4)';
  l3=repmat(gamini([0:L],(L+1)^3)',(L+1)^1,1);
  l4=repmat(gamini([0:L],(L+1)^2)',(L+1)^2,1);
  l5=repmat(gamini([0:L],(L+1)^1)',(L+1)^3,1);
  l6=repmat([0:L]',                (L+1)^4,1);

  % Selection rules - those of the 6j contain those of the 6j
  selx=triangle(l4,l2,l6) & triangle(l4,l5,l3);
  
  % Apply the selection rules
  l2=l2(selx); l3=l3(selx);
  l4=l4(selx); l5=l5(selx); l6=l6(selx);
  
  for ix=1:length(l2)
    % Calculate the Wigner 6j symbols with normalization check
    [w6j,l1,norma]=wigner6j(l2(ix),l3(ix),l4(ix),l5(ix),l6(ix));
    
    % Now remember that l1 returned may exceed L of database, cut
    if l1(end)>L
      w6j=w6j(1:L+1);
      l1=l1(1:L+1);
    end

    % disp(sprintf('Passed normalization with %8.3e',norma-1))
    nontriv=w6j~=0;
    if sum(nontriv)>0
      % Now construct running indices for the non-zeroes only
      % This is the 3-D equivalent for SUB2IND which you may use
      C=1+l1(nontriv)+l2(ix)*(L+1)^1+...
	l3(ix)*(L+1)^2+...
	l4(ix)*(L+1)^3+...
	l5(ix)*(L+1)^4+...
	l6(ix)*(L+1)^5; 
      % Write this out on the fly
      fwrite(fid1,C,fmt1); % The column number, "C"
      fwrite(fid2,w6j(nontriv),fmt2); % The actual element, "S"
    end
  end
end

% Close files
if simi==0
  fclose(fid1);
  fclose(fid2);
end

disp(sprintf('WIGNERCYCLE: Elapsed CPU time %8.3f',cputime-t))

