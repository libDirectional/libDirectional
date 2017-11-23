function [s,C,S,L]=threej(l1,l2,l3,m1,m2,m3,L,meth,C,S)
% [s,C,S,L]=THREEJ(l1,l2,l3,m1,m2,m3,L,meth,C,S)
%
% Wigner 3j-symbol from a database precomputed by WIGNERCYCLE. Tries to
% find files which take account of the Regge symmetries as described by
% Rasch and Yu (2003) first.
%
% INPUT:
%
% l1,l2,l3     Top row of the Wigner 3j symbol [may be same-length vectors]
% m1,m2,m3     Bottom row of the Wigner 3j symbol [may be same-length vectors]
%              Their default is 0,0,0
% L            The bandwidth of the database [default: best available]
% meth         0 Uses sparse matrices [elegant, but slow]
%              1 Performs linear search on unsorted array [slow]
%              2 Performs binary search on presorted array [default]
% C,S          The column and element vectors resulting from a previous load
%
% OUTPUT:
%
% s            The (vector of) Wigner 3j symbols
% C,S          The column and element vectors good for the next load, OR
%              The entire Wigner 3j symmetric database and an additional empty 
% L            The L of the database that is loaded, should you not know
%
% EXAMPLE:
%
% threej('demo1') % Should return nothing if it all works
%
% SEE ALSO: WIGNER3JM, GAUNT, WIGNER0J, GUSEINOV, ZEROJ
%
% Last modified by kwlewis-at-princeton.edu, 06/14/2010
% Last modified by plattner-at-princeton.edu, 05/24/2011
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2011
% Last modified by Florian Pfaff for libDirectional, 11/04/2016

% Sometimes it returns zero if the database is corrupted, i.e. if it's
% not what it purports to be.
persistent wigner3jCell
if ~ischar(l1)
    switch nargin
        case 3
            m1=0;m2=0;m3=0;
            C=[];S=[];
            meth=2;
        case 4
            m2=0;m3=0;
            C=[];S=[];
            meth=2;
        case 5
            m3=0;
            C=[];S=[];
            meth=2;
        case {6,7}
            C=[];S=[];
            meth=2;
        case 8
            C=[];S=[];
        case 9
            S=[];
    end
  % Method
  % disp(sprintf('Using method %i',meth))

  
  % All saved values must be integers
  if sum(mod(l1,1)) || sum(mod(l2,1)) || sum(mod(l3,1)) ...
	sum(mod(m1,1)) || sum(mod(m2,1)) || sum(mod(m3,1))
    error('All degrees and orders must be integers for the database')
  end
  if any(m1+m2+m3)
    error('Sum of the orders must be zero')
  end

  if isempty(C) && isempty(S)
    % Got something to do
    
    % These are the modifications by Kevin W. Lewis
    % First check for saved files which take account of Regge symmetries
    if numel(wigner3jCell)<(max([l1(:)',l2(:)',l3(:)'])+1)||isempty(wigner3jCell{max([l1(:)',l2(:)',l3(:)'])+1}) % Do not look for files if already in memory
        try
          mfn=mfilename('fullpath');
          fnpl=fullfile(mfn(1:end-7),'WIGNER','WIGNER3J-*-sym.mat'); %modified for libDirectional 
          Els=ls2cell(fnpl);
        catch
          Els=[];
          EL=[];
        end
        for index=1:length(Els)
          EL(index)=str2num(rindeks(parse(Els{index},'-'),2));
        end
        EL=sort(EL);
        % Bandwidth of the database; keep this as low as possible
        % Need to provide for empties if not found
        fmax=find(max([l1(:)' l2(:)' l3(:)'])<=EL);
    else
        fmax=1;
        EL=(max([l1(:)',l2(:)',l3(:)']));
    end
    if ~isempty(fmax)
      if nargin<7
          L=EL(fmax(1));
      end
      % Check, identify and load database
      if any([l1(:)' l2(:)' l3(:)']>L)
        error('Insufficient bandwidth for database')
      end
      if numel(wigner3jCell)<(L+1)||isempty(wigner3jCell{L+1})
        mfn=mfilename('fullpath');
        fnpl=fullfile(mfn(1:end-7),'WIGNER',sprintf('WIGNER3J-%i-sym',L)); %modified for libDirectional 
        fprintf('Loading %s\n',fnpl); 
        load(fnpl) 
        wigner3jCell{L+1}=w3js; %#ok<NODEF>
      else
        w3js=wigner3jCell{L+1};
      end
      
      % These must be same-length vectors but two pairs of them can be
      % scalars 
      l1=l1(:); l2=l2(:); l3=l3(:);
      m1=m1(:); m2=m2(:); m3=m3(:);
      % But if ONE of them is a vector and the others aren't that's fine
      % too, as in the old version. Just take care of this one case now
      if length(l1)~=1 && length(l2)==1 && length(l2)==1
        l2=repmat(l2,length(l1),1);
        l3=repmat(l3,length(l1),1);
      elseif length(l2)~=1 && length(l1)==1 && length(l3)==1
        l1=repmat(l1,length(l2),1);
        l3=repmat(l3,length(l2),1);
      elseif length(l3)~=1 && length(l1)==1 && length(l2)==1
        l1=repmat(l1,length(l3),1);
        l2=repmat(l2,length(l3),1);
      end
      
      % Find where they're being kept
      [C,oddperm,phasefix]=wignersort(l1,l2,l3,m1,m2,m3);

      % Do the initial evaluation from the loaded variable
      s=full(w3js(C));
      % Fix the phase
      s(oddperm)=s(oddperm).*phasefix;
      % Now fix the triangle condition violations
      s=s.*triangle(l1,l2,l3);
      % Now fix the order violations
      s=s.*~[l1<abs(m1) | l2<abs(m2) | l3<abs(m3)];
      % There may be other violations that aren't being tested! The
      % database only stores the nonzero ones, for whichever reason 
      s=s(:)';
      
      % If you want this to return the database, you need to supply variables 
      C=w3js;
      S=[];
      return
    else
      % Try for files which do not take account of symmetries
      
      % What is the lowest of the available database bandwidths?
      try 
        mfn=mfilename('fullpath');
        fnpltmp=fullfile(mfn(1:end-7),'WIGNER','WIGNER3JCS-*-C'); %modified for libDirectional 
        Els=ls2cell(fnpltmp);
%         Els=ls2cell(fullfile(getenv('IFILES'),...
%                              'WIGNER','WIGNER3JCS-*-C'));
      catch
        Els=[];
        EL=[];
      end
      for index=1:length(Els)
        EL(index)=str2num(rindeks(parse(Els{index},'-'),2));
      end
      EL=sort(EL);
      % Bandwidth of the database; keep this as low as possible
      % Need to provide for empties if not found
      fmax=find(max([l1(:)' l2(:)' l3(:)'])<=EL);
      
      if ~isempty(fmax)
        if nargin<7
            L=EL(fmax(1));
        end
        % Check, identify and load database
        if any([l1(:)' l2(:)' l3(:)']>L)
          error('Insufficient bandwidth for database')
        end
%         fnpl1=sprintf('%s/WIGNER3JCS-%i-C',...
%                       fullfile(getenv('IFILES'),'WIGNER'),L);
%         fnpl2=sprintf('%s/WIGNER3JCS-%i-S',...
%                       fullfile(getenv('IFILES'),'WIGNER'),L);
        mfn=mfilename('fullpath');
        fnpl1=fullfile(mfn(1:end-7),'WIGNER',sprintf('WIGNER3JCS-%i-C',L));
        fnpl2=fullfile(mfn(1:end-7),'WIGNER',sprintf('WIGNER3JCS-%i-S',L));  
        if exist(fnpl1,'file')==2 && exist(fnpl2,'file')==2
          fmt1='uint64'; fmt2='float64';
          disp(sprintf('Loading %s',fnpl1))
          C=eval(sprintf('loadb(''%s'',''%s'')',fnpl1,fmt1));
          disp(sprintf('Loading %s',fnpl2))
          S=eval(sprintf('loadb(''%s'',''%s'')',fnpl2,fmt2));
          
          % Whatever happens, this better be sorted; check some entries
          randin=unique(ceil(rand(min(100,length(C)),1)*length(C)));
          if ~all(unique(C(randin))==C(randin))
            disp('Column arrays were not properly sorted')
            [C,j]=sort(C);
            writeb(C,fnpl1,fmt1)
            S=S(j);
            writeb(S,fnpl2,fmt2); clear j
            disp('Column arrays now sorted once and for all')
          end
        else
          if nargin<7
              L=-1;
          end
          % Precompute the database, this will fail though, on purpose
          wignercycle(L);
          % And have a go again
          s=threej(l1,l2,l3,m1,m2,m3,L,meth);
        end
      else
        if nargin<7
          L=-1;
        end
        % Precompute the database, this will fail though, on purpose
        wignercycle(L);
        % And have a go again
        s=threej(l1,l2,l3,m1,m2,m3,L,meth);
      end
    end
  else
    if isempty(L)
      error('If supplying vectors with coefficients must also specify bandwidth')
    end
    % Else have C and S from a previous load and do nothing, but check
    % There is NO check on whether the L supplied is appropriate for the
    % C & S combination that is indeed supplied
    if any([l1(:)' l2(:)' l3(:)']>L)
      error('Insufficient bandwidth for database')
    end
  end

  if isempty(S)
    % This is the symmetric storage scheme
    % These must be same-length vectors but two pairs of them can be
    % scalars 
    l1=l1(:); l2=l2(:); l3=l3(:);
    m1=m1(:); m2=m2(:); m3=m3(:);
    % But if ONE of them is a vector and the others aren't that's fine
    % too, as in the old version. Just take care of this one case now
    if length(l1)~=1 && length(l2)==1 && length(l2)==1
      l2=repmat(l2,length(l1),1);
      l3=repmat(l3,length(l1),1);
    elseif length(l2)~=1 && length(l1)==1 && length(l3)==1
      l1=repmat(l1,length(l2),1);
      l3=repmat(l3,length(l2),1);
    elseif length(l3)~=1 && length(l1)==1 && length(l2)==1
      l1=repmat(l1,length(l3),1);
      l2=repmat(l2,length(l3),1);
    end
    
    % Find where they're being kept
    [CC,oddperm,phasefix]=wignersort(l1,l2,l3,m1,m2,m3);

    % Do the initial evaluation from the loaded variable
    s=full(C(CC));
    % Fix the phase
    s(oddperm)=s(oddperm).*phasefix;
    % Now fix the triangle condition violations
    s=s.*triangle(l1,l2,l3);
    % Now fix the order violations
    s=s.*~[l1<abs(m1) | l2<abs(m2) | l3<abs(m3)];
    % There may be other violations that aren't being tested! The
    % database only stores the nonzero ones, for whichever reason 
    s=s(:)';
    
    return
  else
    % This is the old-school storage scheme
    % Now to the calculation proper if you hadn't stopped before
    % Find running index into this matrix
    if length(m1)~=1 || length(m2)~=1 || length(m3)~=1
      % Note there must be the same number of l1, l2, l3, m1, m2, m3
      % Initialize index vector
      rind=repmat(NaN,1,length(l1));
      for index=1:length(l1)
	rind(index)=addmabout(L,l1(index),m1(index))+...
	    (L+1)^2*(addmabout(L,l2(index),m2(index))-1)+...
	    (L+1)^4*(addmabout(L,l3(index),m3(index))-1);
      end
    else
      % So you can make a vector of l1 l2 l3 and have all zero m's
      % Turns out you can still have ONE that is a vector and all the
      % other ones that are singles, that would work too
      rind=addmabout(L,l1,m1)+(L+1)^2*(addmabout(L,l2,m2)-1)+...
	   (L+1)^4*(addmabout(L,l3,m3)-1);
    end
    % Initialize results vector
    s=repmat(NaN,1,length(rind));
    
    switch meth
     case 0
      % Turn it into the properly indexed sparse matrix
      % This step takes some of time but most of all, memory
      % HOLD ON - SHOULDN't THIS BE TO THE POWER SIX, RATHER?
      W=sparse(1,C,S,1,(L+1)^8);
      
      % Extract the Wigner 3j-symbol
      s=W(1,rind);
     case 1
      for index=1:length(rind)
	posi=find(C==rind(index));
	if ~isempty(posi)
	  s(index)=S(posi);
	else 
	  s(index)=0;
	end
      end
     case 2
      % Binary search algorithm on sorted arrays
      for index=1:length(rind)
	posi=binsearch(C,rind(index));
	if ~isempty(posi)
	  s(index)=S(posi);
	else 
	  s(index)=0;
	end
      end
     otherwise
      error('Specify valid method')
    end
  end
elseif strcmp(l1,'demo1')
  difer(full(threej([8 7 6 5 4],6,5,3,2,-5))-...
	     [sqrt(42/4199) 35*sqrt(2/138567) ...
	      7*sqrt(1/2431) sqrt(35/2431) sqrt(5/858)])
  % From Lai et al. based on Rotenberg's book
  difer(threej(16,16,6,0,0,0)-rotenberg([4 -2 2 0 -1 -1 1 1 0 -1 -1 -1],1))
  % For database of a bandwidth at least 3 
  % Illustrates vector capabilities
  difer(threej([0 0],[1 1],[1 1],[0 0],[1 0],[-1 0])-[sqrt(1/3) -sqrt(1/3)])
  difer(threej([2 0],[1 1],[1 1],0,0,0)-[sqrt(2/15) -sqrt(1/3)])
  % Illustrates straight scalar capabilities
  difer(threej(0,1,1,0,1,-1)-sqrt(1/3))
  difer(threej(0,1,1,0,0,0)--sqrt(1/3))
  difer(threej(1,1,1,0,1,-1)-sqrt(1/6))
  difer(threej(2,1,1,0,1,-1)-sqrt(1/30))
  difer(threej(2,1,1,0,0,0)-sqrt(2/15))
  difer(threej(2,2,3,-2,0,2)--sqrt(1/14))
  difer(threej(3,2,2,1,0,-1)--sqrt(1/35))
  difer(threej(3,2,3,-2,0,2)-0)
  difer(threej(3,2,3,3,0,-3)-(1/2)*sqrt(5/21))
  difer(threej(3,6,5,3,2,-5)-sqrt(1/1001))
  % For database of a bandwidth at least 6
  difer(threej(6,2,4,-4,1,3)-4*sqrt(1/429))
  % For database of a bandwidth at least 10
  difer(threej(10,10,12,9,3,-12)--(1/15)*sqrt(4199/9889))
  % For database of a bandwidth at least 20
  difer(threej(20,15,9,-3,2,1)+(311/2)*sqrt(115/1231496049))
  difer(threej(20,10,10,0,0,0)-indeks(wigner0j(20,10,10),'end'))
  difer(threej(20,8,12,-1,-3,4)- -(3/5)*sqrt(2261/1363783))
  % For database of a bandwidth at least 30
  difer(threej(30,21,12,-22,19,3)--(151/2)*sqrt(527/1277814153))
  difer(threej(27,19,8,-21,18,3)-(4/5)*sqrt(38/208131))
  % For database of a bandwidth at least 40
  difer(threej(0:40,10,13,0,0,0)-wigner0j(40,10,13))
  % For database of a bandwidth at least 60
  difer(threej(60,60,60,-60,40,20)-indeks(wigner3jm(60,60,60,-60,40,20),'end'))
  difer(threej(59,60,40,-20,10,10)-indeks(wigner3jm(59,60,40,-20,10,10),'end'))
  % Look at symmetry
  threej(repmat(20,1,19),repmat(15,1,19),repmat(9,1,19),...
	 repmat(-3,1,19),repmat(2,1,19),[-9:9])
end
