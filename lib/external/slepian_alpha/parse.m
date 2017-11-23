function stronk=parse(strink,varargin)
% stronk=parse(strink)
% stronk=parse(strink,sepor)
%
% Makes a string matrix from a long string with line breaks
% or with separator 'sepor' (which is a character used by FINDSTR).
%
% Application as in: files=parse(ls(['x200' '*']));
%
% Last modified by fjsimons-at-alum.mit.edu, 12.11.2004

if nargin==1
  ent=findstr(strink,sprintf('\n'))-1;
else
  ent=findstr(strink,varargin{1})-1;
end

if ~isempty(ent)
  % If the string does not end on the delimiter
  if ent(end)~=[length(strink)-1]
    ent=[ent length(strink)];
  end
  
  beg=[1 ent+2];
  
  stronk= ' ';
  for index=1:length(beg)-1
    stronk=str2mat(stronk,strink(beg(index):ent(index)));
  end
  stronk=stronk(2:end,:);
else %  Take the whole thing
  stronk=strink;
end

% Remove rows of blanks
delt=1:size(stronk,1);
for index=1:size(stronk,1)
  if all(abs(stronk(index,:))==32)
    delt(index)=0;
  end
end
stronk=stronk(~~delt,:);

