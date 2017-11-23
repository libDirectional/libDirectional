function [path,name,ext]=star69
% [path,name,ext]=STAR69
%
% Returns calling function name if exists.
%
% Last modified by fjsimons-at-alum.mit.edu, 04/13/2007

str=dbstack;

switch str2num(indeks(version,1))
 case 5
  if length(str) > 2
    str=str(3).name;
  else
    str='';
  end
  % A ridiculous holdout for version 5
 otherwise
  % Works for version 6 and version 7
  str=str(min(length(str),3)).name;
end

% This used to have a fourth output argument which is no longer supported 
[path,name,ext]=fileparts(str);

