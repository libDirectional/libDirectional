function cls=ls2cell(ddir,fullpath)
% cls=LS2CELL(ddir,fullpath)
%
% Makes cell array with the names in a directory
%
% INPUT:
%
% ddir      A directory name string, (some) wildcards possible
% fullpath  1 if you would like the full file path output [default: 0]
%             in this case include the filespe in the directory name
%           0 returns just the file names
%
% OUTPUT:
%
% cls     Cell structure with file names, with full file path if you want 
%
% NOTE: Remember that directories end with a "filesep" and files do not.
% So if you want a fullpath and submit a directory, it needs to end in
% filesep to be correct. Otherwise, you get the "non-exist" error. 
%
% EXAMPLE:
%
% names=ls2cell(fullfile(pwd,'*SAC'));
% names=ls2cell([pwd '/*SAC'])
% names=ls2cell(fullfile(pwd,'*SAC'),1);
% names=ls2cell([pwd '/*SAC'],1)
%
% Last modified by fjsimons-at-alum.mit.edu, 09/10/2014
% Last modified by charig-at-princeton.edu, 09/23/2014
% Last modified by Florian Pfaff for libDirectional, 10/04/2016
if nargin==1
    fullpath=0;
end

% Get the directory-content information
d=dir(ddir);

if ~numel(d)
   error('This directory or file does not exist')
end 
files=str2mat(d.name);

mino=3;
% When using a wildcard that refers directly to a file
if findstr(ddir,'*')
  mino=1;
end

 % Full file paths
if fullpath
    if strcmp(ddir(end),filesep)
        % It's really a directory
        dpath = ddir;
    else
        % It's really some filename wildcard. Get the actual path.
	dpath = [fileparts(ddir) filesep];
    end
    for index=mino:size(files,1)
        cls{index-mino+1}=[dpath deblank(files(index,:))];
    end
else % Just filenames
    for index=mino:size(files,1)
      cls{index-mino+1}=deblank(files(index,:));
    end
end
