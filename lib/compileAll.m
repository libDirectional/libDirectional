% This script compiles everything.

clear mex %#ok<CLMEX> % % make sure mex files are not open

% switch to path of this file
cd(fileparts(mfilename('fullpath')));

% verify that we are in the correct folder
[~,dirName]=fileparts(pwd);
if ~isequal(dirName,'lib') 
    error('This command has to be run within the lib directory of libDirectional');
end

% util
cd util
compileUtil
cd ..

% externals
cd external


cd Faddeeva_MATLAB
Faddeeva_build
cd ..

cd mhg13
mex mhg.c
mex logmhg.c
mex mhg.c
mex mhgi.c
cd ..

cd TessellateS3
compileTessellate
cd ..

%cd back
cd ..
