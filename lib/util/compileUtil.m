% This script compiles everything in util.

clear mex %#ok<CLMEX> % make sure mex files are not open

% switch to path of this file
cd(fileparts(mfilename('fullpath')));

% verify that we are in the correct folder
[~,dirName]=fileparts(pwd);
if isequal(dirName,'util') 
    utilDir='.';
else
    error('This command has to be run within the util directory in libDirectional');
end

if ispc
    % Windows
    
    % enable AVX if the CPU and OS support it
    % optimflags = 'OPTIMFLAGS=$OPTIMFLAGS /openmp /Ox /arch:AVX';
    
    % otherwise disable AVX
    optimflags = 'OPTIMFLAGS=$OPTIMFLAGS /openmp /Ox';
    options = {optimflags};
elseif ismac
    % Mac
    options = {};
else
    % Linux
    cxxFlags = '-std=c++0x -Wall -Wfatal-errors -march=native -fopenmp';
    ldFlags = '-fopenmp';
    options  = { ['CXXFLAGS=$CXXFLAGS ' cxxFlags], ...
                     ['LDFLAGS=$LDFLAGS ' ldFlags] };
end   

%options = {'-v', options{:} };

mex(options{:}, fullfile(utilDir,'besselratio.cpp'))
mex(options{:}, fullfile(utilDir,'numericalSaddlepointWithDerivatives.cpp'), fullfile(utilDir,'binghamNormalizationConstant.cpp'))
mex(options{:}, fullfile(utilDir,'numericalBinghamMLE.cpp'),fullfile(utilDir,'binghamNormalizationConstant.cpp'))
mex(options{:},fullfile(utilDir,'wnpdf.cpp'))
mex(options{:},fullfile(utilDir,'toroidalwnpdf.cpp'))
mex(options{:},fullfile(utilDir,'mvnpdffast.cpp'))

mex(options{:},fullfile(utilDir,'hypergeometricRatioInverse.cpp'));

mex(options{:},fullfile(utilDir,'glover.cpp'));


for d = 2:6
    fprintf('Generate helper functions for d = %d\n', d);
    generateSymbolic(d);
end
