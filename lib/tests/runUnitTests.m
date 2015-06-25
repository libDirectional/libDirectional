function results = runUnitTests()
    % Execute all unit tests.
    
    % switch to path of this file
    cd(fileparts(mfilename('fullpath')));
    
    [~,dirName]=fileparts(pwd);
    if isequal(dirName,'tests')
        results = runtests(pwd(), 'Recursively', true);
    else
        error('This command has to be run within the test directory in libDirectional.');
    end
end
