function results = runLibDirectionalUnitTests(enableExpensiveTests)
    % Execute all unit tests.

    global enableExpensive
    if nargin>0
        enableExpensive=enableExpensiveTests;
    else
        enableExpensive=false;
    end
    
    warning on % enable all warnings in case they were previously disabled
        
    % switch to path of this file
    cd(fileparts(mfilename('fullpath')));
    
    [~,dirName]=fileparts(pwd);
    if isequal(dirName,'tests')
        suite = matlab.unittest.TestSuite.fromFolder(pwd, 'IncludingSubfolders', true);

        runner = matlab.unittest.TestRunner.withTextOutput();
        tapFile = fullfile(pwd, '../../testResults.tap');
        fprintf('writing tap file: %s\n', tapFile);
        if exist(tapFile, 'file')
            delete(tapFile)
        end
        runner.addPlugin(matlab.unittest.plugins.TAPPlugin.producingOriginalFormat(matlab.unittest.plugins.ToFile(tapFile)));
        results = runner.run(suite);
    else
        error('This command has to be run within the test directory in libDirectional.');
    end
end
