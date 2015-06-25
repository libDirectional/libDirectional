function result = runTest(stringOrMetaClass, varargin)
    % Universal front end for different possibilities to run a test.
    %
    % :param stringOrMetaClass: determines a folder path, package name, file
    %   path or TestClass to execute
    % :type stringOrMetaClass: character array or matlab.unittest.meta.class
    % :param method: determines whether just a single method of a Test Class 
    %   should be executed; if empty all methods are executed; optional 
    %   name-value (default is '')
    % :type method: character array
    % :param tapFileName: determines whether results should be stored to
    %   tap-File; if empty results are not stored; optional name-value 
    %   (default is '')
    % :type tapFileName: character array
    % :return: result (matlab.unittest.TestResult) - result of Test Classes or 
    %   single method
    %
    % Calling test class via meta class:
    %
    % .. code-block:: matlab
    %   
    %   test.run(?package.TestClass)
    %
    % Calling test classes of a package via string:
    %
    % .. code-block:: matlab
    %   
    %   test.run('packageName')
    %
    % or
    %
    % .. code-block:: matlab
    %   
    %   test.run('path/to/+packageName')
    %
    % Calling method of test class:
    %
    % .. code-block:: matlab
    %   
    %   test.run(?package.TestClass, 'method', 'methodName')
    %
    % or
    %
    % .. code-block:: matlab
    %   
    %   test.run('path/to/TestClass.m', 'method', 'methodName')
    %
    % .. note: if method is transfered, but stringOrMetaClass does not specify a 
    %   single Test Class, method is going to be ignored.
    %
    % Store results to tap-File:
    %   
    % .. code-block:: matlab
    %
    %   test.run(stringOrMetaClass, 'tapFileName', 'path/to/file.tap')
    %
    %
    % Consult the corresponding test case for further examples.

    import matlab.unittest.TestSuite;
    import matlab.unittest.TestRunner;
    import matlab.unittest.plugins.TAPPlugin;
    import matlab.unittest.plugins.ToFile;
    
    p = inputParser;
    
    addRequired(p, 'stringOrMetaClass', @(x) isa(x,'char') || isa(x, 'matlab.unittest.meta.class') )
    addParameter(p, 'method', '', @ischar)
    addParameter(p, 'tapFileName', '', @ischar)
    parse(p, stringOrMetaClass, varargin{:});
    
    if isa(p.Results.stringOrMetaClass, 'char')
        % stringOrMetaClass is string
        if ~isempty(meta.package.fromName(p.Results.stringOrMetaClass))
            % stringOrMetaClass is string specifying package Name
            if ~isempty(p.Results.method)
                warning('since stringOrMetaClass is a package, parameter method is not used') 
            end
            suite = TestSuite.fromPackage(p.Results.stringOrMetaClass, 'IncludingSubpackages', true);
            
        elseif exist(p.Results.stringOrMetaClass, 'dir')
            % stringOrMetaClass is string specifying folder
            packageInd = strfind(p.Results.stringOrMetaClass, '+');
            if any(packageInd)
                % stringOrMetaClass is string specifying package Folder
                packageName = p.Results.stringOrMetaClass(packageInd(1)+1:end);
                packageName = strrep(packageName, '/+', '.');
                SlashInd = strfind(packageName, '/');
                if ~isempty(SlashInd)
                    if (size(SlashInd, 2) > 1 || SlashInd ~= size(packageName, 2))
                        error('string might contain folder in package')
                    end
                    packageName = strrep(packageName, '/', '');
                end
                if ~isempty(p.Results.method)
                    warning('since stringOrMetaClass is a package, parameter method is not used') 
                end
                suite = TestSuite.fromPackage(packageName, 'IncludingSubpackages', true);
                
            else
                % stringOrMetaClass is string specifying Folder (non-package)
                if ~isempty(p.Results.method)
                    warning('since stringOrMetaClass is a folder, parameter method is not used') 
                end                
                suite = TestSuite.fromFolder(p.Results.stringOrMetaClass);
            end
            
        elseif(exist(fullfile(pwd, p.Results.stringOrMetaClass), 'file') == 2)
            % stringOrMetaClass is string specifying file
            packageInd = strfind(p.Results.stringOrMetaClass, '+');
            if any(packageInd)
                % stringOrMetaClass is string specifying file within package Folder
                className = p.Results.stringOrMetaClass(packageInd(1)+1:end);
                className = strrep(className, '/+', '.');
                if size(strfind(className, '/'), 2) > 1
                    error('string might contain folder in package')
                end
                className = strrep(className, '/', '.');
                className = strrep(className, '.m', '');
                metaClass = meta.class.fromName(className);
                if ~isa(metaClass, 'matlab.unittest.meta.class')
                    error('given file is no TestCase')
                end
                % stringOrMetaClass is string specifying TestClass file within 
                % package Folder 
                if isempty(p.Results.method)
                    % stringOrMetaClass is string specifying TestClass file 
                    % within package Folder, method is not specified.
                    suite = TestSuite.fromClass(metaClass);
                else
                    % stringOrMetaClass is string specifying TestClass file 
                    % within package Folder, method is specified.
                    suite = TestSuite.fromMethod(metaClass, p.Results.method);
                end
                
            else
                % stringOrMetaClass is string specifying file in non-package 
                % Folder
                suite = TestSuite.fromFile(p.Results.stringOrMetaClass);
            end
            
        else
            error('parameter stringOrMetaClass does not fit any package, folder or file')
        end
        
    elseif isa(p.Results.stringOrMetaClass, 'matlab.unittest.meta.class')
        % stringOrMetaClass is meta class
        if isempty(p.Results.method)
            % stringOrMetaClass is meta class, method is not specified.
            suite = TestSuite.fromClass(p.Results.stringOrMetaClass);
        else
            % stringOrMetaClass is meta class, method is specified.
            suite = TestSuite.fromMethod(p.Results.stringOrMetaClass, p.Results.method);
        end
        
    end
    
    % delete old tap-File
    if (~isempty(p.Results.tapFileName) && exist(fullfile(pwd, p.Results.tapFileName), 'file'))
        [status] = system(['rm ' p.Results.tapFileName]);
        if status ~= 0
            error('Failed to delete old tap file.')
        end
    end
    
    % run Test suite by adequate runner
    if isempty(p.Results.tapFileName)
        result = run(suite);
    else
        runner = TestRunner.withTextOutput;
        runner.addPlugin(TAPPlugin.producingOriginalFormat(ToFile(p.Results.tapFileName)));
        result = runner.run(suite);
        disp(fileread(p.Results.tapFileName));
    end
end
