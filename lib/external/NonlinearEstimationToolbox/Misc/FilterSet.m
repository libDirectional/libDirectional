
classdef FilterSet < handle
    % A set of filters.
    %
    % FilterSet Methods:
    %   FilterSet           - Class constructor.
    %   add                 - Add a filter to the set.
    %   remove              - Remove a filter from the set
    %   get                 - Get a filter from the set.
    %   getIndex            - Get the index of a filter from the set.
    %   getNumFilters       - Get the number of filters in the set.
    %   getNames            - Get the names of all filters.
    %   setStates           - Set the system states of all filters.
    %   getStates           - Get the system states of all filters.
    %   setStatesMeanAndCov - Set the system states of all filters by means of mean and covariance matrix.
    %   getStatesMeanAndCov - Get state means and state covariance matrices of all filters.
    %   getStateDim         - Get the dimension of the system state.
    %   predict             - Perform a state prediction for all filters.
    %   predictSingle       - Perform a state prediction for a particular filter.
    %   update              - Perform a measurement update for all filters.
    %   updateSingle        - Perform a measurement update for a particular filter.
    %   step                - Perform a combined state prediction and measurement update for all filters.
    %   stepSingle          - Perform a combined state prediction and measurement update for a particular filter.
    
    % >> This function/class is part of the Nonlinear Estimation Toolbox
    %
    %    For more information, see https://bitbucket.org/nonlinearestimation/toolbox
    %
    %    Copyright (C) 2015-2017  Jannik Steinbring <nonlinearestimation@gmail.com>
    %
    %    This program is free software: you can redistribute it and/or modify
    %    it under the terms of the GNU General Public License as published by
    %    the Free Software Foundation, either version 3 of the License, or
    %    (at your option) any later version.
    %
    %    This program is distributed in the hope that it will be useful,
    %    but WITHOUT ANY WARRANTY; without even the implied warranty of
    %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %    GNU General Public License for more details.
    %
    %    You should have received a copy of the GNU General Public License
    %    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    methods
        function obj = FilterSet()
            % Default constructor.
            %
            % Returns:
            %   << obj (FilterSet)
            %      A new FilterSet instance.
            
            obj.filterMap = containers.Map('KeyType', 'char', ...
                                           'ValueType', 'any');
            
            obj.filters    = obj.filterMap.values();
            obj.numFilters = obj.filterMap.length();
            obj.dimState   = 0;
        end
        
        function add(obj, filter)
            % Add a filter to the set.
            %
            % Parameters:
            %   >> filter (Subclass of Filter)
            %      The filter to be added.
            
            if ~Checks.isClass(filter, 'Filter')
                error('FilterSet:InvalidFilter', ...
                      'filter must be a subclass of Filter.');
            end
            
            filterName = filter.getName();
            
            index = obj.getMapIndex(filterName);
            
            if ~isempty(index)
                error('FilterSet:InvalidFilterName', ...
                      'Filter name "%s" already exists.', filterName);
            end
            
            obj.filterMap(filterName) = filter;
            
            obj.filters    = obj.filterMap.values();
            obj.numFilters = obj.filterMap.length();
        end
         
        function remove(obj, filterName)
            % Remove a filter from the set.
            %
            % Parameters:
            %   >> name (Char)
            %      Name of the filter to be removed.
            
            if ~ischar(filterName)
                error('FilterSet:InvalidFilterName', ...
                      'filterName must be a char.');
            end
            
            index = obj.getMapIndex(filterName);
            
            if isempty(index)
                error('FilterSet:NoFilter', ...
                      'There exists no filter with the name "%s".', filterName);
            end
            
            obj.filterMap.remove(filterName);
            
            obj.filters    = obj.filterMap.values();
            obj.numFilters = obj.filterMap.length();
        end
        
        function filter = get(obj, id)
            % Get a filter from the set.
            %
            % Parameters:
            %   >> id (Char or positive scalar)
            %      The filter id can be either the filter name or its index in the set.
            %
            % Returns:
            %   << filter (Subclass of Filter)
            %      The filter corresponding to the given id.
            
            if ischar(id)
                try
                    filter = obj.filterMap(id);
                catch ex
                    if strcmp(ex.identifier, 'MATLAB:Containers:Map:NoKey')
                        error('FilterSet:NoFilter', ...
                              'There exists no filter with the name "%s".', id);
                    else
                        rethrow(ex);
                    end
                end
            elseif Checks.isPosScalar(id)
                try
                    filter = obj.filters{id};
                catch ex
                    if strcmp(ex.identifier, 'MATLAB:badsubscript')
                        error('FilterSet:NoFilter', ...
                              'Filter index out of range.');
                    else
                        rethrow(ex);
                    end
                end
            else
                error('FilterSet:InvalidIdentifier', ...
                      'id must be a filter name or a positive index.');
            end
        end
        
        function index = getIndex(obj, filterName)
            % Get the index of a filter from the set.
            %
            % Note: the index of a filter may have changed when a filter was added or removed.
            %
            % Parameters:
            %   >> name (Char)
            %      Name of the filter to be removed.
            %
            % Returns:
            %   << index (Scalar)
            %      Current index of the filter in the set.
            
            if ischar(filterName)
                index = obj.getMapIndex(filterName);
                
                if isempty(index)
                    error('FilterSet:NoFilter', ...
                          'There exists no filter with the name "%s".', filterName);
                end
            else
                error('FilterSet:InvalidFilterName', ...
                      'filterName must be a char.');
            end
        end
        
        function numFilters = getNumFilters(obj)
            % Get the number of filters in the set.
            %
            % Returns:
            %   << numFilters (Scalar)
            %      The current number of filters in the set.
            
            numFilters = obj.numFilters;
        end
        
        function names = getNames(obj)
            % Get the names of all filters.
            %
            % Returns:
            %   << names (Cell array of chars)
            %      Each cell array element contains the name of a filter.
            
            names = obj.filterMap.keys();
        end
        
        function setStates(obj, state)
            % Set the system states of all filters.
            %
            % Parameters:
            %   >> state (Subclass of Distribution)
            %      The new system state.
            
            obj.forAllFilters(@setState, state);
            
            obj.dimState = state.getDim();
        end
        
        function states = getStates(obj)
            % Get the system states of all filters.
            %
            % Returns:
            %   << states (Cell array of Distributions)
            %      Each cell array element contains the system state of a filter.
            
            states = cell(1, obj.numFilters);
            
            for i = 1:obj.numFilters
                states{i} = obj.filters{i}.getState();
            end
        end
        
        function setStatesMeanAndCov(obj, stateMean, stateCov, stateCovSqrt)
            % Set the system states of all filters by means of mean and covariance matrix.
            %
            % Note: this method does not perform input validation like setStates()!
            % It is intended for fastly setting the system state without creating
            % a temporary Gaussian distribution, e.g., in order to assign the
            % Gaussian state estimate of a filter to another one.
            %
            % Parameters:
            %   >> stateMean (Column vector)
            %      The new mean vector of the system state.
            %
            %   >> stateCov (Positive definite matrix)
            %      The new covariance matrix of the system state.
            %
            %   >> stateCovSqrt (Square matrix)
            %      Lower Cholesky decomposition of the new system state covariance matrix.
            %      If no square root is passed, it will be computed.
            
            if nargin < 4
                stateCovSqrt = chol(stateCov, 'Lower');
            end
            
            obj.forAllFilters(@setStateMeanAndCov, stateMean, stateCov, stateCovSqrt);
            
            obj.dimState = size(stateMean, 1);
        end
        
        function [stateMeans, stateCovs, stateCovSqrts] = getStatesMeanAndCov(obj)
            % Get state means and state covariance matrices of all filters.
            %
            % Returns:
            %   << stateMeans (Matrix)
            %      Column-wise arranged state means.
            %
            %   << stateCovs (3D matrix of positive definite matrices)
            %      State covariance matrices arranged along the 3rd dimension.
            %
            %   << stateCovSqrts (3D matrix)
            %      Lower Cholesky decompositions of the system state covariance matrices
            %      arranged along the 3rd dimension.
            
            if nargout == 1
                stateMeans = nan(obj.dimState, obj.numFilters);
                
                for i = 1:obj.numFilters
                    stateMeans(:, i) = obj.filters{i}.getStateMeanAndCov();
                end
            elseif nargout == 2
                stateMeans = nan(obj.dimState, obj.numFilters);
                stateCovs  = nan(obj.dimState, obj.dimState, obj.numFilters);
                
                for i = 1:obj.numFilters
                    [stateMeans(:, i), stateCovs(:, :, i)] = obj.filters{i}.getStateMeanAndCov();
                end
            else
                stateMeans    = nan(obj.dimState, obj.numFilters);
                stateCovs     = nan(obj.dimState, obj.dimState, obj.numFilters);
                stateCovSqrts = nan(obj.dimState, obj.dimState, obj.numFilters);
                
                for i = 1:obj.numFilters
                    [stateMeans(:, i), stateCovs(:, :, i), stateCovSqrts(:, :, i)] = obj.filters{i}.getStateMeanAndCov();
                end
            end
        end
        
        function dim = getStateDim(obj)
            % Get the dimension of the system state.
            %
            % Returns:
            %   << dim (Scalar)
            %      The dimension of the system state.
            
            dim = obj.dimState;
        end
        
        function runtimes = predict(obj, sysModel)
            % Perform a state prediction for all filters.
            %
            % Parameters:
            %   >> sysModel (Arbitrary class; filter dependent)
            %      System model that describes the temporal behavior of the system state.
            %
            % Returns:
            %   << runtimes (Row vector)
            %      Contains the state prediction runtimes for all filters.
            
            if nargout == 1
                runtimes = obj.forAllFiltersRuntime(@predict, sysModel);
            else
                obj.forAllFilters(@predict, sysModel);
            end
        end
        
        function runtime = predictSingle(obj, id, sysModel)
            % Perform a state prediction for a particular filter.
            %
            % Parameters:
            %   >> id (Char or positive scalar)
            %      The filter id can be either the filter name or its index in the set.
            %
            %   >> sysModel (Arbitrary class; filter dependent)
            %      System model that describes the temporal behavior of the system state.
            %
            % Returns:
            %   << runtime (Scalar)
            %      Time needed to perform the state prediction.
            
            if nargout == 1
                runtime = obj.executeFilterRuntime(id, @predict, sysModel);
            else
                obj.executeFilter(id, @predict, sysModel);
            end
        end
        
        function runtimes = update(obj, measModel, measurement)
            % Perform a measurement update for all filters.
            %
            % Parameters:
            %   >> measModel (Arbitrary class; filter dependent)
            %      Measurement model that describes the relationship between system state and measurement.
            %
            %   >> measurement (Arbitrary data)
            %      Measurement data that has to be processed by the measurement update.
            %      What type of data is supported depends on the passed measurement model
            %      and implemented filter. Usually, this is a column vector.
            %
            % Returns:
            %   << runtimes (Row vector)
            %      Contains the measurement update runtimes for all filters.
            
            if nargout == 1
                runtimes = obj.forAllFiltersRuntime(@update, measModel, measurement);
            else
                obj.forAllFilters(@update, measModel, measurement);
            end
        end
        
        function runtime = updateSingle(obj, id, measModel, measurement)
            % Perform a measurement update for a particular filter.
            %
            % Parameters:
            %   >> id (Char or positive scalar)
            %      The filter id can be either the filter name or its index in the set.
            %
            %   >> measModel (Arbitrary class; filter dependent)
            %      Measurement model that describes the relationship between system state and measurement.
            %
            %   >> measurement (Arbitrary data)
            %      Measurement data that has to be processed by the measurement update.
            %      What type of data is supported depends on the passed measurement model
            %      and implemented filter. Usually, this is a column vector.
            %
            % Returns:
            %   << runtime (Scalar)
            %      Time needed to perform the measurement update.
            
            if nargout == 1
                runtime = obj.executeFilterRuntime(id, @update, measModel, measurement);
            else
                obj.executeFilter(id, @update, measModel, measurement);
            end
        end
        
        function runtimes = step(obj, sysModel, measModel, measurement)
            % Perform a combined state prediction and measurement update for all filters.
            %
            % Parameters:
            %   >> sysModel (Arbitrary class; filter dependent)
            %      System model that describes the temporal behavior of the system state.
            %
            %   >> measModel (Arbitrary class; filter dependent)
            %      Measurement model that describes the relationship between system state and measurement.
            %
            %   >> measurement (Arbitrary data)
            %      Measurement data that has to be processed by the measurement update.
            %      What type of data is supported depends on the passed measurement model
            %      and implemented filter. Usually, this is a column vector.
            %
            % Returns:
            %   << runtimes (Row vector)
            %      Contains the runtimes for all filters.
            
            if nargout == 1
                runtimes = obj.forAllFiltersRuntime(@step, sysModel, measModel, measurement);
            else
                obj.forAllFilters(@step, sysModel, measModel, measurement);
            end
        end
        
        function runtime = stepSingle(obj, id, sysModel, measModel, measurement)
            % Perform a combined state prediction and measurement update for a particular filter.
            %
            % Parameters:
            %   >> sysModel (Arbitrary class; filter dependent)
            %      System model that describes the temporal behavior of the system state.
            %
            %   >> measModel (Arbitrary class; filter dependent)
            %      Measurement model that describes the relationship between system state and measurement.
            %
            %   >> measurement (Arbitrary data)
            %      Measurement data that has to be processed by the measurement update.
            %      What type of data is supported depends on the passed measurement model
            %      and implemented filter. Usually, this is a column vector.
            %
            % Returns:
            %   << runtime (Scalar)
            %      Time needed to perform the combined state prediction and measurement update.
            
            if nargout == 1
                runtime = obj.executeFilterRuntime(id, @step, sysModel, measModel, measurement);
            else
                obj.executeFilter(id, @step, sysModel, measModel, measurement);
            end
        end
    end
    
    methods (Access = 'protected')
        function index = getMapIndex(obj, name)
            index = find(strcmp(obj.filterMap.keys(), name), 1);
        end
        
        function forAllFilters(obj, func, varargin)
            cellfun(@(f) func(f, varargin{:}), obj.filters);
        end
        
        function runtimes = forAllFiltersRuntime(obj, func, varargin)
            runtimes = cellfun(@(f) func(f, varargin{:}), obj.filters);
        end
        
        function executeFilter(obj, id, func, varargin)
            filter = obj.get(id);
            
            func(filter, varargin{:});
        end
        
        function runtimes = executeFilterRuntime(obj, id, func, varargin)
            filter = obj.get(id);
            
            runtimes = func(filter, varargin{:});
        end
    end
    
    properties (Access = 'protected')
        % Key (filter name) / value (filter itself) storage.
        filterMap;
        
        % Cell array of filters.
        filters;
        
        % Number of filters hold by the set.
        numFilters;
        
        % State dimension of all filters.
        dimState;
    end
end
