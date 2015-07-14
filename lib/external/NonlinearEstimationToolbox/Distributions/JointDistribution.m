
classdef JointDistribution < Distribution
    % This class represents a joint distribution of other distributions.
    %
    % It is assumed that all distributions are mutually independent.
    %
    % JointDistribution Methods:
    %   JointDistribution    - Class constructor.
    %   getDimension         - Get the dimension of the distribution.
    %   getMeanAndCovariance - Get mean and covariance of the distribution.
    %   drawRndSamples       - Draw random samples from the distribution.
    %   logPdf               - Evaluate the logarithmic probability density function (pdf) of the distribution.
    %   getDistributions     - Get the distributions forming the joint distribution.
    
    % >> This function/class is part of the Nonlinear Estimation Toolbox
    %
    %    For more information, see https://bitbucket.org/nonlinearestimation/toolbox
    %
    %    Copyright (C) 2015  Jannik Steinbring <jannik.steinbring@kit.edu>
    %
    %                        Institute for Anthropomatics and Robotics
    %                        Chair for Intelligent Sensor-Actuator-Systems (ISAS)
    %                        Karlsruhe Institute of Technology (KIT), Germany
    %
    %                        http://isas.uka.de
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
        function obj = JointDistribution(dists)
            % Class constructor.
            %
            % Parameters:
            %   >> dists (Sublass of Distribution or cell array containing subclasses of Distribution)
            %      The distributions forming the joint distribution.
            
            if ~JointDistribution.areDists(dists)
                error('JointDistribution:InvalidDistributions', ...
                      ['dists must be a sublass of Distribution or ' ...
                       'a cell array containing subclasses of Distribution.']);
            end
            
            if ~iscell(dists)
                obj.dists = { dists };
            else
                obj.dists = dists;
            end
            
            obj.numDists = numel(obj.dists);
            
            obj.dimDists = cellfun(@getDimension, obj.dists);
            
            obj.jointDim = sum(obj.dimDists);
            
            obj.jointMean       = [];
            obj.jointCovariance = [];
            obj.jointCovSqrt    = [];
        end
        
        function dim = getDimension(obj)
            dim = obj.jointDim;
        end
        
        function [mean, covariance, covSqrt] = getMeanAndCovariance(obj)
            if isempty(obj.jointMean)
                obj.jointMean       = zeros(obj.jointDim, 1);
                obj.jointCovariance = zeros(obj.jointDim, obj.jointDim);
                obj.jointCovSqrt    = zeros(obj.jointDim, obj.jointDim);
                a = 1;
                b = 0;
                
                for i = 1:obj.numDists
                    b = b + obj.dimDists(i);
                    
                    [obj.jointMean(a:b), ...
                     obj.jointCovariance(a:b, a:b), ...
                     obj.jointCovSqrt(a:b, a:b)] = obj.dists{i}.getMeanAndCovariance();
                    
                    a = b + 1;
                end
            end
            
            mean       = obj.jointMean;
            covariance = obj.jointCovariance;
            
            if nargout >= 3
                if isempty(obj.jointCovSqrt)
                    obj.jointCovSqrt = zeros(obj.jointDim, obj.jointDim);
                    a = 1;
                    b = 0;
                    
                    for i = 1:obj.numDists
                        b = b + obj.dimDists(i);
                        
                        [~, ~, obj.jointCovSqrt(a:b, a:b)] = obj.dists{i}.getMeanAndCovariance();
                        
                        a = b + 1;
                    end
                end
                
                covSqrt = obj.jointCovSqrt;
            end
        end
        
        function samples = drawRndSamples(obj, numSamples)
            if ~Checks.isPosScalar(numSamples)
                error('JointDistribution:InvalidNumberOfSamples', ...
                      'numSamples must be a positive scalar.');
            end
            
            samples = nan(obj.jointDim, numSamples);
            
            a = 1;
            b = 0;
            
            for i = 1:obj.numDists
                b = b + obj.dimDists(i);
                
                samples(a:b, :) = obj.dists{i}.drawRndSamples(numSamples);
                
                a = b + 1;
            end
        end
        
        function logValues = logPdf(obj, values)
            obj.checkValues(values);
            
            numValues = size(values, 2);
            logValues = zeros(1, numValues);
            
            a = 1;
            b = 0;
            
            for i = 1:obj.numDists
                b = b + obj.dimDists(i);
                
                logV = obj.dists{i}.logPdf(values(a:b, :));
                
                logValues = logValues + logV;
                
                a = b + 1;
            end
        end
        
        function [dists, numDists, dimDists] = getDistributions(obj)
            % Get the distributions forming the joint distribution.
            %
            % Returns:
            %   << dists (Cell array containing subclasses of Distribution)
            %      The distributions forming the joint distribution.
            %
            %   << numDists (Scalar)
            %      The number of distributions forming the joint distribution.
            %
        	%   << dimDists (Row vector)
            %      Row-wise arranged dimensions of the distributions forming the joint distribution.
            
            dists    = obj.dists;
            numDists = obj.numDists;
            dimDists = obj.dimDists;
        end
    end
    
    methods (Static, Access = 'private')
        function ret = areDists(dists)
            if Checks.isClass(dists, 'Distribution')
                ret = true;
            elseif iscell(dists)
                L = numel(dists);
                
                if L <= 0
                    ret = false;
                    return;
                end
                
                for i = 1:L
                    if ~Checks.isClass(dists{i}, 'Distribution')
                        ret = false;
                        return;
                    end
                end
                
                ret = true;
            else
                ret = false;
            end
        end
    end
    
    properties (Access = 'private')
        dists;
        numDists;
        dimDists;
        jointDim;
        
        jointMean;
        jointCovariance;
        jointCovSqrt
    end
end

