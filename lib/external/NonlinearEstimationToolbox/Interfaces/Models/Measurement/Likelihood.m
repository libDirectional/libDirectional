
classdef Likelihood < handle
    % Abstract base class for likelihood functions to represent measurement models.
    %
    % Likelihood Methods:
    %   logLikelihood - Evaluate the logarithmic likelihood function.
    
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
    
    methods (Abstract)
        % Evaluate the logarithmic likelihood function.
        %
        % Parameters:
        %   >> stateSamples (Matrix)
        %      L column-wise arranged state samples.
        %
        %   >> measurement (Arbitrary data)
        %      The measurement data required by the likelihood function.
        %      Usually, this is a column vector.
        %
        % Returns:
        %   << logValues (Row vector)
        %      L column-wise arranged logarithmic likelihood function values.
        logValues = logLikelihood(obj, stateSamples, measurement);
    end
end
