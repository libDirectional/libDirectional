
classdef CubatureLinearGaussianFilter < SampleBasedLinearGaussianFilter
    % Abstract base class for linear Gaussian filters that are based on the
    % fifth-degree cubature Gaussian sampling technique.
    %
    % CubatureLinearGaussianFilter Methods:
    %   CubatureLinearGaussianFilter - Class constructor.
    %   copy                         - Copy a Filter instance.
    %   copyWithName                 - Copy a Filter instance and give the copy a new name/description.
    %   getName                      - Get the filter name/description.
    %   setColor                     - Set the filter color/plotting properties.
    %   getColor                     - Get the filter color/plotting properties.
    %   setState                     - Set the system state.
    %   getState                     - Get the system state.
    %   setStateMeanAndCov           - Set the system state by means of mean and covariance matrix.
    %   getStateMeanAndCov           - Get mean and covariance matrix of the system state.
    %   getStateDim                  - Get the dimension of the system state.
    %   predict                      - Perform a state prediction.
    %   update                       - Perform a measurement update.
    %   step                         - Perform a combined state prediction and measurement update.
    %   setStateDecompDim            - Set the dimension of the unobservable part of the system state.
    %   getStateDecompDim            - Get the dimension of the unobservable part of the system state.
    %   setPredictionPostProcessing  - Set a post-processing method for the state prediction.
    %   getPredictionPostProcessing  - Get the post-processing method for the state prediction.
    %   setUpdatePostProcessing      - Set a post-processing method for the measurement update.
    %   getUpdatePostProcessing      - Get the post-processing method for the measurement update.
    %   setMeasGatingThreshold       - Set the measurement gating threshold.
    %   getMeasGatingThreshold       - Get the measurement gating threshold.
    
    % Literature:
    %   Bin Jia, Ming Xin, and Yang Cheng,
    %   High-Degree Cubature Kalman Filter,
    %   Automatica, vol. 49, no. 2, pp. 510-518, Feb. 2013.
    
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
        function obj = CubatureLinearGaussianFilter(name)
            % Class constructor.
            %
            % Parameters:
            %   >> name (Char)
            %      An appropriate filter name / description of the implemented
            %      filter. The Filter subclass should set this during its
            %      construction to a meaningful default value (e.g., 'EKF'),
            %      or the user should specify an appropriate name (e.g.,
            %      'PF (10k Particles)').
            %
            % Returns:
            %   << obj (CubatureLinearGaussianFilter)
            %      A new CubatureLinearGaussianFilter instance.
            
            % Call superclass constructor
            obj = obj@SampleBasedLinearGaussianFilter(name);
            
            obj.sampling = GaussianSamplingCKF();
        end
    end
    
    methods (Sealed, Access = 'protected')
        function [stdNormalSamples, ...
                  weights, numSamples] = getStdNormalSamplesPrediction(obj, dim)
            [stdNormalSamples, ...
             weights, numSamples] = obj.sampling.getStdNormalSamples(dim);
        end
        
        function [stdNormalSamples, ...
                  weights, numSamples] = getStdNormalSamplesUpdate(obj, dim)
            [stdNormalSamples, ...
             weights, numSamples] = obj.sampling.getStdNormalSamples(dim);
        end
    end
    
    methods (Access = 'protected')
        % Copy Gaussian sampling technique correctly
        function cpObj = copyElement(obj)
            cpObj = obj.copyElement@SampleBasedLinearGaussianFilter();
            
            cpObj.sampling = obj.sampling.copy();
        end
    end
    
    properties (Access = 'private')
        % Gaussian sampling technique used for the state prediction and measurement update.
        sampling;
    end
end
