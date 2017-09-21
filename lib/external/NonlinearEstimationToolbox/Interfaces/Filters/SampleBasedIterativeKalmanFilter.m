
classdef SampleBasedIterativeKalmanFilter < IterativeKalmanFilter & SampleBasedLinearGaussianFilter
    % Abstract base class for sample-based iterative Kalman filters.
    %
    % SampleBasedIterativeKalmanFilter Methods:
    %   SampleBasedIterativeKalmanFilter - Class constructor.
    %   copy                             - Copy a Filter instance.
    %   copyWithName                     - Copy a Filter instance and give the copy a new name/description.
    %   getName                          - Get the filter name/description.
    %   setColor                         - Set the filter color/plotting properties.
    %   getColor                         - Get the filter color/plotting properties.
    %   setState                         - Set the system state.
    %   getState                         - Get the system state.
    %   setStateMeanAndCov               - Set the system state by means of mean and covariance matrix.
    %   getStateMeanAndCov               - Get mean and covariance matrix of the system state.
    %   getStateDim                      - Get the dimension of the system state.
    %   predict                          - Perform a state prediction.
    %   update                           - Perform a measurement update.
    %   step                             - Perform a combined state prediction and measurement update.
    %   setStateDecompDim                - Set the dimension of the unobservable part of the system state.
    %   getStateDecompDim                - Get the dimension of the unobservable part of the system state.
    %   setPredictionPostProcessing      - Set a post-processing method for the state prediction.
    %   getPredictionPostProcessing      - Get the post-processing method for the state prediction.
    %   setUpdatePostProcessing          - Set a post-processing method for the measurement update.
    %   getUpdatePostProcessing          - Get the post-processing method for the measurement update.
    %   setMeasGatingThreshold           - Set the measurement gating threshold.
    %   getMeasGatingThreshold           - Get the measurement gating threshold.
    %   setMaxNumIterations              - Set the maximum number of iterations that will be performed by a measurement update.
    %   getMaxNumIterations              - Get the maximum number of iterations that will be performed by a measurement update.
    %   getNumIterations                 - Get number of iterations performed by the last measurement update.
    %   setConvergenceCheck              - Set a convergence check to determine if no further iterations are required.
    %   getConvergenceCheck              - Get the convergence check.
    
    % Literature:
    %   Ángel F. Garcı́a-Fernández, Lennart Svensson, Mark Morelande, and Simo Särkkä,
    %   Posterior Linearisation Filter: Principles and Implementation Using Sigma Points,
    %   IEEE Transactions on Signal Processing, vol. 63, no. 20, pp. 5561-5573, Oct. 2015.
    
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
        function obj = SampleBasedIterativeKalmanFilter(name)
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
            %   << obj (SampleBasedIterativeKalmanFilter)
            %      A new SampleBasedIterativeKalmanFilter instance.
            
            % Call superclass constructors
            obj = obj@IterativeKalmanFilter(name);
            obj = obj@SampleBasedLinearGaussianFilter(name);
        end
    end
    
    methods (Sealed, Access = 'protected')
        function setupMeasModel(obj, measModel, dimMeas)
            [noiseMean, ~, noiseCovSqrt] = measModel.noise.getMeanAndCov();
            
            obj.momentFuncHandle = @(stateMean, stateCovSqrt) ...
                                   obj.momentFuncMeasModel(measModel, dimMeas, ...
                                                           noiseMean, noiseCovSqrt, ...
                                                           stateMean, stateCovSqrt);
        end
        
        function setupAddNoiseMeasModel(obj, measModel, dimMeas)
            [addNoiseMean, addNoiseCov] = measModel.noise.getMeanAndCov();
            dimAddNoise = size(addNoiseMean, 1);
            
            obj.checkAdditiveMeasNoise(dimMeas, dimAddNoise);
            
            obj.momentFuncHandle = @(stateMean, stateCovSqrt) ...
                                   obj.momentFuncAddNoiseMeasModel(measModel, dimMeas, ...
                                                                   addNoiseMean, addNoiseCov, ...
                                                                   stateMean, stateCovSqrt);
        end
        
        function setupMixedNoiseMeasModel(obj, measModel, dimMeas)
            [noiseMean, ~, noiseCovSqrt] = measModel.noise.getMeanAndCov();
            [addNoiseMean, addNoiseCov]  = measModel.additiveNoise.getMeanAndCov();
            dimAddNoise = size(addNoiseMean, 1);
            
            obj.checkAdditiveMeasNoise(dimMeas, dimAddNoise);
            
            obj.momentFuncHandle = @(stateMean, stateCovSqrt) ...
                                   obj.momentFuncMixedNoiseMeasModel(measModel, dimMeas, ...
                                                                     noiseMean, noiseCovSqrt, ...
                                                                     addNoiseMean, addNoiseCov, ...
                                                                     stateMean, stateCovSqrt);
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov] = getMeasMoments(obj, priorMean, ~, priorCovSqrt)
            [measMean, measCov, ...
             stateMeasCrossCov] = obj.momentFuncHandle(priorMean, priorCovSqrt);
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov] = getMeasMomentsIteration(obj, priorStateMean, priorStateCov, priorStateCovSqrt, ...
                                                               updatedStateMean, ~, updatedStateCovSqrt)
            [measMean, measCov, ...
             stateMeasCrossCov] = obj.momentFuncHandle(updatedStateMean, updatedStateCovSqrt);
            
            A = stateMeasCrossCov' / updatedStateCovSqrt';
            H = A / updatedStateCovSqrt;
            P = H * priorStateCovSqrt;
            
            measMean          = measMean + H * (priorStateMean - updatedStateMean);
            measCov           = measCov + P * P' - A * A';
            stateMeasCrossCov = priorStateCov * H';
        end
    end
    
    methods (Access = 'private')
        function [measMean, measCov, ...
                  stateMeasCrossCov] = momentFuncMeasModel(obj, measModel, dimMeas, ...
                                                           noiseMean, noiseCovSqrt, ...
                                                           stateMean, stateCovSqrt)
            [hSamples, weights, ...
             zeroMeanStateSamples] = obj.evaluateMeasModelUncorr(measModel, dimMeas, ...
                                                                 stateMean, stateCovSqrt, ...
                                                                 noiseMean, noiseCovSqrt);
            
            % Compute moments:
            %  - Measurement mean:
            %    E[y] = E[h(x,v)]
            %  - Measurement covariance matrix:
            %    E[(y - E[y])*(y - E[y])'] = E[(h(x,v) - E[h(x,v)])*(h(x,v) - E[h(x,v)])']
            %  - State--measurement cross-covariance matrix:
            %    E[(x - E[x])*(y - E[y])'] = E[(x - E[x])*(h(x,v) - E[h(x,v)])']
            [measMean, measCov, ...
             stateMeasCrossCov] = obj.getMeasModelMomements(hSamples, weights, ...
                                                            zeroMeanStateSamples);
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov] = momentFuncAddNoiseMeasModel(obj, measModel, dimMeas, ...
                                                                   addNoiseMean, addNoiseCov, ...
                                                                   stateMean, stateCovSqrt)
            [hSamples, weights, ...
             zeroMeanStateSamples] = obj.evaluateAddNoiseMeasModel(measModel, dimMeas, ...
                                                                   stateMean, stateCovSqrt);
          	
            % Compute moments:
            %  - E[h(x)]
            %  - E[(h(x) - E[h(x)])*(h(x) - E[h(x)])']
            %  - E[(x - E[x])*(h(x) - E[h(x)])']
            [hMean, hCov, ...
             stateHCrossCov] = obj.getMeasModelMomements(hSamples, weights, ...
                                                         zeroMeanStateSamples);
            
            % Measurement mean:
            % E[y] = E[h(x)] + E[v]
            measMean = hMean + addNoiseMean;
            
            % Measurement covariance matrix:
            % E[(y - E[y])*(y - E[y])'] = E[(h(x) + v - E[h(x)] - E[v])*(h(x) + v - E[h(x)] - E[v])']
            %                           = E[(h(x) - E[h(x)])*(h(x) - E[h(x)])']
            %                           + E[(v - E[v])*(v - E[v])']
            %                           + E[(h(x) - E[h(x)])*(v - E[v])']
            %                           + E[(v - E[v])*(h(x) - E[h(x)])']
            %                           = E[(h(x) - E[h(x)])*(h(x) - E[h(x)])']
            %                           + E[(v - E[v])*(v - E[v])']
            measCov = hCov + addNoiseCov;
            
            % State--measurement cross-covariance matrix
            % E[(x - E[x])*(y - E[y])'] = E[(x - E[x])*(h(x) + v - E[h(x)] - E[v])']
            %                           = E[(x - E[x])*(h(x) - E[h(x)])']
            %                           + E[(x - E[x])*(v - E[v])']
            %                           = E[(x - E[x])*(h(x) - E[h(x)])']
            stateMeasCrossCov = stateHCrossCov;
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov] = momentFuncMixedNoiseMeasModel(obj, measModel, dimMeas, ...
                                                                     noiseMean, noiseCovSqrt, ...
                                                                     addNoiseMean, addNoiseCov, ...
                                                                     stateMean, stateCovSqrt)
            [hSamples, weights, ...
             zeroMeanStateSamples] = obj.evaluateMeasModelUncorr(measModel, dimMeas, ...
                                                                 stateMean, stateCovSqrt, ...
                                                                 noiseMean, noiseCovSqrt);
            
            % Compute moments:
            %  - E[h(x,v)]
            %  - E[(h(x,v) - E[h(x,v)])*(h(x,v) - E[h(x,v)])']
            %  - E[(x - E[x])*(h(x,v) - E[h(x,v)])']
            [hMean, hCov, ...
             stateHCrossCov] = obj.getMeasModelMomements(hSamples, weights, ...
                                                         zeroMeanStateSamples);
            
            % Measurement mean:
            % E[y] = E[h(x,v)] + E[r]
            measMean = hMean + addNoiseMean;
            
            % Measurement covariance matrix:
            % E[(y - E[y])*(y - E[y])'] = E[(h(x,v) + r - E[h(x,v)] - E[r])*(h(x,v) + r - E[h(x,v)] - E[r])']
            %                           = E[(h(x,v) - E[h(x,v)])*(h(x,v) - E[h(x,v)])']
            %                           + E[(r - E[r])*(r - E[r])']
            %                           + E[(h(x,v) - E[h(x,v)])*(r - E[r])']
            %                           + E[(r - E[r])*(h(x,v) - E[h(x,v)])']
            %                           = E[(h(x,v) - E[h(x,v)])*(h(x,v) - E[h(x,v)])']
            %                           + E[(r - E[r])*(r - E[r])']
            measCov = hCov + addNoiseCov;

            % State--measurement cross-covariance matrix
            % E[(x - E[x])*(y - E[y])'] = E[(x - E[x])*(h(x,v) + r - E[h(x,v)] - E[r])']
            %                           = E[(x - E[x])*(h(x,v) - E[h(x,v)])']
            %                           + E[(x - E[x])*(r - E[r])']
            %                           = E[(x - E[x])*(h(x,v) - E[h(x,v)])']
            stateMeasCrossCov = stateHCrossCov;
        end
    end
    
    properties (Access = 'private')
        % Function handle to the currently used moment computation method.
        momentFuncHandle;
    end
end
