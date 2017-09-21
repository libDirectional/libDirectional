
classdef ParticleFilter < Filter
    % Abstract base class for particle filters.
    %
    % ParticleFilter Methods:
    %   ParticleFilter     - Class constructor.
    %   copy               - Copy a Filter instance.
    %   copyWithName       - Copy a Filter instance and give the copy a new name/description.
    %   getName            - Get the filter name/description.
    %   setColor           - Set the filter color/plotting properties.
    %   getColor           - Get the filter color/plotting properties.
    %   setState           - Set the system state.
    %   getState           - Get the system state.
    %   setStateMeanAndCov - Set the system state by means of mean and covariance matrix.
    %   getStateMeanAndCov - Get mean and covariance matrix of the system state.
    %   getStateDim        - Get the dimension of the system state.
    %   predict            - Perform a state prediction.
    %   update             - Perform a measurement update.
    %   step               - Perform a combined state prediction and measurement update.
    
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
        function obj = ParticleFilter(name)
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
            %   << obj (ParticleFilter)
            %      A new ParticleFilter instance.
            
            % Call superclass constructor
            obj = obj@Filter(name);
        end
    end
    
    methods (Access = 'protected')
        function predictedParticles = predictParticlesArbitraryNoise(obj, sysModel, particles, numParticles)
            % Sample system noise
            noise = sysModel.noise.drawRndSamples(numParticles);
            
            % Propagate particles and noise through system equation
            predictedParticles = sysModel.systemEquation(particles, noise);
            
            % Check predicted particles
            obj.checkPredictedStateSamples(predictedParticles, numParticles);
        end
        
        function predictedParticles = predictParticlesAdditiveNoise(obj, sysModel, particles, numParticles)
            % Sample additive system noise
            noise = sysModel.noise.drawRndSamples(numParticles);
            
            dimNoise = size(noise, 1);
            
            obj.checkAdditiveSysNoise(dimNoise);
            
            % Propagate particles through system equation
            predictedParticles = sysModel.systemEquation(particles);
            
            % Check predicted particles
            obj.checkPredictedStateSamples(predictedParticles, numParticles);
            
            % Add system noise
            predictedParticles = predictedParticles + noise;
        end
        
        function predictedParticles = predictParticlesMixedNoise(obj, sysModel, particles, numParticles)
            % Sample system noise
            noise = sysModel.noise.drawRndSamples(numParticles);
            
            % Sample additive system noise
            addNoise = sysModel.additiveNoise.drawRndSamples(numParticles);
            
            dimAddNoise = size(addNoise, 1);
            
            obj.checkAdditiveSysNoise(dimAddNoise);
            
            % Propagate particles and noise through system equation
            predictedParticles = sysModel.systemEquation(particles, noise);
            
            % Check predicted particles
            obj.checkPredictedStateSamples(predictedParticles, numParticles);
            
            % Add system noise
            predictedParticles = predictedParticles + addNoise;
        end
        
        function values = evaluateLikelihood(obj, measModel, measurement, particles, numParticles)
            % Evaluate logarithmic likelihood
            logValues = measModel.logLikelihood(particles, measurement);
            
            obj.checkLogLikelihoodEvaluations(logValues, numParticles);
            
            % For numerical stability
            maxLogValue = max(logValues);
            
            logValues = logValues - maxLogValue;
            
            % Compute likelihohod values
            values = exp(logValues);
        end
    end
end
