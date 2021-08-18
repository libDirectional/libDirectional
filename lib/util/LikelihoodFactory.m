classdef LikelihoodFactory
    % Creates likelihood functions
    
    methods (Static)
        function f = additiveNoiseLikelihood(h, noiseDistribution)
            % Creates a likelihood function for additive noise.
            % Returns f(z,x) = P(z|x) for z = h(x) + v with noise v
            % 
            % Parameters:
            %   h (function handle)
            %      function mapping the d-dimensional state to the measurement (needs to
            %      support vectorized inputs, i.e., d x n matrices)
            %   noiseDistribution (any suitable class with a pdf function)
            %       probability distribution of noise
            % Returns:
            %   f (function handle)
            %       likelihood function f(z,x)
            arguments
                h (1,1) function_handle
                noiseDistribution (1,1) AbstractDistribution
            end
            
            %use repmat if z contains one measurement, but x contains many
            %state vecors
            f = @(z,x) noiseDistribution.pdf(z-h(x)); % Use implicit expansion
        end
    end
    
end

