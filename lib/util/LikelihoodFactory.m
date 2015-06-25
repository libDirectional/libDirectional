classdef LikelihoodFactory
    % Creates likelihood functions
    
    methods (Static)
        function f = additiveNoiseLikelihood(h, noiseDistribution)
            % Creates a likelihood function for additive noise.
            % Returns f(z,x) = P(z|x) for z = h(x) + v with noise v
            % 
            % Parameters:
            %   h (function handle)
            %      function mapping the state to the measurement 
            %   noiseDistribution (any suitable class with a pdf function)
            %       probability distribution of noise
            % Returns:
            %   f (function handle)
            %       likelihood function f(z,x)
            assert(isa(h, 'function_handle'));
            
            f = @(z,x) noiseDistribution.pdf(z-h(x));
        end
    end
    
end

