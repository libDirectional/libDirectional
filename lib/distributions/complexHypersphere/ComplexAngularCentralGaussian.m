classdef ComplexAngularCentralGaussian < AbstractComplexHypersphericalDistribution
    % ComplexAngularCentralGaussian represents a Complex Angular Gaussian
    % distribution

    
    properties
        C       % parameter matrix
    end
    
    methods
        function cAG = ComplexAngularCentralGaussian(C_)
            
            assert(all(all(C_ == C_')), 'C must be hermitian');
            
            cAG.C = C_;
            cAG.dim = size(C_, 1);
        end
        
        function p = pdf(this, za)
            %PDF Evaluates pdf at each column of za
            % Parameters:
            %   za (d x n matrix)
            %       each column represents one of the n points in R^d that the
            %       pdf is evaluated at; can be just a (d x 1) vector as well
            % Returns:
            %   p (1 x n row vector)
            %       values of the pdf at each column of za
            
            p = gamma(this.dim)/(2*pi^this.dim) * abs(dot(za, this.C \ za)).^(-this.dim)./det(this.C);
        end
        
        function Z = sample(this, N)
            % The complex angular centric Gaussian distribution is an
            % alternative to the complex Bingham distribution.
            %
            % Sampling is much faster and the uniform case is equal for both
            % distrinutions.
            %
            % :param N: Number of samples.
            % :type N: Scalar.
            % :return Z: Sampled observation vectors,
            %   complex matrix with size #dimensions times #samples.
            
            R = chol(this.C)';
            a = randn(this.dim, N);
            b = randn(this.dim, N);
            Z = R*complex(a, b);
            
            Z = bsxfun(@rdivide, Z, sqrt(sum(Z .* conj(Z), 1)));
        end
    end
    
    methods (Static)
        function cAG = fit(Z, I)
            if nargin < 2
                I = 100;
            end
            C_ = ComplexAngularCentralGaussian.estimateParameterMatrix(Z, I);
            cAG = ComplexAngularCentralGaussian(C_);
        end
        
        function C_ = estimateParameterMatrix(Z, I)
            if nargin < 2
                I = 100;
            end
            
            D = size(Z, 1);
            N = size(Z, 2);
            
            C_ = eye(D);
            for i = 1:I
                % Z*Z^H for each vector
                Y = bsxfun(@times,permute(Z, [1 3 2]), permute(Z', [3 2 1]));
                
                v = (D - 1)./abs(dot(Z, C_ \ Z));
                v = reshape(v, 1, 1, []);
                
                C_ = sum(bsxfun(@times, v, Y), 3)/N;
            end
        end
    end
end