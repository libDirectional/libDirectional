classdef BinghamFilter < AbstractAxialFilter
    % Implements a recursive filter based on the Bingham distribution.
    %
    % For antipodally symmetric complex numbers (2D) and quaternions (4D)
    % only.
    %
    % Gerhard Kurz, Igor Gilitschenski, Simon Julier, Uwe D. Hanebeck,
    % Recursive Bingham Filter for Directional Estimation Involving 180 Degree Symmetry
    % Journal of Advances in Information Fusion, 9(2):90 - 105, December 2014.
    %
    % Igor Gilitschenski, Gerhard Kurz, Simon J. Julier, Uwe D. Hanebeck,
    % Unscented Orientation Estimation Based on the Bingham Distribution (accepted)
    % IEEE Transactions on Automatic Control, January 2016.
    
    properties
        B BinghamDistribution % state as Bingham distribution
    end
    
    methods
        function this = BinghamFilter()
            % Constructor
            B_ = BinghamDistribution([-1;-1;-1;0], eye(4));
            this.setState(B_);
        end
        
        function setState(this, B_)
            % Sets the current system state
            %
            % Parameters:
            %   B_ (BinghamDistribution)
            %       new state
            assert(isa(B_, 'BinghamDistribution'));
            assert (B_.dim==2 || B_.dim==4, 'Only 2D and 4D distributions are supported.')
            this.B = B_;
            this.setCompositionOperator();
        end
        
        function predictIdentity(this, Bw)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) (+) w(k)    
            % where w(k) is noise given by Bw.
            % The composition operator (+) refers to a complex or quaternion
            % multiplication.
            %
            % Parameters:
            %   Bw (BinghamDistribution)
            %       distribution of noise
            assert(isa(Bw, 'BinghamDistribution'));
            this.B = this.B.compose(Bw);
        end
        
        function predictNonlinear(this, a, Bw)
            % Predicts assuming nonlinear system model with "additive" noise, i.e.,
            % x(k+1) = a(x(k)) (+) w(k)    mod 2pi,
            % where w(k) is Bingham distributed noise represented by Bw.
            % The composition operator (+) refers to a complex or quaternion
            % multiplication.
            %
            % Parameters:
            %   a  (reference to function)
            %       nonlinear system function
            %   Bw (BinghamDistribution)
            %       distribution of noise
          
            [deterministicSamples, sampleWeights] = this.B.sampleDeterministic(0.5);

            % Propagate each sample through the system function.
            for i=1:numel(sampleWeights)
                deterministicSamples(:,i) = a(deterministicSamples(:,i));
            end

            % Vectorized version.
            %deterministicSamples = a(deterministicSamples);

            % The following computes Cov(g(x) (+) w).
            S = deterministicSamples*diag(sampleWeights)*deterministicSamples';
            S = (S+S')/2;
            predictedBingham = BinghamDistribution.fitToMoment(S);
            this.B = predictedBingham.compose(Bw);
        end
        
        function updateIdentity(this, Bv, z)
            % Updates assuming identity measurement model, i.e.,
            % z(k) = x(k) (+) v(k)    mod 2pi,
            % where v(k) is additive noise given by vmMeas.
            % The composition operator (+) refers to a complex or quaternion
            % multiplication.
            %
            % Parameters:
            %   Bv (BinghamDistribution)
            %       distribution of additive noise
            %   z (dim x 1 vector)
            %       measurement on the unit hypersphere
            assert(isa(Bv, 'BinghamDistribution'));
            assert(Bv.dim == this.B.dim);
            assert(isequal(size(z), [this.B.dim, 1]));
            assert(abs(norm(z) - 1) < 1E-10);
            
            conjugate = @(q) [q(1); -q(2:end)];
            for i=1:this.B.dim
                %Bv.M(:,i) = this.compositionOperator(conjugate(Bv.M(:,i)), z');
                Bv.M(:,i) = this.compositionOperator(z,conjugate(Bv.M(:,i)));
            end
            this.B = this.B.multiply(Bv);
        end
        
        function B = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   B (BinghamDistribution)
            %       current estimate
            B = this.B;
        end
            
    end
    
end

