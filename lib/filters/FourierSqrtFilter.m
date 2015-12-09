classdef FourierSqrtFilter < AbstractCircularFilter
    % Fourier filter using the sqrt form.
    %
    % see Florian Pfaff, Gerhard Kurz, and Uwe D. Hanebeck,
    % Multimodal Circular Filtering Using Fourier Series
    % Proceedings of the 18th International Conference on Information Fusion (Fusion 2015), 
    % Washington, D.C., USA, July 2015.
    
    properties
        fd
    end
    
    methods
        function this = FourierSqrtFilter(noOfCoefficients)
            % Constructor
            %
            % Parameters:
            %   noOfCoefficients (integer > 0)
            %       number of Fourier coefficients to use
            this.fd = FourierDistribution.fromDistribution(CircularUniformDistribution(),noOfCoefficients,'sqrt');
        end
        
        function setState(this, fd_)
            % Sets the current system state
            %
            % Parameters:
            %   fd_ (FourierDistribution)
            %       new state
            assert(isa(fd_,'FourierDistribution'));
            this.fd = fd_;
        end
        
        function predictIdentity(this,dSys)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by dSys.
            %
            % Parameters:
            %   dSys (AbstractCircularDistribution)
            %       distribution of additive noise
            assert(isa(dSys,'AbstractCircularDistribution'));
            if ~isa(dSys,'FourierDistribution')
                dSys=FourierDistribution.fromDistribution(dSys,2*length(this.fd.a)-1,'sqrt');
            end
            this.fd=this.fd.convolve(dSys);  
        end
         
        function updateIdentity(this, dMeas, z) 
            % Updates assuming identity measurement model, i.e.,
            % z(k) = x(k) + v(k)    mod 2pi,
            % where v(k) is additive noise given by dMeas.
            %
            % Parameters:
            %   dMeas (AbstractCircularDistribution)
            %       distribution of additive noise
            %   z (scalar)
            %       measurement in [0, 2pi)
            assert(isa(dMeas,'AbstractCircularDistribution'));
            if ~isa(dMeas,'FourierDistribution')
                dMeas=FourierDistribution.fromDistribution(dMeas,2*length(this.fd.a)-1,'sqrt');
            end
            dMeasNew=dMeas.shift(z);
            fdtmp=this.fd.multiply(dMeasNew);
            this.fd=fdtmp.truncate(2*length(this.fd.a)-1);
        end
        
        function fd = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   fd (FourierDistribution)
            %       current estimate
            fd = this.fd;
        end
    end
    
end

