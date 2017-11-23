classdef SphericalHarmonicFilter < AbstractHypersphericalFilter
    % Filter using spherical harmonics
    
    % see Florian Pfaff, Gerhard Kurz, and Uwe D. Hanebeck,
    % Filtering on the Unit Sphere Using Spherical Harmonics
    % Proceedings of the 2017 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI 2017), 
    % Daegu, Korea, November 2017.
    properties
        state
    end
        
    methods 
        function this = SphericalHarmonicFilter(degree,transformation)
            if nargin==1,transformation='identity';end
            assert(strcmp(transformation,'identity'),'Currently only identity transformation is supported.');
            coeffMat=zeros(degree+1,2*degree+1);
            coeffMat(1)=1/sqrt(4*pi);
            this.state=SphericalHarmonicDistributionComplex(coeffMat,transformation);
            this.dim=3;
        end

        function setState(this, state)
            assert(isa(state, 'AbstractSphericalHarmonicDistribution'));
            if isa(state,'SphericalHarmonicDistributionReal')
                state=state.toComplexSphericalHarmonicDistribution;
            end
            this.state=state.truncate(size(this.state.coeffMat,1)-1);
        end
        
        function est = getEstimate(this)
           est = this.state;
        end
        
        function mean = getEstimateMean(this)
            mean=this.state.mean;
        end
        
        function predictIdentity(this, sysNoise)
            % sysNoise has to be zonal.
            assert(isa(sysNoise, 'AbstractSphericalHarmonicDistribution'));
            if isa(sysNoise,'RealSphericalHarmonicDistribution')
                sysNoise=sysNoise.toComplexSphericalHarmonicDistribution;
            end
            this.setState(this.state.convolve(sysNoise));
        end
        
        function updateIdentity(this, measNoise, z)
            % Splitting the noise into a noise term and a measurement is not a
            % universal solution. This function was added for interface
            % compatibility of the VMF Filter. If z is used, measNoise needs to be
            % zonal with rotational symmetry around the z-axis and have a
            % mean of [0;0;1].
            % It is possible to use 0 or [0;0;1] as measurement if no rotation should be
            % performed. Please note that it is usually more efficient to encode it in
            % measNoise than to perform a rotation. 
            
            assert(isa(measNoise, 'AbstractSphericalHarmonicDistribution'));
            if isa(measNoise,'RealSphericalHarmonicDistribution')
                measNoise=measNoise.toComplexSphericalHarmonicDistribution;
            end
            
            if (norm(z)>1E-6)&&norm(z-[0;0;1])>1E-6 
                [phi,theta]=cart2sph(z(1),z(2),z(3));
                theta=-theta+pi/2;
                measNoise=measNoise.rotate(0,theta,phi);
            end
            this.setState(this.state.multiply(measNoise));
        end
        
        function updateNonlinearUsingLikelihood(this, likelihood, z) %measurement z, likelihood(z,x)=P(Z|X)
            % No zonality is required. Can use array valued likelihoods and
            % matrix using cell array of likelihoods and measurements.
            switch this.state.transformation
                case 'identity'
                    % If identity transformation is used and a likelihood is given,
                    % we can simply directly evaluate the likelihood at this
                    % positions.
                    degree=size(this.state.coeffMat,1)-1;
                    lat=linspace(0,2*pi,2*degree+2);
                    lon=linspace(pi/2,-pi/2,degree+2);
                    [latMesh,lonMesh]=meshgrid(lat,lon);
                    [x1,x2,x3]=sph2cart(latMesh,lonMesh,1);
                    fvalCurr=plm2xyz(this.state.plmWithDegreeAndOrder);
                    if ~iscell(likelihood)
                        fvalCurr=fvalCurr.*reshape(likelihood(z,[x1(:)';x2(:)';x3(:)']),size(latMesh));
                    else % if cell array, update for each likelihood and measurement
                        for i=1:numel(likelihood)
                            fvalCurr=fvalCurr.*reshape(likelihood{i}(z{i},[x1(:)';x2(:)';x3(:)']),size(latMesh));
                        end
                    end
                    warningSettings=warning('off','Normalization:notNormalized');
                    this.state=SphericalHarmonicDistributionComplex.fromPlmWithDegreeAndOrder(xyz2plm(fvalCurr),'identity');
                    warning(warningSettings);
                case 'sqrt'
                    likelihoodShd=SphericalHarmonicDistributionComplex.fromFunctionFast(@(x)likelihood(z,x),size(this.state.coeffMat,1)-1,'sqrt');
                    this.setState(this.state.multiply(likelihoodShd));
                otherwise
                    error('Transformation unsupported');
            end
        end
    end
    
end