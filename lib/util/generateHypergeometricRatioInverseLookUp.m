function generateHypergeometricRatioInverseLookUp(limError_dB, splines, MaxDim, MaxKappa)
    % This lookup table is neccessary for a maximum likelihood fit of a complex
    % Watson distribution. Since the generation of the lookup table is slow, the
    % lookup table is already provided.
    
    if nargin == 0
        limError_dB = -10;
        splines = 30;
        MaxDim = 6;
        MaxKappa = 70;
    end
    
    slicesDeltaKappa=zeros(splines, MaxDim-1);
    slicesKappa=zeros(splines, MaxDim-1);
    slicesLambda=zeros(splines, MaxDim-1);
    slicesFaktor=zeros(4, splines, MaxDim-1);
    
    for UpsDim=2:MaxDim
        % Search Spline Range for Limit Approx
        for n=10:100
            slicesKappa(splines, UpsDim-1)=n;
            slicesLambda(splines, UpsDim-1)=hypergeometricRatio(UpsDim, slicesKappa(splines, UpsDim-1));
            e=log10(abs(slicesLambda(splines, UpsDim-1) - (1-(UpsDim-1)/slicesKappa(splines, UpsDim-1))));
            if  e < limError_dB
                break;
            end
        end
        slicesLambda(1, UpsDim-1)=1/UpsDim;
        slicesLambda(splines, UpsDim-1)=1-(UpsDim-1)/slicesKappa(splines, UpsDim-1);
        
        % Divide in equally spaced spines in Kappa domain
        KappaStep=slicesKappa(splines, UpsDim-1)/(splines-1);
        for n=2:(splines-1)
            slicesKappa(n, UpsDim-1)=(n-1)*KappaStep;
            slicesLambda(n, UpsDim-1)=hypergeometricRatio(UpsDim, slicesKappa(n, UpsDim-1));
        end
        
        % Compute derivative
        for n=1:splines
            slicesDeltaKappa(n, UpsDim-1) ...
                = gradHypergeometricRatio(UpsDim, slicesKappa(n, UpsDim-1));
        end
        
        % Compute Cubic Spline Coefficients
        UpsilonMat=zeros(4,4);
        for n=1:splines-1
            UpsilonMat(1, :)=[slicesLambda(n, UpsDim-1)^3 slicesLambda(n, UpsDim-1)^2 slicesLambda(n, UpsDim-1) 1];
            UpsilonMat(2, :)=[slicesLambda(n+1, UpsDim-1)^3 slicesLambda(n+1, UpsDim-1)^2 slicesLambda(n+1, UpsDim-1) 1];
            UpsilonMat(3, :)=[3*slicesLambda(n, UpsDim-1)^2 2*slicesLambda(n, UpsDim-1) 1 0 ];
            UpsilonMat(4, :)=[3*slicesLambda(n+1, UpsDim-1)^2 2*slicesLambda(n+1, UpsDim-1) 1 0 ];
            UpsilonVec=[slicesKappa(n, UpsDim-1); slicesKappa(n+1, UpsDim-1); 1/slicesDeltaKappa(n, UpsDim-1); 1/slicesDeltaKappa(n+1, UpsDim-1)];
            slicesFaktor(:, n, UpsDim-1)=UpsilonMat\UpsilonVec;
        end
    end
    
    %Write to File
    fid = fopen('hypergeometricRatioInverseLookUp.h', 'w');
    fprintf(fid, '// Generated with generateHypergeometricRatioInverseLookUp(%g, %g, %g, %g);\n', limError_dB, splines,  MaxDim, MaxKappa );
    fprintf(fid, 'static unsigned int INVUPSILON__maxDimension=%i;\n', MaxDim);
    fprintf(fid, 'static unsigned int INVUPSILON__splines=%i;\n', splines);
    fprintf(fid, 'static double INVUPSILON__slicesLambda[]={\n\t');
    for c=1:((MaxDim-1)*splines)-1
        fprintf(fid, '%20.16e, ', slicesLambda(c));
        if mod(c,5)==0
            fprintf(fid, '\n\t');
        end
    end
    fprintf(fid, '%20.16e\n\t};\n', slicesLambda(((MaxDim-1)*splines)));
    fprintf(fid, 'static double INVUPSILON__slicesCoeff[]={\n\t');
    for c=1:((MaxDim-1)*splines*4)-1
        fprintf(fid, '%20.16e, ', slicesFaktor(c));
        if mod(c,5)==0
            fprintf(fid, '\n\t');
        end
    end
    fprintf(fid, '%20.16e\n\t};\n', slicesFaktor(((MaxDim-1)*splines*4)));
    fprintf(fid, '\n');
    fclose(fid);
end