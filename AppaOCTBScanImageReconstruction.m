%% APPA OCT - B Scan Image RECONSTRUCTION  
function bScanImage = AppaOCTBScanImageReconstruction(RawDir,nImage)
%% OCT System Parameter Initialization
NPixel = 2048;%
lambda0 = 835E-9;
bandwidth = 50E-9;
lambdaMin = (lambda0 - bandwidth/2);
lambdaMax = (lambda0 + bandwidth/2);
dLambda = bandwidth / (NPixel-1);
lambda = lambdaMin:dLambda:lambdaMax;
c0 = 300000000;
omega0 =  2*pi*c0 / lambda0;
omega = 2*pi*c0 ./ lambda;
nAscanPerBscan = 2000;
%% Read Raw Image Data
for iImage = 1:nImage
    fid = fopen(fullfile(RawDir, sprintf('%05d.raw', iImage)), 'r');
    bScanCameraRawDataStream = fread(fid, 'int16', 'l');
    reshapedBScanCameraRawData = reshape(bScanCameraRawDataStream, [NPixel, nAscanPerBscan, 1]);% before 2048 2000
    fclose(fid);
    
    %% Size of Interference Fringe
    [nCameraPixels,nAscanPerBscan] = size(reshapedBScanCameraRawData);
    
    %% Type Conversion & subtract (2^16)/2
    typeConvertedBScanCameraRawData = double((reshapedBScanCameraRawData) - 32768);
    bScanData = typeConvertedBScanCameraRawData';
    
    %% Background Subtraction
    avgBScanSpectra = mean(bScanData, 1);
    avgBScan = ones(nAscanPerBscan,1) * avgBScanSpectra;
    dcRemovedbScanData = bScanData - avgBScan;
    
    %% Resampling or Interpolation
    preprocessing = load('mn_R&D_OEM_System.txt');
    aScanLength=1:1:nCameraPixels;
    for noAscan=1:nAscanPerBscan
        interpolatedBScanData(noAscan,:) = interp1(aScanLength,dcRemovedbScanData(noAscan,aScanLength),preprocessing(aScanLength),'spline','extrap');
    end
    
    %% Windowing the Data or Spectral Apodization
    window =  (nCameraPixels).*ones(nAscanPerBscan,1)*4096;
    filteredBScanData = interpolatedBScanData.*window;
    
    %% Dispersion Mode Selection
    AutoDispersion = 0;
    
    %% if Select Auto Dispersion means load Auto Dispersive Phase
    if (AutoDispersion == 1)
        dispersivePhase = load('AN_R&D_OCT_System_0010.txt');
    end
    
    %% Manual Dispersive Phase is Generated Using Taylor Series Expansion & Fix a2 & a3 Coefficients
    if AutoDispersion == 0
    a2 = 2500e-30;
    a3 = 1000e-45;
    lambda = lambdaMin:dLambda:lambdaMax;
    omega = 2*pi*c0 ./ lambda;
    dispersivePhase = -a2*((omega - omega0)).^2-a3*((omega - omega0)).^3;
    end
    
    %% Dispersion Compensation Part
    PhaseCorr = ones(nAscanPerBscan,1) * dispersivePhase;
    dispersionCompensatedBscan = filteredBScanData .* exp(-1i*PhaseCorr);
    dispersionCompensatedBscan = dispersionCompensatedBscan';
    
    %% IFFT & Exclude 25 Pixels that near by Zero Delay
    reconstructedBScanImageLinearScale = abs(ifft(dispersionCompensatedBscan));
    dbScaleBScanImage = 20 * log10(reconstructedBScanImageLinearScale);
    bScanImage(:,:,iImage) = dbScaleBScanImage(25:1024,:);
end
end