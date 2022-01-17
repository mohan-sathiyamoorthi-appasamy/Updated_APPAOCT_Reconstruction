%% B-Scan Image Reconstruction with Noise Flattened Process
nImage = 1;
nAscanPerBscan = 2000;
%% Raw Reference or No Eye B-Scan Image Directory
RawDirRef = '.\RawNoEyeData_test';
noEyeBScanImage = AppaOCTBScanImageReconstruction(RawDirRef,nImage);

%% Noise Floor Determination
windowSize=5;
noiseAvg=mean(noEyeBScanImage,2);% dBScalePSF is the reconstructed B-scan (after logarithm step) of a reference scan
noiseFloor=filter(ones(1,windowSize)/windowSize,1,noiseAvg);
noiseFloor(1:4)=noiseFloor(5);
noiseFlatten = max(noiseFloor) - noiseAvg;
noEyeNoiseFloor = noiseFlatten(:,1);
noEyeNoiseFloorBscan = repmat(noEyeNoiseFloor,1,nAscanPerBscan);

%% Raw Human Eye B-Scan Image Directory
RawDirHuman = '.\RawHumanEyeData_test';
bScanImage = AppaOCTBScanImageReconstruction(RawDirHuman,nImage);

%% Add Noise Floor to B-Scan Image
FinalImage = bScanImage + noEyeNoiseFloorBscan;

%% Image Normalization
baseImageFloor = 44;%44
dynamicRange = 55;%55
grayScaleUnitsPerDB = 255/dynamicRange;
normImage = uint8(grayScaleUnitsPerDB*(FinalImage-baseImageFloor));
%figure(4),imagesc(normImage,[20 210]);colormap(gray);
%% Display Image & Plot
figure(1),imshow(FinalImage,[min(FinalImage(:)),max(FinalImage(:))]);
figure(2),plot(noEyeNoiseFloor),xlabel('Pixels'),ylabel('Intensity in dB'),title('Noise Floor in Reference Data');
figure(3),plot(mean(bScanImage,2),'r'),hold on,plot(mean(FinalImage,2),'b'),xlabel('Pixels'),ylabel('Intensity in dB'),title('Comparison of Intensity Profile of Before & After Noise Flatten'),legend('Before Noise Flatten','After Noise Flatten');
