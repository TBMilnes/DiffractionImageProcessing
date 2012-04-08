% Housekeeping
tic; format compact
warning('OFF', 'MATLAB:colon:nonIntegerIndex'); %Temporary...bad form...

% Declare some physical variables
sensorPixelPitch = 8.3*10^-3; %mm
aperturePixelPitch = 26*10^-3 %mm
Si = 300; %mm
wavelength = 500*10^-6 %mm

% Define image stack
apertureHeight = 768;
apertureStartWidth = 150; %pixels
apertureEndWidth = 200; %pixels
apertureStepWidth = 5; %pixels
apertureWidths = apertureStartWidth:apertureStepWidth:apertureEndWidth;
initialImageSize = [580,780];
imageTrim = [100,100,0,0];%[xLeft,xRight,yTop,yBottom]
numImages = length(apertureWidths);

% Declare other variable
gaussFigure = figure;
imageLineFigure = figure;
computedImageFigure = figure;

% Load images into stack and calculate imageSize
imageStack = zeros(initialImageSize(1)-sum(imageTrim(3:4)),initialImageSize(2)-sum(imageTrim(1:2)),numImages,'uint8');
for ii=1:numImages
    anImage = imread(sprintf('DiffractionLimit_2012-05-Apr-02-18/Aperture_%03dx%d.jpg',apertureWidths(ii),apertureHeight));
    anImage = anImage(imageTrim(3)+1:end-imageTrim(4),imageTrim(1)+1:end-imageTrim(2));
    imageStack(:,:,ii) = anImage;
end
imageSize = size(anImage); clear anImage

% Initialize computed image and optimization options
computedImage = zeros(imageSize,'uint8');
options = optimset('TolFun',10^-6);%10000000000*eps)

% Huge-ass 'for' loop to calculate actual intensities for EVERY image line...
for ii=240:imageSize(1)


    % Set up y vector -- this is the measured intensities at each pixel
    measuredIntensities = [];
    for jj=1:numImages
        measuredIntensities = [measuredIntensities,imageStack(ii,:,jj)];
    end
    

    % Set up C matrix -- the matrix of coefficients
    C = zeros(numImages*imageSize(2),imageSize(2));
    for jj=1:numImages    
        
%         % Set up gaussian profile to approximate Airy Disk
%         % WHAT IS THE RIGHT PROFILE??? GAUSSIAN APPEARS TO BE HALF THE
%         % WIDTH IT SHOULD BE BASED ON ABBE DIFFRACTION CALCULATIONS AND THE
%         % QUESTIONABLE FFT APPROACH
%         % I(q) = I0 * e^(-q^2/2w^2), w==gaussianWidth, q==radial distance
%         apertureWidth = apertureWidths(jj) * aperturePixelPitch;
%         fNumber = Si/apertureWidth; %f-Number: N=f/D
%         gaussianWidth = 0.43*wavelength*fNumber/0.5;
%         q = 0:sensorPixelPitch:1; q = [-fliplr(q),q(2:end)];
%         gaussianDist = exp(-(q.^2)/(2*gaussianWidth^2));
%         figure(gaussFigure); subplot(2,1,1); plot(q,gaussianDist)
%         % Trim gaussian distribution
%         gaussianDist = gaussianDist(find(gaussianDist>.01,1):end);
%         gaussianDist = gaussianDist(1:find(gaussianDist<.01,1)-1);
%         gaussianDist = gaussianDist/sum(gaussianDist);
%         subplot(2,1,2); plot(gaussianDist)
%         gaussDistHalfWidth = (length(gaussianDist)-1)/2

        % Get coefficients
        coefficients = ScaledDiffractionPatternGenerator(apertureWidths(jj),apertureHeight);
        coefficientsHalfWidth = (length(coefficients)-1)/2
        
        % Compute and store block of C for this aperture setting
        for kk=1:imageSize(2)
            %for ll = -gaussDistHalfWidth:gaussDistHalfWidth
            for ll = -coefficientsHalfWidth:coefficientsHalfWidth
                if kk+ll>0 && kk+ll<=imageSize(2)
                    %C((jj-1)*imageSize(2)+kk,kk+ll) = gaussianDist(ll+gaussDistHalfWidth+1);
                    C((jj-1)*imageSize(2)+kk,kk+ll) = coefficients(ll+coefficientsHalfWidth+1);
                end
            end
        end
        
        
    end

    % Run the optimization, using largest-aperture values as starting point
    [actualIntensities,RESNORM,RESIDUAL,EXITFLAG,output] = lsqlin(C,double(measuredIntensities),[],[],[],[],zeros(imageSize(2),1),Inf(imageSize(2),1),double(imageStack(ii,:,end)),options);
    output
    figure(imageLineFigure); subplot(2,1,1); plot(measuredIntensities(end-imageSize(2)+1:end)); axis tight
    subplot(2,1,2); plot(actualIntensities); axis tight
    max(actualIntensities)

    
    % Save results to computed image
    computedImage(ii,:) = actualIntensities;
    figure(computedImageFigure); imshow(computedImage);
    
    toc
    
    % Set up x vector -- vector of actual values which are to be solved for

    
    
end


% Finish up
toc