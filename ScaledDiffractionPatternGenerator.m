% From slide 13: https://docs.google.com/viewer?a=v&q=cache:th9OZFCV_FsJ:ece661web.groups.et.byu.net/notes/complex_apertures.ppt+&hl=en&gl=us&pid=bl&srcid=ADGEESgxIHxPcA9G7PBjgsfuZqp5OoHJGSk7QbJ4DErEX4tdfeoEREaBjg-MNDHNJ1mblXN0ApeHozc0rFf5ywj90EXTybwpTsCN_d1RiPac01LaN3K8CjRUyycm9PY6D76dawWBcw8Y&sig=AHIEtbQKpeguoO6dh-zqpLW5U-Gs5vzTUg
% QUESTIONS: Should I be squared?  Should it be abs()'d or real()'d?

function [coefficients] = ScaledDiffractionPatternGenerator(apPixWidth,apPixHeight)

% Slit dimensions in pixels
%apPixWidth  = 200;
%apPixHeight = 768;

% Define physical constants
aperturePlateResolution = [768,1024];
aperturePixelPitch = 26*10^-6;
sensorPixelPitch = 8.3*10^-6;
lambda = 525e-9; %Green
f = 310*10^-3;
paddingMultiple = 3;
%W = apPixWidth*aperturePixelPitch; %Width of aperture image
W = aperturePlateResolution(2)*paddingMultiple*aperturePixelPitch;

% Generate aperture matrix
A=zeros(paddingMultiple*aperturePlateResolution);
[M,N] = size(A);
apSize = size(A);
A(apSize(1)/2+1-apPixHeight/2:apSize(1)/2+apPixHeight/2,...
    apSize(2)/2+1-apPixWidth/2:apSize(2)/2+apPixWidth/2)=1;
%figure; imagesc(A);

% FFT
I = abs(fftshift(fft2(A))).^2;
x = linspace(-lambda*f*N/(2*W), lambda*f*N/(2*W), N);
y = linspace(-lambda*f*M/(2*W), lambda*f*M/(2*W), M);

%Plotting
%figure; pcolor(x, y, I)
trim = [300,428]+(paddingMultiple-1)/2*aperturePlateResolution;
%figure; surf(x, y, I); shading('interp'); axis tight;
%figure; surf(x(trim(2)+1:end-trim(2)), y(trim(1)+1:end-trim(1)),...
%    I(trim(1)+1:end-trim(1),trim(2)+1:end-trim(2))); axis tight; %shading('interp');

% Extract central profile
centralProfile = I(aperturePlateResolution(1)*paddingMultiple/2+1,:);
%figure; plot(x,centralProfile)

% Calculate minima and truncate profile
numLobes = 2; %2 == central lobe and first ring
minima = [false, (centralProfile(2:end-1)<centralProfile(3:end))...
    & (centralProfile(2:end-1)<centralProfile(1:end-2)), false];
%hold on; plot(x,minima*max(centralProfile)/10,'r')
minimaIndices = find(minima==1);
centralProfile = centralProfile(minimaIndices(sum(minima)/2-numLobes+1):...
    minimaIndices(sum(minima)/2+numLobes));

% Calculate centered xCoords for truncated profile
profileWidth = x(minimaIndices(sum(minima)/2+numLobes)) - ...
    x(minimaIndices(sum(minima)/2-numLobes+1));
centralProfileXCoords = linspace(-profileWidth/2,profileWidth/2,length(centralProfile));
%figure; plot(centralProfileXCoords,centralProfile); axis tight

% 
% centralProfileXCoords(end)-centralProfileXCoords(1)
% centralProfileXCoords = x(minimaIndices(sum(minima)/2-numLobes+1):...
%     minimaIndices(sum(minima)/2+numLobes));

% Sanity check
if max(centralProfile) ~= max(max(I))
    error() %Generic error
end

% Calculate number of coefficients for fitting
numCoefficients = floor(profileWidth/sensorPixelPitch);
if (floor(numCoefficients/2) == numCoefficients/2) %We want an odd number
    numCoefficients = numCoefficients + 1;
end

% Calculate coefficient values
%profileWidth = (length(centralProfile)-1)/paddingMultiple*aperturePixelPitch;
%centralProfileXCoords = linspace(-profileWidth/2,profileWidth/2,length(centralProfile));
coefficientWidth = numCoefficients*sensorPixelPitch;
coefficientCenters = linspace(-(coefficientWidth-sensorPixelPitch)/2,...
    (coefficientWidth-sensorPixelPitch)/2,numCoefficients);
coefficientCenterValues = interp1(centralProfileXCoords,centralProfile,coefficientCenters);
%figure; plot(centralProfileXCoords,centralProfile,'-',centralProfileXCoords,centralProfile,'+',coefficientCenters,coefficientCenterValues,'rx')

% Normalize coefficients
coefficients = coefficientCenterValues/sum(coefficientCenterValues);
%figure; plot(coefficients); axis tight
return 



