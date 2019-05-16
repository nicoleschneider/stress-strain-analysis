% process image of necked specimen to determine strain
% Nicole Schneider
% Fall 2017 to Spring 2019

%% Scaling

scale = imread('.01mm.jpg');
%scale = imread('2018_scale.tif');
figure; imshow(scale);
%scale = imcrop(scale,[0 1000 3000 500]);
%%figure; imshow(scale);
%scale = imrotate(scale, -2.5); %%FIX ME HARDCODE
bwScale = im2bw(scale);  % convert scale image to black and white
%figure; imshow(bwScale)
scaleProfile = sum(bwScale);
%figure; plot(scaleProfile);

%%% from scale we have that 40 pixels is 1 mm %%%% (for 2018_Scale.tif)

format long;

%% code to automatically find scale conversion using spacing between peaks in profile of ruler
% %% Apply smooting
% w = gausswin(25);
% y = conv(scaleProfile,w);
% scaleDeriv = diff(y);
% 
% peaks = [];
% for i = 1:length(scaleDeriv)
%     if scaleDeriv(i) > 0 && (scaleDeriv(i+1) < 0 || scaleDeriv(i+1) == 0)
%         peaks = [peaks, i];
%     end
% end
% 
% peakSpacing = [];
% for i = 1 : length(peaks)-1
%     peakSpacing = [peakSpacing, peaks(i+1) - peaks(i)];
% end
% 
scaleFactor = 22.09; %mean(peakSpacing) % pixels per mm (for .01mm.jpg file)

%%
% choose an image to use by uncommenting one of the ones below or adding your own in this format
%filename = '2018.tif';
%filename = '171117 DSC_0284.tif';
%filename = 'DSC_0176_01(Round2).tif';
filename = 'CUspecimen1.tif';
img = imread(filename);
%img = imrotate(img, 90); %% FIX ME for 171117 DSC_0285 only, hardcoded to get rid of background garbage
%img = img(:,1850:length(img(1))); %% FIX ME hardcode cropping the image to get rid of garbage in backgrounud


%%
% new stuff 12/2/2018
% from: https://www.mathworks.com/help/images/detecting-a-cell-using-image-segmentation.html
% the code in this section is to do edge detection which helps find the specimen in the image

I = rgb2gray(img);
[~, threshold] = edge(I, 'sobel');
fudgeFactor = 1;
BWs = edge(I,'sobel', threshold * fudgeFactor);
figure, imshow(BWs), title('binary gradient mask');

se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);

BWsdil = imdilate(BWs, [se90 se0]);
figure, imshow(BWsdil), title('dilated gradient mask');

BWdfill = imfill(BWsdil, 'holes');
figure, imshow(BWdfill);
title('binary image with filled holes');

BWnobord = imclearborder(BWdfill, 4);
figure, imshow(BWnobord), title('cleared border image');

seD = strel('diamond',1);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);
figure, imshow(BWfinal), title('segmented image');

BWoutline = bwperim(BWfinal);
Segout = I; 
Segout(BWoutline) = 255; 
figure, imshow(Segout), title('outlined original image');

%%
% more edge detection stuff
%figure; imshow(img)
%BW = img; %im2bw(img);
img = rgb2gray(img);
BW = edge(img,'approxcanny', 0.5);
BW = imrotate(BW, 90);
figure; imshow(BW)

%%
% Connected component analysis
% this is to extract the specimen as one large group of pixels
BW = im2bw(img);
BW = imcrop(BW,[556 0 3260 4917]);  % hardcode crop to get rid of background garbage
L = bwlabel(imcomplement(BW));
RGB = label2rgb(L);
figure; imshow(RGB)

conComp = bwconncomp(BW);
numPixels = cellfun(@numel,conComp.PixelIdxList);
[area,idx] = max(numPixels);
area
%BW(conComp.PixelIdxList{idx}) = 0;

% show the profiles found
hProf = mean(BW, 2);
figure; plot(hProf)
hSum = sum(BW, 2);
figure; plot(hSum)

vProf = mean(BW);
figure; plot(vProf)
vSum = sum(BW);
figure; plot(vSum)


[m, n] = size(vSum);
idxMinArr = [];
for i = 1:n
    if vSum(i) == 0 % find columns that are all black
        idxMinArr = [idxMinArr i];
    end
end
%[min, idx] = min(vSum);
%vFlip = flip(vSum);

start = idxMinArr(1) % first all black column
last = idxMinArr(size(idxMinArr')) % last all black column
minWidth = last(1) - start(1)  % the minimum diameter of the specimen

%% 
vMean = mean(BW);
[m, n] = size(vMean);
idxMaxArr = [];
for i = 1:n-1
    if vMean(i) == 1 && vMean(i+1) ~= 1% find first column that is all white
        idxMaxArr = [idxMaxArr i];
    end
end
for i = idxMaxArr(1):n-1
    if vMean(i) ~= 1 && vMean(i+1) ~= 1% find first column that is all white
        idxMaxArr = [idxMaxArr i];
    end
end
maxStart = idxMaxArr(1) % first all black column
maxLast = idxMaxArr(size(idxMaxArr')) % last all black column
maxWidth = maxLast(1) - maxStart(1)


%% from 5/15/19
% this code finds the strain at every row of pixels
% you can run only the initial section to import the image and then this to calculate the strain
% without using edge detection or connected component analysis

BW_clean = imcrop(BW,[0 0 3260 4860]); % hardcode crop to get rid of garbage in background
imshow(BW_clean)

img_width = size(BW_clean(1,:))
img_height = size(BW_clean(:,1))

D_0 = 0
for X = 1 : img_height
    D_i_pixCount = 0;
    for I = 1 : img_width(2)
        if BW_clean(X,I) == 0
            D_i_pixCount = D_i_pixCount + 1;
        end
    end
    if D_i_pixCount > D_0
        D_0 = D_i_pixCount
    end
end
D_0

%%
% output epsilon values to a file
fileID = fopen('out.txt','w');
fprintf(fileID,'%6s %12s\r\n','Pixel Row','strain');


D_0 = 2143  % FOUND AFTER RUNNING CODE ONCE
for X = 1 : img_height
    D_i_pixCount = 0;
    for I = 1 : img_width(2)
        if BW_clean(X,I) == 0
            D_i_pixCount = D_i_pixCount + 1;
        end
    end
    
    epsilon = (D_0^2)/(D_i_pixCount^2) - 1
    A = [X; epsilon];
    fprintf(fileID,'%6.0f %22.15f\r\n', A);
end

fclose(fileID);

%%
% calculations using scale factors calculated earlier
% scaledMinWidth = minWidth / scaleFactor;
% scaledMinWidth
% 
% scaledMaxWidth = maxWidth / scaleFactor;
% scaledMaxWidth
