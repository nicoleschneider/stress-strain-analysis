% process image of necked specimen to determine width
% Nicole Schneider

%% Scaling

%scale = imread('.01mm.jpg');
scale = imread('2018_scale.tif');
%figure; imshow(scale);
scale = imcrop(scale,[0 1000 3000 500]);
%%figure; imshow(scale);
%scale = imrotate(scale, -2.5); %%FIX ME HARDCODE
bwScale = im2bw(scale);
%figure; imshow(bwScale)
scaleProfile = sum(bwScale);
%figure; plot(scaleProfile);

%%% from scale we have that 40 pixels is 1 mm %%%%

format long;

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
scaleFactor = 40; %mean(peakSpacing) % pixels per mm

%%
filename = '2018.tif';
%filename = '171117 DSC_0284.tif';
%filename = 'DSC_0176_01(Round2).tif';
%filename = 'CUspecimen1.tif';
img = imread(filename);
%img = imrotate(img, 90); %% FIX ME for 171117 DSC_0285 only
%img = img(:,1850:length(img(1))); %% FIX ME same as above
%%
% new stuff 12/2/2018
% from: https://www.mathworks.com/help/images/detecting-a-cell-using-image-segmentation.html

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

%figure; imshow(img)
%BW = img; %im2bw(img);
img = rgb2gray(img);
BW = edge(img,'approxcanny', 0.5);
BW = imrotate(BW, 90);
figure; imshow(BW)


% Connected component analysis
L = bwlabel(imcomplement(BW));
RGB = label2rgb(L);
figure; imshow(RGB)

conComp = bwconncomp(BW);
numPixels = cellfun(@numel,conComp.PixelIdxList);
[area,idx] = max(numPixels);
area;
%BW(conComp.PixelIdxList{idx}) = 0;

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
minWidth = last(1) - start(1)

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


%%

% scaledMinWidth = minWidth / scaleFactor;
% scaledMinWidth
% 
% scaledMaxWidth = maxWidth / scaleFactor;
% scaledMaxWidth
