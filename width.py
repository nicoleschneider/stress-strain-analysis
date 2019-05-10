import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import cv2

import PIL
from PIL import Image
import imageio



# process image of necked specimen to determine width
# Nicole Schneider

## Scaling

scale = matplotlib.pyplot.imread('.01mm.jpg')   #read in the image of a ruler
#scale = matplotlib.pyplot.imread('2018_scale.tif');   #read in different image of ruler
plt.figure()          #display image of ruler
plt.imshow(scale)
plt.show
scale = matplotlib.pyplot.imcrop(scale,[0 1000 3000 500])    #FIX ME HARDCODED, hard coded cropped part. 
#In future will have to automatically detect ruler part of image. 
plt.figure()          #Display cropped image for testing
plt.imshow(scale)
plt.show

##gray = scale.convert('L')        #convert colored image(scale), save resulting bw image to bwscal
##bw = gray.point(lambda x: 0 if x<128 else 255, '1')
##bw.save("scale_bw.png")

#scale = imrotate(scale, -2.5); ##FIX ME HARDCODE, Purpose: rotate image to unskew the ruler. 
#This will need to be automated instead of the -2.5 degree hard code.
bwScale = im2bw(scale)       #convert colored image(scale), save resulting bw image to bwscale
plt.figure()          #Display bwimage image for testing
plt.imshow(bwScale)
plt.show
scaleProfile = sum(bwScale)  #  crate a profile by summing pixel values
#figure; plot(scaleProfile);
###FIX ME, MISSING CODE. To calculate the actual scale based on the bw scale instead of the hardcoded 40 pixels.
%%% from scale we have that 40 pixels is 1 mm %%%%
#----------------------------------------------------------------------------

% %% Apply smoothing
  blur = []
  blur = cv2.GaussianBlur(bwscale,(25,25),0) #create a guassian blurred object
% w = gausswin(25)  # create a guassian smoothing object
% y = conv(scaleProfile,w)    # convolve scale profile with gaussian smoothing
% scaleDeriv = diff(y)    # find derivative of smoothed profile
% 
% peaks = []  # to hold relative maxima of scale profile
% for i = 1:length(scaleDeriv)    # loop over profile
%     if scaleDeriv(i) > 0 && (scaleDeriv(i+1) < 0 || scaleDeriv(i+1) == 0)   # if point is a peak
%         peaks = [peaks, i]    # add index to peaks array
%     end
% end
% 
% peakSpacing = []   # to hold values of differences between peaks
% for i = 1 : length(peaks)-1    # loop over peaks array but dont run over end (stay one before end)
%     peakSpacing = [peakSpacing, peaks(i+1) - peaks(i)]   # add difference value to peakSpacing array
% end
% 
scaleFactor = 40; %mean(peakSpacing) % pixels per mm   # hardcoded as 40, but should be calculated using mean(peakSpacing)

%%   # list of different specimen images to test on
filename = '2018.tif'
%filename = '171117 DSC_0284.tif'
%filename = 'DSC_0176_01(Round2).tif'
%filename = 'CUspecimen1.tif'
img = imread(filename)    # read in image from file
%img = imrotate(img, 90) %% FIX ME for 171117 DSC_0285 only   # this image was skewed so rotate it, should be automated in future.
%img = img(:,1850:length(img(1))) %% FIX ME same as above    # hardcoded cropping that should be fixed in future
%%
% new stuff 12/2/2018
% from: https://www.mathworks.com/help/images/detecting-a-cell-using-image-segmentation.html

    # Edge detection code
I = rgb2gray(img)   # convert image to grayscale
[~, threshold] = edge(I, 'sobel')    # sobel type of edge detection
fudgeFactor = 1
BWs = edge(I,'sobel', threshold * fudgeFactor)    # run edge detector
figure, imshow(BWs), title('binary gradient mask')    # display result containing the edges it found

# list of possible sturucturing elemnts to choose from
se90 = strel('line', 3, 90)  
se0 = strel('line', 3, 0)

BWsdil = imdilate(BWs, [se90 se0])     # dilate image using structuring elemnt (metamorphic operations?)
figure, imshow(BWsdil), title('dilated gradient mask')   # display results

BWdfill = imfill(BWsdil, 'holes')    # fill in spaces
figure, imshow(BWdfill)
title('binary image with filled holes')

BWnobord = imclearborder(BWdfill, 4)    # get rid of borders
figure, imshow(BWnobord), title('cleared border image')

seD = strel('diamond',1)
BWfinal = imerode(BWnobord,seD)   # do erosion operation
BWfinal = imerode(BWfinal,seD)
figure, imshow(BWfinal), title('segmented image')

BWoutline = bwperim(BWfinal)
Segout = I
Segout(BWoutline) = 255
figure, imshow(Segout), title('outlined original image')

%%
%figure; imshow(img)
%BW = img; %im2bw(img)
img = rgb2gray(img)
BW = edge(img,'approxcanny', 0.5)
BW = imrotate(BW, 90)
figure; imshow(BW)


% Connected component analysis
L = bwlabel(imcomplement(BW))
RGB = label2rgb(L)
figure; imshow(RGB)

conComp = bwconncomp(BW)
numPixels = cellfun(@numel,conComp.PixelIdxList)
[area,idx] = max(numPixels)
area
%BW(conComp.PixelIdxList{idx}) = 0

hProf = mean(BW, 2)  # horizontal profile
figure; plot(hProf)
hSum = sum(BW, 2)
figure; plot(hSum)

vProf = mean(BW)   # vertical profile
figure; plot(vProf)
vSum = sum(BW)
figure; plot(vSum)


[m, n] = size(vSum)
idxMinArr = []
for i = 1:n
    if vSum(i) == 0 # find columns that are all black
        idxMinArr = [idxMinArr i]
    end
end
%[min, idx] = min(vSum)
%vFlip = flip(vSum)

start = idxMinArr(1) # find index of first all black column (ie column with min average value)
last = idxMinArr(size(idxMinArr')) # last all black column
minWidth = last(1) - start(1)

%% 
vMean = mean(BW)
[m, n] = size(vMean)
idxMaxArr = []
for i = 1:n-1
    if vMean(i) == 1 && vMean(i+1) ~= 1 # find first column that is all white
        idxMaxArr = [idxMaxArr i]
    end
end
for i = idxMaxArr(1):n-1
    if vMean(i) ~= 1 && vMean(i+1) ~= 1  # find first column that is all white
        idxMaxArr = [idxMaxArr i]
    end
end
maxStart = idxMaxArr(1) # first all black column
maxLast = idxMaxArr(size(idxMaxArr')) # last all black column
maxWidth = maxLast(1) - maxStart(1)


%%
% scaledMinWidth = minWidth / scaleFactor
% scaledMinWidth
% 
% scaledMaxWidth = maxWidth / scaleFactor
% scaledMaxWidth

