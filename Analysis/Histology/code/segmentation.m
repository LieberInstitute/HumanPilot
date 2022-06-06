% Read in histology image with imread function %
imgRGB = imread('/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/10X/20190925_JHU_Lieber_HumanBrain_LP/');

% convert image datatype to double %
imgRGB = im2double(imgRGB);

% Applied gaussian smoothening and increased the contrast in image to get rid of the noise in our images, you might not need this step if your images are high quality %
imgRGB_smooth = imgaussfilt(imgRGB,4);
imgRGB_smooth_adj = imadjust(imgRGB_smooth, [.2 .3 0; .6 .7 1],[]); % These numbers are manually adjusted for each image until desired contrast is obtained %

%Applying Kmeans color based segmentation on the above resulting image
%following the technique mentioned in https://www.mathworks.com/help/images/color-based-segmentation-using-k-means-clustering.html%
he = imgRGB_smooth_adjus;
lab_he = rgb2lab(he); %convert image from rgb color space to lab color space%
ab = lab_he(:,:,2:3);
ab = im2single(ab);
nColors = 5; %user defined number, number of colors visually seen by the user in a image% 
pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3); %apply Kmeans

%resulting clusters from Kmeans
mask1 = pixel_labels==1; %mask* is a binary segmented image per color
cluster1 = he .* uint8(mask1); %cluster* is a colored segmented image per color
mask2 = pixel_labels==2;
cluster2 = he .* uint8(mask2);
mask3 = pixel_labels==3;
cluster3 = he .* uint8(mask3);
mask4 = pixel_labels==4;      
cluster4 = he .* uint8(mask4);
mask5 = pixel_labels==5;      
cluster5 = he .* uint8(mask5);

%one of these five cluster/binary images is the segmented nuclei image
%Usually nuclei in histology images have a distinct dark color and the
%nuclei are clearly segmented from the background

imwrite(cluster1,'/dcl01/lieber/ajaffe/Maddy/RNAscope/Histology/10Ximages/imgRGB_smooth_adjus_cluster1.png')
imwrite(cluster2,'/dcl01/lieber/ajaffe/Maddy/RNAscope/Histology/10Ximages/imgRGB_smooth_adjus_cluster2.png')
imwrite(cluster3,'/dcl01/lieber/ajaffe/Maddy/RNAscope/Histology/10Ximages/imgRGB_smooth_adjus_cluster3.png')
imwrite(cluster4,'/dcl01/lieber/ajaffe/Maddy/RNAscope/Histology/10Ximages/imgRGB_smooth_adjus_cluster4.png') 
imwrite(cluster5,'/dcl01/lieber/ajaffe/Maddy/RNAscope/Histology/10Ximages/imgRGB_smooth_adjus_cluster5.png') 

