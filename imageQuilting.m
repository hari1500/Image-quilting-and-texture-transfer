%%MainScript

%Input image
img = 'inputs/3.gif';

[original_img,map] = imread(img);
original_img = im2double(ind2rgb(original_img,map));
texture_img = rgb2gray(original_img);

[H,W,D] = size(texture_img);
imshow(texture_img);

% Parameters

block_size = 60;
overlap = 10;

% Ouput image
temp = block_size - overlap;
H_out = 2*temp*floor(H/temp) + overlap;
W_out = 2*temp*floor(W/temp) + overlap;

texture_out = zeros([H_out,W_out,D]);
