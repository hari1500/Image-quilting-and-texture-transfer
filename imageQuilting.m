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
tolerance = 

% Ouput image
temp = block_size - overlap;
H_out = 2*temp*floor(H/temp) + overlap;
W_out = 2*temp*floor(W/temp) + overlap;

texture_out = zeros([H_out,W_out]);
texture_out3D = zeros([H_out, W_out,D]);

for i=1:H_out
    for j=1:W_out
        
        if i==1 && j == 1
            xind = randi(H-block_size,1);
            yind = randi(W-block_size,1);
            texture_out(1,1) = texture_img(xind,yind);
            texture_out3D(1,1,:) = original_img(xind,yind,:);
        end
    end
end

            