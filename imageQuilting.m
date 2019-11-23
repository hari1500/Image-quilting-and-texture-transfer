%%MainScript

%Input image
img = 'inputs/quilting/S2.gif';

[original_img,map] = imread(img);
original_img = im2double(ind2rgb(original_img,map));
texture_img = rgb2gray(original_img);

% original_img = original_img(1:60, 1:60);
[H,W,D] = size(original_img);

imshow(texture_img);

figure;
imshow(original_img);
% Parameters

blocksize = 90;
o = 20;
tolerance = 1.1;

% Ouput image
temp = blocksize - o;
H_out = 2*temp*floor(H/temp) + o;
W_out = 2*temp*floor(W/temp) + o;

texture_out = zeros([H_out,W_out]);
texture_out3D = zeros([H_out, W_out,D]);

net_patch = blocksize - o;
for i=1:net_patch:H_out-blocksize+1
    for j=1:net_patch:W_out-blocksize+1
        
        j_prev = j - net_patch;
        i_prev = i - net_patch;
        
        i_inds = i:i+blocksize-1;
        j_inds = j:j+blocksize-1;
        b_inds = blocksize-o+1:blocksize;
        
        if i==1 && j == 1
            xind = randi(H-blocksize,1);
            yind = randi(W-blocksize,1);
            texture_out(i_inds,j_inds) = texture_img(xind:xind+blocksize-1,yind:yind+blocksize-1);
            texture_out3D(i_inds,j_inds,:) = original_img(xind:xind+blocksize-1,yind:yind+blocksize-1,:);
        
        elseif i == 1
            block_h = texture_out(i_inds,j_prev:j_prev+blocksize-1);
            block_h3D = texture_out3D(i_inds,j_prev:j_prev+blocksize-1,:);
            [xind, yind] = getPatch(block_h,-1,texture_img,tolerance,o,blocksize,'h');
            
            curr_patch = texture_img(xind:xind+blocksize-1,yind:yind+blocksize-1);
            curr_patch3D = original_img(xind:xind+blocksize-1,yind:yind+blocksize-1,:);
            
            [mask_h,mask_v] = getBoundary(block_h,-1,curr_patch,o,blocksize,'h');
            mask_h3D = repmat(mask_h,[1,1,3]);
            
            curr_patch(:,1:o) = block_h(:,b_inds).*(1-mask_h) + curr_patch(:,1:o).*mask_h;
            curr_patch3D(:,1:o,:) = block_h3D(:,b_inds,:).*(1-mask_h3D) + curr_patch3D(:,1:o,:).*mask_h3D;
            
            texture_out(i_inds,j_inds) = curr_patch;
            texture_out3D(i_inds,j_inds,:) = curr_patch3D;
            
        elseif j == 1
            block_v = texture_out(i_prev:i_prev+blocksize-1,j_inds);
            block_v3D = texture_out3D(i_prev:i_prev+blocksize-1,j_inds,:);
            [xind, yind] = getPatch(-1,block_v,texture_img,tolerance,o,blocksize,'v');
            
            curr_patch = texture_img(xind:xind+blocksize-1,yind:yind+blocksize-1);
            curr_patch3D = original_img(xind:xind+blocksize-1,yind:yind+blocksize-1,:);
            
            [mask_h,mask_v] = getBoundary(-1,block_v,curr_patch,o,blocksize,'v');
            mask_v3D = repmat(mask_v,[1,1,3]);
            
            curr_patch(1:o,:) = block_v(b_inds,:).*(1-mask_v) + curr_patch(1:o,:).*mask_v;
            curr_patch3D(1:o,:,:) = block_v3D(b_inds,:,:).*(1-mask_v3D) + curr_patch3D(1:o,:,:).*mask_v3D;
            
            texture_out(i_inds,j_inds) = curr_patch;
            texture_out3D(i_inds,j_inds,:) = curr_patch3D;
            
        else
            block_h = texture_out(i_inds,j_prev:j_prev+blocksize-1);
            block_h3D = texture_out3D(i_inds,j_prev:j_prev+blocksize-1,:);
            block_v = texture_out(i_prev:i_prev+blocksize-1,j_inds);
            block_v3D = texture_out3D(i_prev:i_prev+blocksize-1,j_inds,:);
            [xind, yind] = getPatch(block_h,block_v,texture_img,tolerance,o,blocksize,'m');
            
            curr_patch = texture_img(xind:xind+blocksize-1,yind:yind+blocksize-1);
            curr_patch3D = original_img(xind:xind+blocksize-1,yind:yind+blocksize-1,:);
            
            [mask_h,mask_v] = getBoundary(block_h,block_v,curr_patch,o,blocksize,'m');
            mask_h3D = repmat(mask_h,[1,1,3]);
            mask_v3D = repmat(mask_v,[1,1,3]);
            
            curr_patch(:,1:o) = block_h(:,b_inds).*(1-mask_h) + curr_patch(:,1:o).*mask_h;
            curr_patch3D(:,1:o,:) = block_h3D(:,b_inds,:).*(1-mask_h3D) + curr_patch3D(:,1:o,:).*mask_h3D;
            curr_patch(1:o,:) = block_v(b_inds,:).*(1-mask_v) + curr_patch(1:o,:).*mask_v;
            curr_patch3D(1:o,:,:) = block_v3D(b_inds,:,:).*(1-mask_v3D) + curr_patch3D(1:o,:,:).*mask_v3D;
            
            texture_out(i_inds,j_inds) = curr_patch;
            texture_out3D(i_inds,j_inds,:) = curr_patch3D;
        end
        
    end
end

figure;
imshow(texture_out);
figure;
imshow(texture_out3D);
