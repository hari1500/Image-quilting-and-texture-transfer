%%MainScript
tic;

% Input texture image
img1 = 'inputs/transfer/rice.png';
% Input target image
img2 = 'inputs/transfer/bill.png';

original_img = im2double(imread(img1));
target_img = im2double(imread(img2));

if(length(size(original_img)) ~= 3)
    original_img = repmat(original_img,[1 1 3]);
end

if(length(size(target_img)) ~= 3)
    target_img = repmat(target_img,[1 1 3]);
end

texture_img = rgb2gray(original_img);
target_img = rgb2gray(target_img);

[H,W,D] = size(original_img);

blocksize = 32;
o = round(blocksize/3);
tolerance = 1.1;

% Ouput image
[H_out,W_out] = size(target_img);
net_patch = blocksize - o;
H_out = net_patch*floor((H_out-o)/net_patch) + o;
W_out = net_patch*floor((W_out-o)/net_patch) + o;

texture_out = zeros([H_out,W_out]);
texture_out3D = zeros([H_out, W_out,D]);

nIter = 2;
% alpha = 0.4;
for k = 1:nIter
    alpha = 0.8*((k-1)/(nIter-1))+0.1;
    foo = waitbar(0,'Texture Transfer');
    for i=1:net_patch:H_out-blocksize+1
        for j=1:net_patch:W_out-blocksize+1
            j_prev = j - net_patch;
            i_prev = i - net_patch;

            i_inds = i:i+blocksize-1;
            j_inds = j:j+blocksize-1;
            b_inds = blocksize-o+1:blocksize;

            if i==1 && j == 1
                xind = 1;
	            yind = 1;
                texture_out(i_inds,j_inds) = texture_img(xind:xind+blocksize-1,yind:yind+blocksize-1);
                texture_out3D(i_inds,j_inds,:) = original_img(xind:xind+blocksize-1,yind:yind+blocksize-1,:);

            elseif i == 1
                block_h = texture_out(i_inds,j_prev:j_prev+blocksize-1);
                block_h3D = texture_out3D(i_inds,j_prev:j_prev+blocksize-1,:);
                
                targetPatch = target_img(i_inds,j_inds);
                [xind, yind] = getPatchTransfer(block_h,-1,texture_img,tolerance,o,blocksize,'h', targetPatch, alpha);

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
                
                targetPatch = target_img(i_inds,j_inds);
                [xind, yind] = getPatchTransfer(-1,block_v,texture_img,tolerance,o,blocksize,'v', targetPatch, alpha);

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
                
                targetPatch = target_img(i_inds,j_inds);
                [xind, yind] = getPatchTransfer(block_h,block_v,texture_img,tolerance,o,blocksize,'m', targetPatch, alpha);

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
        waitbar(i/((H_out-blocksize+1)),foo);
    end
    close(foo);
    
	figure;
    imshow(texture_out3D);
    [filepath,name,ext] = fileparts(img2);
	saveas(gcf,strcat('outputs/transfer/', strcat(name,int2str(k)), '.png'));
    
    target_img = texture_out;
    blocksize = round(blocksize/3);
    o = round(blocksize/3);
    net_patch = blocksize - o;
	H_out = net_patch*floor((H_out-o)/net_patch) + o;
	W_out = net_patch*floor((W_out-o)/net_patch) + o;

	texture_out = zeros([H_out,W_out]);
	texture_out3D = zeros([H_out, W_out,D]);
%     alpha = 0.8*((k)/(nIter-1))+0.1;
end

toc;