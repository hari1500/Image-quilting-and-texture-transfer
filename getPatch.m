function [xi,xj] = getPatch(block_h, block_v, texture, tolerance, o, blocksize, flag)
    [H,W] = size(texture);
    errors = zeros([H - blocksize+1, W - blocksize+1]);
    overlap_bh = block_h(:,blocksize-o+1:blocksize);
    overlap_bv = block_v(blocksize - o+1:blocksize,:);
    
    for i=1:H-blocksize+1
        for j=1:W-blocksize+1
            patch = texture(i:i+blocksize-1, j:j+blocksize-1);
            overlap_th = patch(:,blocksize - o+1:blocksize); 
            overlap_tv = patch(blocksize - o+1:blocksize,:);
            
            diff_h = (overlap_th - overlap_bh).^2;
            diff_v = (overlap_tv - overlap_bv).^2;
            if flag == 'h'
                errors(i,j) = sum(diff_h(:));
            elseif flag == 'v'
                errors(i,j) = sum(diff_v(:));
            else
                errors(i,j) = sum(diff_h(:)) + sum(diff_v(:));
            end
        end
    end
    
    minError = min(min(errors));
    [indi, indj] = find(errors < tolerance*minError);
    
    ind = randi(len(indi),1);
    xi = indi(ind);
    xj = indj(ind);
%     block = texture(xi:xi+blocksize-1, xj:xj+blocksize-1);
    
end