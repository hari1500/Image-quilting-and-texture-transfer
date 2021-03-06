function [leftCut,topCut] = getBoundary(leftBlock, topBlock, curBlock, o, blocksize, flag)
    % leftCut ------------------------------------------
    leftCut = zeros([blocksize,o]);
    topCut = zeros([o,blocksize]);
    if flag ~= 'v' 
        E_left = zeros([blocksize, o]);
        Back_left = zeros([blocksize, o]);
        E_left(1,:) = (leftBlock(1,blocksize-o+1:blocksize)-curBlock(1,1:o)).^2;
        for i = 2:blocksize
            [M,I] = min(E_left(i-1, 1:2));
            E_left(i,1) = (leftBlock(i,blocksize-o+1)-curBlock(i,1)).^2 + M;
            Back_left(i,1) = I; 
            for j = 2:o-1
                [M,I] = min(E_left(i-1, j-1:j+1));
                E_left(i,j) = (leftBlock(i,blocksize-o+j)-curBlock(i,j)).^2 + M;
                Back_left(i,j)= j-2+I;
            end
            [M,I] = min(E_left(i-1, j-1:j));
            E_left(i,o) = (leftBlock(i,blocksize)-curBlock(i,o)).^2 + M;
            Back_left(i,o) = j-2+I;
        end
        
        [~,I] = min(E_left(blocksize, :));
        leftCut(blocksize, I:o) = 1;
        I = Back_left(blocksize, I);
        for i = blocksize-1:-1:1
            leftCut(i,I:o) = 1;
            I = Back_left(i, I);
        end
    end    
    % topCut -----------------------------------------------
    if flag ~= 'h'
        E_top = zeros([o,blocksize]);
        Back_top = zeros([o,blocksize]);
        E_top(:,1) = (topBlock(blocksize-o+1:blocksize,1)-curBlock(1:o,1)).^2;
        for i = 2:blocksize
            [M,I] = min(E_top(1:2,i-1));
            E_top(1,i) = (topBlock(blocksize-o+1,i)-curBlock(1,i)).^2 + M;
            Back_top(1,i) = I; 
            for j = 2:o-1
                [M,I] = min(E_top(j-1:j+1,i-1));
                E_top(j,i) = (topBlock(blocksize-o+j,i)-curBlock(j,i)).^2 + M;
                Back_top(j,i)= j-2+I;
            end
            [M,I] = min(E_top(j-1:j,i-1));
            E_top(o,i) = (topBlock(blocksize,i)-curBlock(o,i)).^2 + M;
            Back_top(o,i) = j-2+I;
        end
        [~,I] = min(E_top(:,blocksize));
        topCut(I:o,blocksize) = 1;
        I = Back_top(I,blocksize);
        for i = blocksize-1:-1:1
            topCut(I:o,i) = 1;
            I = Back_top(I, i);
        end
    end    
end