function decompositionMatrix = doolittle(orginalMatrix)
    [height, width] = size(orginalMatrix);
    orginalMatrix(2:height,1) = orginalMatrix(2:height,1)/orginalMatrix(1);
    for j = 2:width
        for i = 2:height
            minIdx = min(i, j) - 1;
            orginalMatrix(i, j) = orginalMatrix(i, j) - ...
                orginalMatrix(i, 1:minIdx)*orginalMatrix(1:minIdx, j);
            if i>j
                orginalMatrix(i, j) = orginalMatrix(i, j) / orginalMatrix(j, j);
            end
        end
    end
    decompositionMatrix = orginalMatrix;
end
