function averageMatrix = averMeanFields(allMatrices)
    % averMeanFields - Compute the average matrix from a cell array of 3D matrices.
    %
    % Syntax:
    %   averageMatrix = computeAverageMatrix(allMatrices)
    %
    % Inputs:
    %   allMatrices - Cell array containing 3D matrices (each matrix can have NaN values).
    %
    % Output:
    %   averageMatrix - 3D matrix representing the average value of each element
    %                   over all input matrices, ignoring NaN values.
    %

    % Convert the cell array of matrices into a 4D array
    stackedMatrices = cat(4, allMatrices{:});

    % Compute the average value of each element over all matrices, ignoring NaN
    averageMatrix = nanmean(stackedMatrices, 4);
end
