function reducedMat=reduce2DMatSum(iniMat,reductionAmount)

% Get the subscripts for each element of the input matrix
[subsAi, subsAj] = ndgrid(1:size(iniMat,1), 1:size(iniMat,2));
% Compute the corresponding output subscripts
subsBi = 1 + floor((subsAi-1)/reductionAmount);
subsBj = 1 + floor((subsAj-1)/reductionAmount);
% Carry out the summation
 reducedMat = accumarray([subsBi(:) subsBj(:)], iniMat(:));   % b is reduced matrix
% See if the overall summation was correct
% disp(abs(sum(iniMat(:)) - sum(reducedMat(:))));    % should be small
end