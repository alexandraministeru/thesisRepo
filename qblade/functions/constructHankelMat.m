function mat = constructHankelMat(data,i,s,Nbar)
%constructHankelMat(data,i,s,N) Function used to construct a block Hankel
%matrix using values from the provided data matrix.
%
% Input arguments
%----------------
% data    : matrix of dimension N-by-d, where N is the amount of data 
%           points, and d represents the number of input/output channels.
% i       : time index of the first input/ouput used for the Hankel matrix.
% s       : Hankel matrix block size.
% Nbar    : amount of data samples used to construct the block-hankel
%           matrix.
%
% Output arguments:
%------------------
% mat     : resulting block-Hankel matrix.
%==========================================================================

data = data.';
d = size(data,1);
mat = zeros(s*d,Nbar);

% construct first column
col = reshape(data(:,i:i+s-1),[],1);
mat(:,1) = col;

% construct next columns iteratively based on the previous columnn and
% given data

for idx= i+1:i+Nbar-1
    col = [col(d+1:end); data(:,s+idx-1)];
    mat(:,idx-i+1) = col;
end

end