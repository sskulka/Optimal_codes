function x = exclude(x,indx,rc)
% exclude Exclude array elements
% Y = exclude(X,indx,rc) takes an input X, which can be eithr an array or
% 2D matrix.
% 
% indx is the index of the row(s) or column(s) or elements to omit from the
% matrix or array.
% 
% rc is used to define whether the rows or columns to be excluded. For row
% exclusion use: 'row' and for column exlusion use: 'col' or 'column'. If
% nothing is defined then the default 'col' is set automatically. However,
% rc will automatically be set to 'col' or 'row' if X is a row or column
% vector, respectively.
% 
% Example 1
% x = 1:10
% y1 = exclude(x,[1 5])
% y2 = exclude(x',5)
% 
% Example 2
% X = magic(3)
% y1 = exclude(X,2,'row')
% y2 = exclude(X,2,'col')
% 
% Copyright @ Md Shoaibur Rahman (shaoibur@bcm.edu)

if nargin < 3
    rc = 'col';
end
if size(x,1)==1
    rc = 'col';
elseif size(x,2)==1
    rc = 'row';
end
if strcmp(rc,'col') || strcmp(rc,'column')
    x(:,indx) = [];
elseif strcmp(rc,'row')
    x(indx,:) = [];
end