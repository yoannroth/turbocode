function [ B ] = myde2bi( varargin )
%MYDE2BI Convert decimal input (in double) to binary
%   B = myde2bi(D, n)  where D is the decimal input matrix and n the number
%   of digits (columns) for the binary words (right-msb orientation is 
%   assumed)
%
%   B = myde2bi(D, n, flag)    where B is the decimal input matrix, n the 
%   number of bits for the binary words and flag is equal to right-msb or 
%   left-msb orientation

% Arg check
switch length(varargin)
    case 2
        D       = varargin{1};
        n       = varargin{2};
        msbFlag = 'right-msb';
    case 3
        D       = varargin{1};
        n       = varargin{2};
        msbFlag = varargin{3};
    case {0, 1}
        error('Not enough input argument')
    otherwise
        error('Too many input arguments')
end

% D size check
if size(D, 1)<size(D, 2);
    D = D';
end
if size(D, 2)~=1
    error('Input decimal values must be a vector')
end

% n value check
if ~isa(n, 'double')
    error('Wrong type of argument for number of digit')
end
if max(D)>(2^n-1)
    error('Not enough bits to represent the maximum value')
end


% Base
base    = 2.^((1-n):0);
% Convertion
B       = rem(floor(D*base),2);
% Bits orientation
if strcmp(msbFlag, 'right-msb')
    B   = fliplr(B);
elseif strcmp(msbFlag, 'left-msb')
    % Nothing to do
else
    error('Wrong argument value for flag')
end


end

