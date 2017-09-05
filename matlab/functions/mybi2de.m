function [ D ] = mybi2de( varargin )
%MYBI2DE Convert binary input (in double) to decimal
%   D = mybi2de(B)  where B is the binary input matrix (right-msb is 
%   assumed)
%
%   D = mybi2de(B, flag)    where B is the binary input matrix and flag is
%   equal to right-msb or left-msb

% Arg check
switch length(varargin)
    case 1
        B       = varargin{1};
        flag    = 'right-msb';
    case 2
        B       = varargin{1};
        flag    = varargin{2};
    case 0
        error('Not enough input argument')
    otherwise
        error('Too many input arguments')
end

% Compute base depending on the flag
if strcmp(flag, 'right-msb')
    base        = 2.^(0:(size(B, 2)-1))';
elseif strcmp(flag, 'left-msb')
    base        = 2.^((size(B, 2)-1):-1:0)';
else
    error('Wrong argument value for flag')
end

% Matrix product to get the decimal output
D = B*base;


end

