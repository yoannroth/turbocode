clear all;
close all;
addpath(genpath('.'))

% Compile mex function
mex cmex/ForwardBackward.c

% Create object tc with specific parameters, other parameters to default
tc = TurboCode(...
    'blkLength', 1000, ...
    'algorithm', 'MAP', ...
    'numIter', 10, ...
    'poly', [13 15], ...
    'k', 4 ...
);

% Choose a range of ebno values
ebno = -1:0.3:2;

% Compute 3 trajectories
ber = tc.computeBer(ebno);

% Plot the result (EXIT chart is created if not existing)
tc.plotBer(ebno, ber);
