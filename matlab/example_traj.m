clear all;
close all;
addpath(genpath('.'))

% Compile mex function
mex cmex/ForwardBackward.c

% Create object tc with specific parameters, other parameters to default
tc = TurboCode(...
    'blkLength', 100000, ...
    'numIter', 50 ...
);

% Choose a value of ebno
ebno = 0;

% Compute 3 trajectories
[ia, ie] = tc.computeTraj(ebno, 3);

% Plot the result (EXIT chart is created if not existing)
tc.plotTraj(ebno, ia, ie);
