clear all;
close all;
addpath(genpath('.'))

% Compile mex function
mex cmex/ForwardBackward.c

% Create object tc with specific parameters, other parameters to default
tc = TurboCode(...
    'blkLength', 100000, ...
    'algorithm', 'MAP' ...
);

% Choose a value of ebno
ebno = -1;

% Generate an EXIT chart (saved to file), 10 simulations
[ia, ie] = tc.computeExit(ebno, 10);

% Plot EXIT chart
tc.plotExit(ebno, ia, ie);