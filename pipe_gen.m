close all
clear all
r = 0.05;

xu = [0,1;1,0;2,0;3,1];
[su, N] = pipe_test(xu, r);

% x = [0,3;1,1;3,2;2,-1;3,-3;0,-2;-3,-3;-2,-1;-3,2;-1,1];
% [ss,N] = polygon_test(x,r);
% ss;