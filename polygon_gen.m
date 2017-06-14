close all
clear all
r = 0.2;

x = [0,0;2,0;1,0.4;0.5,1];
[ss,N] = polygon_test(x,r);

x = [0,3;1,1;3,2;2,-1;3,-3;0,-2;-3,-3;-2,-1;-3,2;-1,1];
[ss,N] = polygon_test(x,r);
ss;
