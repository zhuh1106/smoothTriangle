clear all
close all
x1 = [0,0];
x2 = [3,0.5];
% x3 = [1,1];
x3 = [2,-2];

N = 256;
qtype = 'p';
qntype = 'C';

t1 = 0; t2 = 2*pi/3; t3 = 4*pi/3; t4 = 2*pi;
r = 0.04;    % radius



Zt1 = cInfBoundary(x1, x2, x3, r, [0,2*pi/9,4*pi/9,2*pi/3]);
Zt2 = cInfBoundary(x2, x3, x1, r, [2*pi/3,8*pi/9,10*pi/9,4*pi/3]);
Zt3 = cInfBoundary(x3, x1, x2, r, [4*pi/3,14*pi/9,16*pi/9,2*pi]);


%% assembly
ss.Z = @(t) (t<=2*pi/3) .* Zt1(t) + ((t>2*pi/3).* (t<=4*pi/3)) .* Zt2(t)...
        + (t>4*pi/3) .* Zt3(t);

% ss.Z = @(t) (t<=4*pi/3) .* Zt1(t) + (t>4*pi/3) .* (1 + 1i*(t/(2*pi)+0.2));

%% quadrature pts get  
% ss.Z = Zt1;
[ss, N, np] = quadr_pan(ss, N, qtype, qntype);

%% figure of cInf approx to triangle
figure()
plot(real(ss.x),imag(ss.x))
daspect([1 1 1])
hold on
