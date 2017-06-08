function [ss,N] = triangle_test(x1, x2, x3, t, N, r)

qtype = 'p';
qntype = 'C';

Zt1 = cInfBoundary(x1, x2, x3, r, [t(1),t(2),t(3),t(4)]);
Zt2 = cInfBoundary(x2, x3, x1, r, [t(4),t(5),t(6),t(7)]);
Zt3 = cInfBoundary(x3, x1, x2, r, [t(7),t(8),t(9),t(10)]);


%% assembly
ss.Z = @(t) (t<=2*pi/3) .* Zt1(t) + ((t>2*pi/3).* (t<=4*pi/3)) .* Zt2(t)...
        + (t>4*pi/3) .* Zt3(t);

%% quadrature pts get  
[ss, N, np] = quadr_pan(ss, N, qtype, qntype);

%% figure of cInf approx to triangle
figure()
plot(real(ss.x),imag(ss.x))
daspect([1 1 1])
hold on


% x1 = [0,0];
% x2 = [1,1];
% % x3 = [1,1];
% x3 = [1,0];

% N = 128*4;
% t1 = 0; t2 = 2*pi/3; t3 = 4*pi/3; t4 = 2*pi;
