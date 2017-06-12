function [ss,N] = polygon_test(x,r)

% x = [0,0;2,0;1,0.4;0.5,1];
% x = [0,0;2,0;1,1;0.5,1];
[polySize, ~] = size(x);
% r = 0.1;

qtype = 'p';
qntype = 'C';
bdType = 'polygon';
idNum = 'middle';
N = 48*polySize;

f = @(s) 0;
for i = 0:polySize - 1
    x1 = x(mod(i,polySize)+1,:);
    x2 = x(mod(i+1,polySize)+1,:);
    x3 = x(mod(i+2,polySize)+1,:);
    x4 = x(mod(i+3,polySize)+1,:);
    t = linspace(i*(2*pi)/polySize,(i+1)*(2*pi)/polySize,4);
    Zt = cInfBoundary2([x1; x2; x3; x4], r, t, bdType, idNum);
%     ss.Z = @(s) (s<=pi).*Zt1(s) + (s>pi).*Zt2(s);
    f = @(s) f(s) + (s>t(1)).*(s<=t(end)).*Zt(s);
end

% % x1 = [-1,1];
% x2 = [0,0];
% x3 = [1,0.5];
% x4 = [2,-4];
% t = [0,2*pi/3,pi];
% % N = 64;
% 
% % r = 0.3;
% 
% % Zt1 = cInfBoundary2([x1; x2; x3; x4], r, [t(1),t(2),t(3),t(4)], bdType, idNum);
% Zt1 = cInfBoundary2([x2; x3; x4], r, t, bdType, idNum);
% 
% 
% %%
% 
% 
% 
% 
% x1 = [0,0];
% x2 = [1,0.5];
% x3 = [2,-4];
% % x4 = [3,-2];
% 
% t = [pi,pi+pi/3,2*pi];
% N = 64;
% 
% r = 0.3;
% idNum = 'last';
% 
% % Zt2 = cInfBoundary2([x1; x2; x3; x4], r, [t(1),t(2),t(3),t(4)], bdType, idNum);
% Zt2 = cInfBoundary2([x1; x2; x3], r, t, bdType, idNum);
% 
%% assembly
ss.Z = @(s) f(s);
% ss.Z = Z2;

%% quadrature pts get  
[ss, N, np] = quadr_pan(ss, N, qtype, qntype);

%% figure of cInf approx to triangle
figure()
% plot(real(ss.x),imag(ss.x),'.')
plot(real(ss.x),imag(ss.x))
daspect([1 1 1])
% hold on
% figure()
% plot(real(ss.x(1:end/12)),imag(ss.x(1:end/12)),'.')
% daspect([1 1 1])
% hold on
% plot(real(ss.x(end/12:end/6)),imag(ss.x(end/12:end/6)),'*')
% plot(real(ss.x(end/6:end/4)),imag(ss.x(end/6:end/4)),'o')
% figure()
% plot(real(ss.x(end/4:end/3)),imag(ss.x(end/4:end/3)),'.')
% daspect([1 1 1])
% hold on
% plot(real(ss.x(end/3:5*end/12)),imag(ss.x(end/3:5*end/12)),'*')
% plot(real(ss.x(5*end/12:end/2)),imag(ss.x(5*end/12:end/2)),'o')
% 
% 
% 
