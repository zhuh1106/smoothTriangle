function [su, N] = pipe_test(xu, r)

[uWallSize, ~] = size(xu); uWallSize = uWallSize - 1;
% [lWallSize, ~] = size(xl); lWallSize = lWallSize - 1;
% r = 0.1;

qtype = 'p';
qntype = 'C';
bdType = 'polygon';
% idNum = 'middle';
N = 48*uWallSize;

%% predefine function, and then add each piece
f = @(s) 0;

%% first piece of upper wall
idNum = 'first';
x1 = xu(1,:); x2 = xu(2,:); x3 = xu(3,:);
% t = linspace(0,2*pi/uWallSize,3);
t = [0,4*pi/(3*uWallSize),2*pi/uWallSize];
Zt = cInfBoundary2([x1; x2; x3], r, t, bdType, idNum);
f = @(s) f(s) + (s>=t(1)).*(s<t(end)).*Zt(s);

%% middle pieces of upper wall
idNum = 'middle';
for i = 2:uWallSize - 1 % 2nd to 2nd to last
    x1 = xu(i-1,:); x2 = xu(i,:); x3 = xu(i+1,:); x4 = xu(i+2,:);
    t = linspace((i-1)*(2*pi)/uWallSize,i*(2*pi)/uWallSize,4);
    Zt = cInfBoundary2([x1; x2; x3; x4], r, t, bdType, idNum);
%     ss.Z = @(s) (s<=pi).*Zt1(s) + (s>pi).*Zt2(s);
    f = @(s) f(s) + (s>=t(1)).*(s<t(end)).*Zt(s);
end

%% last piece of upper wall
idNum = 'last';
x1 = xu(end-2,:); x2 = xu(end-1,:); x3 = xu(end,:);
% t = linspace((uWallSize-1)*(2*pi)/uWallSize,2*pi,3);
t = [(uWallSize-1)*(2*pi)/uWallSize,2*pi-4*pi/(3*uWallSize),2*pi];
Zt = cInfBoundary2([x1; x2; x3], r, t, bdType, idNum);
f = @(s) f(s) + (s>=t(1)).*(s<=t(end)).*Zt(s);

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
su.Z = @(s) f(s);
% ss.Z = Z2;

%% quadrature pts get  
[su, N, np] = quadr_pan(su, N, qtype, qntype);

%% figure of cInf approx to triangle
figure()
% plot(real(ss.x),imag(ss.x),'.')
plot(real(su.x),imag(su.x),'.')
daspect([1 1 1])

figure()
plot(linspace(0,2*pi,N),su.cur)


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
