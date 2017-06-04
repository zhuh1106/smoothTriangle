clear all
close all
x1 = [0,0];
x2 = [1,0];
x3 = [1,1];

N = 128;
qtype = 'p';
qntype = 'C';

t1 = 0; t2 = 2*pi/3; t3 = 4*pi/3;
r = 0.02;    % radius

% Zt1 = cInfBoundary(x1, x2, x3, r, [t1,t2,t3,2*pi], pi, pi);
% Zt2 = cInfBoundary(x2, x3, x1, r, [t1,t2,t3,2*pi], pi/2, pi/2);
% Zt3 = cInfBoundary(x3, x1, x2, r, [t1,t2,t3,2*pi], 7*pi/4, 7*pi/4);
% no issue when run them one by one

Zt1 = cInfBoundary(x1, x2, x3, r, [0,2*pi/9,4*pi/9,2*pi/3], pi, pi);
Zt2 = cInfBoundary(x2, x3, x1, r, [2*pi/3,8*pi/9,10*pi/9,4*pi/3], pi/2, pi/2);
Zt3 = cInfBoundary(x3, x1, x2, r, [4*pi/3,14*pi/9,16*pi/9,2*pi], 7*pi/4, 7*pi/4);

% Zt1 = cInfBoundary(x1, x2, x3, r, [0,2*pi/9,2*pi-2*pi/9,2*pi]);
% Zt2 = cInfBoundary(x2, x3, x1, r, [2*pi/3,8*pi/9,10*pi/9,4*pi/3]);
% Zt3 = cInfBoundary(x3, x1, x2, r, [4*pi/3,14*pi/9,16*pi/9,2*pi]);

% %% starting piece
% v = x2 - x1;
% u = x3 - x1;
% alpha1 = 1/2*acos((u*v')/(norm(u)*norm(v))); % angle formed at x1
% l1 = r/tan(alpha1);   % length from x1 to point D on x1--x2 inscribed by circle
% xd1 = l1/norm(v)*x2 + (norm(v)-l1)/norm(v)*x1; % coordinates of this point D
% xe1 = l1/norm(u)*x3 + (norm(u)-l1)/norm(u)*x1; % similarly, point E on x1--x3
% xf1 = 1/2*(xd1 + xe1);     % point going to be used to determine center of circle
% xo1 = x1 + [(xf1(1)-x1(1))/cos(alpha1)^2,(xf1(2)-x1(2))/cos(alpha1)^2];
% theta1 = pi/2 - alpha1;   % angle formed by cInf function
% fOption = 0;
% beta = pi;
% 
% Zt1 = cInfFunModified(beta, fOption, theta1, r, t1, t2);
% 
% 
% %% end piece   
% v = x1 - x2;
% u = x3 - x2;
% alpha2 = 1/2*acos((u*v')/(norm(u)*norm(v))); % angle formed at x1
% l2 = r/tan(alpha2);   % length from x1 to point D on x1--x2 inscribed by circle
% xd2 = l2/norm(v)*x1 + (norm(v)-l2)/norm(v)*x2; % coordinates of this point D
% xe2 = l2/norm(u)*x3 + (norm(u)-l2)/norm(u)*x2; % similarly, point E on x1--x3
% xf2 = 1/2*(xd2 + xe2);     % point going to be used to determine center of circle
% xo2 = x2 + [(xf2(1)-x2(1))/cos(alpha2)^2,(xf2(2)-x2(2))/cos(alpha2)^2];
% theta2 = pi/2 - alpha2;   % angle formed by cInf function
% fOption = 1;
% beta = pi;
% 
% Zt3 = cInfFunModified(beta, fOption, theta2, r, t3, 2*pi);
% 
% %% middle piece
% Zt2 = @(t) xd1(1) + (t-t2)/(t3-t2)*(xd2(1)-xd1(1)) ...
%         + 1i*( xd1(2) + (t-t2)/(t3-t2)*(xd2(2)-xd1(2)) );
% 
% 
% 
% 
%% assembly
ss.Z = @(t) (t<=2*pi/3) .* Zt1(t) + ((t>2*pi/3).* (t<=4*pi/3)) .* Zt2(t)...
        + (t>4*pi/3) .* Zt3(t);

% ss.Z = @(t) (t<=2*pi/3) .* Zt1(t) + (t>2*pi/3) .* (1 + 1i*(t/(2*pi)+0.2));

%% quadrature pts get  
% ss.Z = Zt1;
[ss, N, np] = quadr_pan(ss, N, qtype, qntype);

%% figure of cInf approx to triangle
figure()
plot(real(ss.x),imag(ss.x),'*')
daspect([1 1 1])
hold on
