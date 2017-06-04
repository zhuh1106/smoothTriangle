function Zt = cInfBoundary(x1, x2, x3, r, tvec, beta1, beta2)
% provide a curved line for x1--x2, x3 is also needed for triangle, another
% point is going to be needed 

t1 = tvec(1); t2 = tvec(2); t3 = tvec(3); t4 = tvec(4);


%% starting piece
v = x2 - x1;
u = x3 - x1;
alpha1 = 1/2*acos((u*v')/(norm(u)*norm(v))); % angle formed at x1
l1 = r/tan(alpha1);   % length from x1 to point D on x1--x2 inscribed by circle
xd1 = l1/norm(v)*x2 + (norm(v)-l1)/norm(v)*x1; % coordinates of this point D
xe1 = l1/norm(u)*x3 + (norm(u)-l1)/norm(u)*x1; % similarly, point E on x1--x3
xf1 = 1/2*(xd1 + xe1);     % point going to be used to determine center of circle
xo1 = x1 + [(xf1(1)-x1(1))/cos(alpha1)^2,(xf1(2)-x1(2))/cos(alpha1)^2];
theta1 = pi/2 - alpha1;   % angle formed by cInf function
fOption1 = 0;
% beta1 = pi;

Zt1 = cInfFunModified(beta1, fOption1, theta1, r, t1, t2);


%% end piece   
v = x1 - x2;
u = x3 - x2;
alpha2 = 1/2*acos((u*v')/(norm(u)*norm(v))); % angle formed at x1
l2 = r/tan(alpha2);   % length from x1 to point D on x1--x2 inscribed by circle
xd2 = l2/norm(v)*x1 + (norm(v)-l2)/norm(v)*x2; % coordinates of this point D
xe2 = l2/norm(u)*x3 + (norm(u)-l2)/norm(u)*x2; % similarly, point E on x1--x3
xf2 = 1/2*(xd2 + xe2);     % point going to be used to determine center of circle
xo2 = x2 + [(xf2(1)-x2(1))/cos(alpha2)^2,(xf2(2)-x2(2))/cos(alpha2)^2];
theta2 = pi/2 - alpha2;   % angle formed by cInf function
fOption2 = 1;
% beta2 = pi;

Zt3 = cInfFunModified(beta2, fOption2, theta2, r, t3, 2*pi);

%% middle piece
Zt2 = @(t) xd1(1) + (t-t2)/(t3-t2)*(xd2(1)-xd1(1)) ...
        + 1i*( xd1(2) + (t-t2)/(t3-t2)*(xd2(2)-xd1(2)) );




%% assembly
Zt = @(t) (t<=t2) .* (xo1(1) + 1i*xo1(2) + Zt1(t)) + ((t>t2).* (t<=t3)) .* (Zt2(t))...
        + (t>t3) .* (xo2(1) + 1i*xo2(2) + Zt3(t));