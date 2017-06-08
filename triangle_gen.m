close all
clear all

x = randi([-10 10],1,4);
x1 = [0,0];
x2 = [x(1),x(2)];
x3 = [x(3),x(4)];

r = 0.1;
% vec = [x1;x2;x3];
% area = 1/2*norm(cross(vec(:,1),vec(:,2)));
% R = 2*area/(norm(x2-x1)+norm(x3-x1)+norm(x3-x2));
% r = R - 0.05;
result = radius_test(x1,x2,x3,r);

if result == 1
    N = 16*36*2;
    t = linspace(0,2*pi,10);
    [ss,N] = triangle_test(x1, x2, x3, t, N, r);
    figure()
    plot(linspace(0,2*pi,N),ss.cur)
    title('curvature')
    hold on
    plot(t,zeros(size(t)),'*')
    hold off
    
    fprintf('smooth triangle generated\n')
else
    fprintf('radius is too large\n')
end