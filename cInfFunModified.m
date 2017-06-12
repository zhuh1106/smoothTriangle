function Zt = cInfFunModified(beta, fOption, theta, r, tlo, thi)
% beta: ratation angle clockwise
% fOption: flip option
% theta = pi/6;
% r = 0.5;
% beta = 7*pi/6;
% % beta = 0;
% fOption = 1;
% [xt, yt] = cInfFun( theta, r);
[xt, yt] = cInfFun2( theta, r);

if fOption == 0
%     ss.Z = @(t) xt(t/(4*pi)) + 1i*yt(t/(4*pi));
    xxt = @(t) cos(2*pi - beta)*xt(t) - sin(2*pi - beta)*yt(t);
    yyt = @(t) sin(2*pi - beta)*xt(t) + cos(2*pi - beta)*yt(t);
    Zt = @(t) xxt((t-tlo)/(2*(thi-tlo))) + 1i*yyt((t-tlo)/(2*(thi-tlo)));
else
%     ss.Z = @(t) -xt(1/2 - t/(4*pi)) + 1i*yt(1/2 - t/(4*pi));
    xxt = @(t) cos(beta)*xt(t) - sin(beta)*yt(t);
    yyt = @(t) sin(beta)*xt(t) + cos(beta)*yt(t);
    Zt = @(t) -xxt(1/2 - (t-tlo)/(2*(thi-tlo))) + 1i*yyt(1/2 - (t-tlo)/(2*(thi-tlo)));
end

% [ss, N, np] = quadr_pan(ss, N, qtype, qntype);
% 
% % figure()
% plot(real(ss.x),imag(ss.x),'*')
% hold on
% % if fOption == 0
% %     xlim([0 0.6])
% % else
% %     xlim([-0.6 0])
% % end
% % ylim([0 0.6])
% xlim([-1 1])
% ylim([-1 1])
