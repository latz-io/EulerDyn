% The random timestep Euler method and its continuous dynamics
% Jonas Latz, University of Manchester, 2024 -
%
% Plot paths of the deterministic Euler dynamics for the linear ODE
% u' = au, u(0)=u0

c = 1;
a = 1;
h_s = [0.0625, 0.125, 0.25, 0.5, 1, 2];
t_end = 40;

for i = 1:6
h = h_s(i);

u = @(t) exp(-a*t)*c;
if abs(1-4*a*h) > eps
w =@(t) (c* exp(-(((1+sqrt(1-4*a*h))*t)/(2*h))).*(-1+2* a* h+sqrt(1-4*a*h)+exp((sqrt(1-4*a*h)*t)/h)*(1-2*a*h+sqrt(1-4* a*h))))/(2*sqrt(1-4*a*h));
w_bar = @(t) (c*exp(-(((1+sqrt(1-4*a*h))*t)/(2*h))).*(-1+sqrt(1-4*a*h)+exp((sqrt(1-4*a*h)*t)/h)*(1+sqrt(1-4*a*h))))/(2*sqrt(1-4*a*h));
else

w = @(t) c* exp(-2*t).*(1+t);
w_bar = @(t) c* exp(-2*t).*(1+2*t);
end


t_plot = 0:0.01:t_end;
u_plot=real(u(t_plot));
w_plot=real(w(t_plot));
w_bar_plot=real(w_bar(t_plot));


figure(1)
subplot(2,3,i)
plot(t_plot,u_plot,'-',"Color",[0.75 0.75 0.75],'LineWidth',2)
hold on
plot(t_plot,w_plot,'-',"Color",[0 0 0],'LineWidth',1.5)
plot(t_plot,w_bar_plot,'-',"Color","red",'LineWidth',1.5)
hold off
if true
    xlabel('$t$','Interpreter','latex')
end
ylim([-0.75,1])
title(['$h=',num2str(h),'$'],'Interpreter','latex')

end

legend('$u(t)$','$w(t)$','$\overline{w}(t)$','Location','northeast','Interpreter','latex')
