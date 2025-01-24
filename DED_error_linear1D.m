% The random timestep Euler method and its continuous dynamics
% Jonas Latz, University of Manchester, 2024 -
%
% For a linear ODE u' = au, u(0) = u0, we study the convergence of the
% deterministic Euler dynamics to the ODE solution as h -> 0.

clear all
close all
c = 1;
a = 1;
ts = [0.01, 0.1, 1];
hs = 0.0001:0.00005:0.5;
Nts = length(ts);
Nh = length(hs);

u = @(t) exp(-a*t)*c;

w =@(t,h) (c* exp(-(((1+sqrt(1-4*a*h))*t)/(2*h))).*(-1+2* a* h+sqrt(1-4*a*h)+exp((sqrt(1-4*a*h)*t)/h)*(1-2*a*h+sqrt(1-4* a*h))))/(2*sqrt(1-4*a*h));
w_bar = @(t,h) (c*exp(-(((1+sqrt(1-4*a*h))*t)/(2*h))).*(-1+sqrt(1-4*a*h)+exp((sqrt(1-4*a*h)*t)/h)*(1+sqrt(1-4*a*h))))/(2*sqrt(1-4*a*h));


for i = 1:3
    for j = 1:Nh
        error_uw(i,j) = abs(u(ts(i))-w(ts(i),hs(j)));
        error_wwbar(i,j) = abs(w(ts(i),hs(j))-w_bar(ts(i),hs(j)));
    end
end




for i = 1:Nts
figure(1)
subplot(2,3,i)
loglog(hs,error_uw(i,:),'-',"Color",[0.75 0.75 0.75],'LineWidth',1.5)
hold on
loglog(hs,error_wwbar(i,:),'--',"Color",'red','LineWidth',1.5)
loglog(hs,hs,'-',"Color",'blue','LineWidth',1.5)
hold off
if true
    xlabel('$h$','Interpreter','latex')
end
ylim([1E-6,1])

title(['$t=',num2str(ts(i)),'$'],'Interpreter','latex')

end
legend('$\|u(t)-w(t)\|$','$\|w(t)-\overline{w}(t)\|$','h','Location','southeast','Interpreter','latex')
