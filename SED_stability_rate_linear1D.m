% The random timestep Euler method and its continuous dynamics
% Jonas Latz, University of Manchester, 2024 -
%
% For a linear ODE u' = au, u(0) = u0, study the longtime behaviour of the 
% stochastic Euler dynamics with different stepsize parameters
clear
u0 = 1;
a = 1;
u = @(t) exp(-a*t)*u0;
logu = @(t) -t + log(u0);

T_end = 60; 
t_range = 0:0.1:T_end;
u_plot = u(t_range);

hs = [0.125, 0.25, 0.5, 1, 2]/a;
Js = 1000000 * ones(6);

Ts = [0:4:T_end]; 
N_Ts = size(Ts,2);

U_all{1} = zeros(Js(1),N_Ts);
U_all{2} = zeros(Js(2),N_Ts);
U_all{3} = zeros(Js(3),N_Ts);
U_all{4} = zeros(Js(4),N_Ts);
U_all{5} = zeros(Js(5),N_Ts);

for i = 1:5
   
    h = hs(i);
    J = Js(i);
disp(['h = ',num2str(h)])
U_diff = zeros(J,N_Ts);
U_all{i}(:,1) = u0*ones(J,1);

for j = 1:J

    if mod(j,50000) == 0

        disp(['Iteration j = ',num2str(j)])
    end
    for nt = 2:N_Ts
    T = 0;
    DT = [];
    ijk = 1;
    DT_new = exprnd(h,ceil(10*T_end/h),1);
        while T(end) < Ts(nt)
            T = [T T(end)+DT_new(ijk)];
            DT = [DT DT_new(ijk)];
            ijk = ijk+1;
        end
        T(end) = Ts(nt);
        DT(end) = Ts(nt) - T(size(T,2)-1);
        N_T = size(DT,2);
        
        U = zeros(N_T+1,1);
        U(1) = u0;
        for k = 1:N_T
            U(k+1) = U(k) - a*DT(k)*U(k);
        end
        U_all{i}(j,nt) = U(end);
        U_diff(j,nt) = (U(end) - U(end-1))^2;
    end

end

plot_means = mean(U_all{i}.^2);
plot_std = std(U_all{i}.^2);

plot_diff_means = mean(U_diff);


t_Euler = 0:h:T_end;
N_T_det = length(t_Euler);

U_Eul = zeros(N_T_det,1);
U_Eul(1) = u0;


figure(1)
subplot(2,3,i)

semilogy(Ts, plot_means,'ko-','LineWidth',1.5,'MarkerSize',3,...
        'MarkerEdgeColor','black',...
        'MarkerFaceColor','black')
hold on

semilogy(t_range, u_plot.^2,'b-','LineWidth',1.5)
semilogy(t_range, exp(-t_range/(2*h))*(u0^2),'g-','LineWidth',1.5)
semilogy(Ts, plot_means+plot_std,'k--','LineWidth',0.75)

hold off
xlabel('$t$','Interpreter','latex')
title(['$h=',num2str(a*h),'$'],'Interpreter','latex')
end
subplot(2,3,5)
legend('$\widehat{\mathrm{E}}[V(t)^2]$',  '$\exp(-2t)$', '$\exp(-t/(2h))$','$(\widehat{\mathrm{E}}+ \widehat{\mathrm{SD}})[V^2(t)]$','Interpreter','latex','Location','eastoutside')



function y = f(t,x)
    y = - x;
end

function y = f_pot(t,x)
    y = 0.5*x^2;
end

