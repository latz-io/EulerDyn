% The random timestep Euler method and its continuous dynamics
% Jonas Latz, University of Manchester, 2024 -
%
% Visualise convergence of the stochastic Euler dynamics as h is reduced
% for an underdamped oscillator

u0 = [1;0];
u = @(t) [(1/3)*exp(-t/2).*(sqrt(3)*sin((sqrt(3)*t)/2) + 3* cos((sqrt(3)* t)/2));-(2*exp(-t/2).*sin((sqrt(3)*t)/2))/sqrt(3)];
T_end = 10;
t_range = 0:0.01:T_end;
u_plot = u(t_range);
N_dim = 2;
close all
% Sampling waiting times
h = [0.25, 0.125, 0.0625];

figure(1)
subplot(2,2,1)
hold on
plot(t_range,u_plot(1,:),'-',"Color",[0.75 0.75 0.75],'LineWidth',2)
hold off
ylabel('position $u_1$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
subplot(2,2,2) 
hold on
plot(t_range,u_plot(2,:),'-',"Color",[0.75 0.75 0.75],'LineWidth',2)
hold off
ylabel('velocity $u_2$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
for i=1:length(h)
    lambda = 1/h(i);
    T = 0;
    DT = [];
    while T(end) < T_end
        DT_new = exprnd(1/lambda);
        T = [T T(end)+DT_new];
        DT = [DT DT_new];
    end
    T(end) = T_end;
    DT(end) = T_end - T(size(T,2)-1);
    N_T = size(DT,2);
    
    
    U = zeros(N_dim,N_T+1);
    U(:,1) = u0;
    for k = 1:N_T
        U(:,k+1) = U(:,k) + DT(k)*f(T(k),U(:,k));
    end
    
    figure(1)
    subplot(2,2,1)
    hold on
    plot(T, U(1,:),'-','LineWidth',1.)
    
    subplot(2,2,2)
    hold on
    plot(T, U(2,:),'-','LineWidth',1.)

    
end
figure(1)
legend('truth','$h=0.25$','$h=0.125$','$h=0.0625$','Location','southeast','Interpreter','latex')



function y = f(t,x)
    y = [x(2); -x(1)-x(2)];
end
