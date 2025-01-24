% The random timestep Euler method and its continuous dynamics
% Jonas Latz, University of Manchester, 2024 -
%
% For an underdamped oscillator, we sample associated stochastic Euler
% dynamics and compare them to the analytical solution of the ODE.

u0 = [1;0];
u = @(t) [(1/3)*exp(-t/2).*(sqrt(3)*sin((sqrt(3)*t)/2) + 3* cos((sqrt(3)* t)/2));-(2*exp(-t/2).*sin((sqrt(3)*t)/2))/sqrt(3)];
T_end = 400;
t_range = 0:0.01:T_end;
u_plot = u(t_range);
N_dim = 2;

h = [0.2, 0.6, 2/3, 0.7];
for i = 1:4
figure(1)
subplot(2,4,i)
plot(t_range,u_plot(1,:),'-',"Color",[0.75 0.75 0.75],'LineWidth',2)
if i == 1
ylabel('position $u_1$','Interpreter','latex')
end

 if i ~= 3
        title(['$h = ',num2str(h(i)),'$'],'Interpreter','latex')
    else
        title('$h = 2/3$','Interpreter','latex')
    end
subplot(2,4,i+4) 

plot(t_range,u_plot(2,:),'-',"Color",[0.75 0.75 0.75],'LineWidth',2)
if i == 1
    ylabel('velocity $u_2$','Interpreter','latex')
end
xlabel('$t$','Interpreter','latex')
end
for i=1:length(h)
    for ijk = 1:5
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
    subplot(2,4,i)
    hold on
    plot(T, U(1,:),'k-','LineWidth',0.5)
    hold off
    subplot(2,4,4+i)
    hold on
    plot(T, U(2,:),'k-','LineWidth',0.5)
    hold off
    end
end
figure(1)
subplot(2,4,5)
legend('truth','$V(t)$','Location','southeast','Interpreter','latex')

function y = f(t,x)
    y = [x(2); -x(1)-x(2)];
end
