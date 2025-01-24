% The random timestep Euler method and its continuous dynamics
% Jonas Latz, University of Manchester, 2024 -
%
% Visualise stochastic Euler dynamics for a logistic and a linear ODE

u0 = 0.25;
u1 = @(t) exp(t)./(3+exp(t));
u2 = @(t) 0.25*exp(-t);
T_end = 6;
t_range = 0:0.01:T_end;
u_plot(1,:) = u1(t_range);
u_plot(2,:) = u2(t_range);
N_dim = 2;
close all
% Sampling waiting times
h = 0.8;

for i =1:2
    if i ==1 
        f = @(t,x) f1(t,x);
    else
        f = @(t,x) f2(t,x);
    end
    lambda = 1/h;
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
    subplot(2,2,i)
     plot(t_range,u_plot(i,:),'-',"Color",[0.75 0.75 0.75],'LineWidth',2)
     xlabel('$t$','Interpreter','latex')
       
hold on
 plot(T, U(1,:),'k-','LineWidth',1.5,'MarkerSize',4)
    stairs(T, U(1,:),'o-','LineWidth',1.5,'MarkerSize',4,...
        'MarkerEdgeColor','red',...
        'MarkerFaceColor','red',...
        'Color','red')
   
hold off

end

figure(1)
legend('$u(t)$','$V(t)$','$\overline{V}(t)$','Location','northeast','Interpreter','latex')


function y = f1(t,x)
    y = (1-x).*x;
end

function y = f2(t,x)
    y = -x;
end
