% The random timestep Euler method and its continuous dynamics
% Jonas Latz, University of Manchester, 2024 -
%
% For a harmonic osciallator, we simulate the deterministic Euler dynamics.

hs = [0.2, 0.6, 2/3, 0.7];
A = [0 1; -1 -1];
u = @(t) [(1/3)*exp(-t/2).*(sqrt(3)*sin((sqrt(3)*t)/2) + 3* cos((sqrt(3)* t)/2));-(2*exp(-t/2).*sin((sqrt(3)*t)/2))/sqrt(3)];
t = 0:0.01:40;
u_plot = u(t);

B = @(h) [zeros(2), A; (1/h)*eye(2), -(1/h)*eye(2)];


Nt = length(t);


u0 = [1;0;1;0];
for i = 1:4
    U_all{i} = zeros(4,Nt);
    U_all{i}(:,i) = [1;0;1;0];
    for k = 1:Nt
        U_all{i}(:,k) = expm(t(k)*B(hs(i)))*u0;
    end
    figure(1)
    subplot(2,4,i)
    plot(t,u_plot(1,:),'-',"Color",[0.75 0.75 0.75],'LineWidth',2) 
    hold on
    plot(t, U_all{i}(1,:),'-',"Color",[0 0 0],'LineWidth',1.5)
    plot(t, U_all{i}(3,:),'-',"Color","red",'LineWidth',1.5)
    hold off
    if i == 1
        ylabel('position $u_1$','Interpreter','latex')
    end
    if i ~= 3
        title(['$h = ',num2str(hs(i)),'$'],'Interpreter','latex')
    else
        title('$h = 2/3$','Interpreter','latex')
    end
    subplot(2,4,4+i)
    plot(t,u_plot(2,:),'-',"Color",[0.75 0.75 0.75],'LineWidth',2)
    hold on
    plot(t, U_all{i}(2,:),'-',"Color",[0 0 0],'LineWidth',1.5)
    plot(t, U_all{i}(4,:),'-',"Color","red",'LineWidth',1.5)
    hold off
    if i == 1
        ylabel("velocity $u_2$",'Interpreter','latex')
    end
    xlabel("$t$",'Interpreter','latex')
end

figure(1)
subplot(2,4,5)
legend('truth','$w(t)$','$\overline{w}(t)$','Location','southeast','Interpreter','latex')

