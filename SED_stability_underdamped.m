% The random timestep Euler method and its continuous dynamics
% Jonas Latz, University of Manchester, 2024 -
%
% For an underdamped oscillator, we study the stability of the stochastic
% Euler dynamics

clear
u0 = 1;
a = 1;
T_end = 600; 
t_range = 0:0.1:T_end;


% Sampling waiting times
hs = [0.2, 0.6, 2/3, 0.7]/a;
Js = 100000 * ones(6);

Ts = [0:50:T_end]; 
N_Ts = size(Ts,2);

U_all{1} = zeros(2,Js(1),N_Ts);
U_all{2} = zeros(2,Js(2),N_Ts);
U_all{3} = zeros(2,Js(3),N_Ts);
U_all{4} = zeros(2,Js(4),N_Ts);
U_all{5} = zeros(2,Js(5),N_Ts);
U_all{6} = zeros(2,Js(6),N_Ts);
for i = 1:4
   
    h = hs(i);
    J = Js(i);
disp(['h = ',num2str(h)])
U_diff = zeros(J,N_Ts);
U_all{i}(1,:,1) = ones(J,1);


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
        
        U = zeros(2,N_T+1);
        U(1) = u0;
        for k = 1:N_T
            U(:,k+1) = U(:,k) + DT(k)*f(0,U(:,k));
        end
        U_all{i}(:,j,nt) = U(:,end);
        
    end

end

%for i = 1:4


plot_means = reshape(mean(U_all{i}(1,:,:).^2+U_all{i}(2,:,:).^2),N_Ts,1);
plot_std = reshape(std(U_all{i}(1,:,:).^2+U_all{i}(2,:,:).^2),N_Ts,1);



figure(1)
subplot(2,4,i)

semilogy(Ts, plot_means,'ko-','LineWidth',1.5,'MarkerSize',3,...
        'MarkerEdgeColor','black',...
        'MarkerFaceColor','black')
hold on

semilogy(Ts, plot_means+plot_std,'k--','LineWidth',0.75)

hold off
xlabel('$t$','Interpreter','latex')
 if i ~= 3
        title(['$h = ',num2str(hs(i)),'$'],'Interpreter','latex')
    else
        title('$h = 2/3$','Interpreter','latex')
    end
end
subplot(2,4,1)
legend('$\widehat{\mathrm{E}}[\|V(t)\|^2]$','$(\widehat{\mathrm{E}}+ \widehat{\mathrm{SD}})[\|V(t)\|^2]$','Interpreter','latex','Location','southwest')


function y = f(t,x)
    y = [x(2); -x(1)-x(2)];
end
