% The random timestep Euler method and its continuous dynamics
% Jonas Latz, University of Manchester, 2024 -
%
% Estimate the root mean square truncation error of stochastic Euler
% dynamics and second order stochastic Euler dynamics for u' = au, u(0)=u0

clear
u0 = 1;
a = 1;
u = @(t) exp(-a*t)*u0;
logu = @(t) -t + log(u0);
%log_T_end = 9;
T_end = 1; %2^log_T_end;
t_range = 0:0.001:T_end;
u_plot = u(t_range);


hs = [0.1, 1];
Js = 100000 * ones(6);

Ts = 2.^(-8:0);
N_Ts = size(Ts,2);
U_all{1} = zeros(Js(1),N_Ts);
U_all{2} = zeros(Js(2),N_Ts);
U_all{3} = zeros(Js(3),N_Ts);
U_all{4} = zeros(Js(1),N_Ts);
U_all{5} = zeros(Js(2),N_Ts);
U_all{6} = zeros(Js(3),N_Ts);


% First order
for i = 1:3
  
    J = Js(i);

U_diff = zeros(J,N_Ts);
U_all{i}(:,1) = u0*ones(J,1);

for j = 1:J
    for nt = 1:N_Ts
    if i == 3
        h = Ts(nt);
    else
        h = hs(i);
    end
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
        U_all{i}(j,nt) = U(end)-exp(-Ts(nt));
       
    end

end
end

% second order
for i = 1:3
  
    J = Js(i);
U_diff = zeros(J,N_Ts);
U_all{i+3}(:,1) = u0*ones(J,1);

for j = 1:J
    for nt = 1:N_Ts
    if i == 3
        h = Ts(nt);
    else
        h = hs(i);
    end
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
            
            U(k+1) = U(k) - a*DT(k)*U(k) + a^2 * DT(k)^2*U(k)/2;
       

        end
        U_all{i+3}(j,nt) = U(end)-exp(-Ts(nt));
      
    end

end
end

for i = 1:3

plot_rmste = sqrt(mean(U_all{i}.^2));
plot_std = std(U_all{i}.^2);
plot_rmste_sec = sqrt(mean(U_all{i+3}.^2));
plot_std_sec = std(U_all{i+3}.^2);

plot_diff_means = mean(U_diff);


t_Euler = 0:h:T_end;
N_T_det = length(t_Euler);

U_Eul = zeros(N_T_det,1);
U_Eul(1) = u0;


figure(1)
subplot(2,3,i)

loglog(Ts, plot_rmste,'ko-','LineWidth',1.5,'MarkerSize',3,...
        'MarkerEdgeColor','black',...
        'MarkerFaceColor','black')
hold on

loglog(Ts, plot_rmste_sec,'o-','Color','[0.8500 0.3250 0.0980]','LineWidth',1.5,'MarkerSize',3,...
        'MarkerEdgeColor','[0.8500 0.3250 0.0980]',...
        'MarkerFaceColor','[0.8500 0.3250 0.0980]')
loglog(t_range,t_range.^2,'k-.','LineWidth',1.0)
loglog(t_range,t_range.^3,'-.','Color','[0.8500 0.3250 0.0980]','LineWidth',1.0)

hold off
xlabel('$\varepsilon$','Interpreter','latex')
if i == 3
    title(['$h=\varepsilon$'],'Interpreter','latex')
else
    title(['$h=',num2str(hs(i)),'$'],'Interpreter','latex')
end
end
legend('$\widehat{\mathrm{E}}[\|V(\varepsilon)-u(\varepsilon)\|^2]^{1/2}$','$\widehat{\mathrm{E}}[\|Y(\varepsilon)-u(\varepsilon)\|^2]^{1/2}$',  '$\varepsilon^2$', '$\varepsilon^3$','Interpreter','latex','Location','northwest','NumColumns',2)



function y = f(t,x)
    y = - x;
end

function y = f_pot(t,x)
    y = 0.5*x^2;
end

