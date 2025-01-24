% The random timestep Euler method and its continuous dynamics
% Jonas Latz, University of Manchester, 2024 -
%
% For a linear ODE u' = au, u(0) = u0, study the convergence of the 
% stochastic Euler dynamics to the deterministic Euler dynamics through
% a Monte Carlo simulation.

clear

%Initialisation
u0 = 1;
a = 1;
u = @(t) exp(-a*t)*u0;

T_end = 10; 
t_range = 0:0.1:T_end;

hs = [0.0625,0.125, 0.25, 0.5, 1, 2]/a;
Js = 100000 * ones(6);

Ts = [0:0.5:10];
N_Ts = size(Ts,2);

U_all{1} = zeros(Js(1),N_Ts);
U_all{2} = zeros(Js(2),N_Ts);
U_all{3} = zeros(Js(3),N_Ts);
U_all{4} = zeros(Js(4),N_Ts);
U_all{5} = zeros(Js(5),N_Ts);
U_all{6} = zeros(Js(6),N_Ts);


for i = 1:6 %Loop stepsize parameters
   
    h = hs(i);

    % eval deterministic Euler dynamics
c = u0;
if abs(1-4*a*h) > eps
w =@(t) (c* exp(-(((1+sqrt(1-4*a*h))*t)/(2*h))).*(-1+2* a* h+sqrt(1-4*a*h)+exp((sqrt(1-4*a*h)*t)/h)*(1-2*a*h+sqrt(1-4* a*h))))/(2*sqrt(1-4*a*h));
w_bar = @(t) (c*exp(-(((1+sqrt(1-4*a*h))*t)/(2*h))).*(-1+sqrt(1-4*a*h)+exp((sqrt(1-4*a*h)*t)/h)*(1+sqrt(1-4*a*h))))/(2*sqrt(1-4*a*h));
else

w = @(t) c* exp(-2*t).*(1+t);
w_bar = @(t) c* exp(-2*t).*(1+2*t);
end

u_plot{i} = w(t_range);

    J = Js(i);
disp(['h = ',num2str(h)])
U_diff = zeros(J,N_Ts);
U_all{i}(:,1) = u0*ones(J,1);

for j = 1:J %J Monte Carlo samples

    if mod(j,50000) == 0

        disp(['Iteration j = ',num2str(j)])
    end
    for nt = 2:N_Ts %sample waiting times
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
        for k = 1:N_T %sample trajectory of stochastic Euler dynamics
            U(k+1) = U(k) - a*DT(k)*U(k);
        end
        U_all{i}(j,nt) = U(end);
        U_diff(j,nt) = (U(end) - U(end-1))^2;
    end

end
% Compute mean/standard deviation at final time
plot_means = mean(U_all{i});
plot_std = std(U_all{i});

plot_diff_means = mean(U_diff);


t_Euler = 0:h:T_end;
N_T_det = length(t_Euler);

U_Eul = zeros(N_T_det,1);
U_Eul(1) = u0;

%Plot
figure(1)
subplot(2,3,i)

plot(t_range, u_plot{i},'r-','LineWidth',1.5) 

hold on 
plot(Ts, plot_means,'ko','LineWidth',1.5,'MarkerSize',3,...
        'MarkerEdgeColor','black',...
        'MarkerFaceColor','black')
%semilogy(t_Euler, U_Eul.^2)

%semilogy(t_range, exp(-t_range/(2*h))*(u0^2),'g-','LineWidth',1.5)
plot(Ts, plot_means+plot_std,'k--','LineWidth',0.75)
plot(Ts, plot_means-plot_std,'k--','LineWidth',0.75)
%semilogy(Ts, plot_means-plot_std,'k:','LineWidth',1)
%semilogy(Ts, plot_diff_means)
hold off
xlabel('$t$','Interpreter','latex')
title(['$h=',num2str(a*h),'$'],'Interpreter','latex')

end
subplot(2,3,1)
legend(  '$w(t)$', '$\widehat{\mathrm{E}}[V(t)]$', '$(\widehat{\mathrm{E}} \pm \widehat{\mathrm{SD}})[V(t)]$','Interpreter','latex','Location','eastoutside')




function y = f(t,x)
    y = - x;
end

function y = f_pot(t,x)
    y = 0.5*x^2;
end

