function Uniform_Solver_S(strain_rate, time_final, time_res, space_res, max_k)
% Function produces solutions for uniform growth cases (solves the PDE & numerical modes & generates instability regions)

%% Parameters used in the reaction-diffusion system
a=0;
b=1.1;
d_1 = 0.004;
d_2 = 0.1;
Tf=time_final;

% Schnakenberg reaction kinetics
function zvec = kinetics(w)
    zvec = [a-w(1)+w(1)^2*w(2); b-w(1)^2*w(2)];
end

% Jacobian matrix of rection kinetics
function M = J_kinetics(w)
    M = [-1+2*w(1)*w(2), w(1)^2; -2*w(1)*w(2), -w(1)^2];
end

% Homogenous equilibrium of reaction kinetics
IC_nopb = [a+b; b/(a+b)^2];

%% Growth functions for uniform evolutions
S=strain_rate;
time_for_pre = linspace(0,Tf,time_res*10);

% Precaculate r(t) for lots of values
r_pre = exp(cumtrapz(time_for_pre, S(time_for_pre)));

% Interpolate between precalculated values
function z = r(tt)
    z = interp1(time_for_pre, r_pre,tt);
end

%% Solve the reaction-diffusion PDEs

function w0 = pdeic(x)
    w0 = IC_nopb + 0.05*(rand(1)-0.5)*[1;1];
end

function [c,f,s] = pdefun(x,t,w,dw)
    c = [1;1];
    f = (1/r(t).^2)*[d_1; d_2].* dw;
    s = kinetics(w)- S(t)*w;
end

function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
    pl = [0;0];
    pr = [0;0];
    ql = [1;1];
    qr = [1;1];
end

x = linspace(0,1,space_res);
t = linspace(0,Tf,time_res);

sol = pdepe(0, @pdefun, @pdeic, @pdebc, x, t);
u = sol(:,:,1);
v = sol(:,:,2);

%% Plot of the u component from the reaction-diffusion system

figure('color','white')

[X,T] = meshgrid(x,t);

% Display domain evolution
r_vec = zeros(size(t));
for ii =1:length(t)
    r_vec(ii) = r(t(ii));
end
X = r_vec'.*X;

ax1 = subplot(1,3,1);
pcolor(X,T,u)
colorbar('eastoutside');
clim([0 3]);
xlabel('Position $x$', Interpreter='latex', FontSize=18)
ylabel('Time $t$', Interpreter='latex', FontSize=18)
title('System evolution', Interpreter='latex', FontSize=20)
hold on
colormap(ax1,parula)
plot(r_vec,t,'red',LineStyle='-',LineWidth=2)
plot(zeros(size(t)),t,'red',LineStyle='-',LineWidth=2)
hold off
shading interp

%% Plot coefficients of the spatial modes over time
krange = 1:max_k;
mag_matrix = zeros(length(t),max_k);
mag_maxs = zeros(length(t),1);
[K,T] = meshgrid(krange,t);

% Solve for base state, used to extract correct modes of (u-u_0,v-v_0)
[~,base_specific] = ode45(@(t,y) kinetics(y)- S(t)*y, t, IC_nopb);

for tt = 1:length(t)
    u_coef = extract_gen_fourier_coff(u(tt,:)-base_specific(tt,1), x, max_k);
    v_coef = extract_gen_fourier_coff(v(tt,:)-base_specific(tt,2), x, max_k);
    mags = sqrt(u_coef.^2 + v_coef.^2);
    mag_matrix(tt,:) = mags;
    [~,ind] = max(mags);
    mag_maxs(tt) = krange(ind);
end

ax2 = subplot(1,3,2);
pcolor(K,T,mag_matrix)
xlabel('Wavenumber $k$', Interpreter='latex', FontSize=18)
ylabel('Time $t$', Interpreter='latex', FontSize=18)
title('Numerically extracted modes', Interpreter='latex', FontSize=20)
shading flat

%% Analytic instability regions (instability regions implied from Corollary 3.1)

% Solve for the time-varying spatially homogenous base state (u_0(t),v_0(t))
[t_base,y_base] = ode45(@(t,y) kinetics(y)- S(t)*y, t, IC_nopb);

% Use Corollary 3.1 to predict regions of instability
stability_matrix = uniform_condition_calc(t_base, y_base, S(t_base), r(t_base));

% Plot
ax3 = subplot(1,3,3);
[K,T] = meshgrid(krange,t_base);
hold on
pcolor(K,T,1-double(stability_matrix))
plot(mag_maxs(10:end),t(10:end),'red',LineStyle='-',LineWidth=2);
hold off
xlabel('Wavenumber $k$', Interpreter='latex', FontSize=18)
ylabel('Time $t$', Interpreter='latex', FontSize=18)
axis([1,max_k,0,Tf]);
title('Analytic unstable modes', Interpreter='latex', FontSize=20)
colormap(ax3, gray(2))
shading flat

%% Helper functions (Appying Corollary 3.1 function in uniform_condition_calc)

function M = uniform_condition_calc(t, y_base, Sb, r)
% Implementation of Corollary 3.1
% Returns a matrix representing the modes we expect to be unstable at each
% time step, based on linearisation about the steady-state solution y_base

% Init matricies
LHS_matrix = zeros(length(t), max_k);
U=y_base(:,1);
V=y_base(:,2);

fu = zeros(length(t), 1);
gu = zeros(length(t), 1);
fv = zeros(length(t), 1);
gv = zeros(length(t), 1);
for jj=1:length(t)
    Jac = J_kinetics([U(jj), V(jj)]);
    fu(jj) = Jac(1,1);
    gu(jj) = Jac(1,2);
    fv(jj) = Jac(2,1);
    gv(jj) = Jac(2,2);
end

% Calculate LHS of inequality in Corollary 3.1
for k=1:max_k
    %% 
    LHS_matrix(:,k) = deter(U,V) + Sb.*(Sb-trac(U,V)-(k^2*pi^2./(r.^2)) * (d_1+d_2))...
        - (k^2*pi^2./(r.^2)) .* (d_1 * gv + d_2 * fu)...
        + (k^2*pi^2./(r.^2)).^2 * d_1 * d_2;
end

% Calculate RHS of inequality in Corollary 3.1
d1 = n_Deriv((fu-Sb)./fv, t);
d2 = n_Deriv(1./(r.^2.*fv), t);
d3 = n_Deriv((gv-Sb)./gu, t);
d4 = n_Deriv(1./(r.^2.*gu), t);


RHS_matrixA = (fv.*d1)*ones(1,max_k) + d_1*(fv.*d2)*((1:max_k).^2*pi^2);
RHS_matrixB = (gu.*d3)*ones(1,max_k) + d_2*(gu.*d4)*((1:max_k).^2*pi^2);

RHS_matrix = max(RHS_matrixA,RHS_matrixB);

% Return True/False matrix for where condition holds
M = (RHS_matrix - LHS_matrix) >0 ;

end

% Helper functions
function v = n_Deriv(x,t) % Numerical calculation of derivative
    xshift = [0;x(1:end-1)];
    tshift = [-1;t(1:end-1)];
    v = (x-xshift)./(t-tshift);
    v(1) = v(2);
end

function d=deter(u,v) % Det of jacobian
    Jac = J_kinetics([u v]);
    d= Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1);
end

function d=trac(u,v) % Tr of jacobian
    Jac = J_kinetics([u v]);
    d= Jac(1,1)+Jac(2,2);
end

end
