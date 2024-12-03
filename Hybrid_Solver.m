function Hybrid_Solver
% Function produces the plots for the hybrid growths (solve the PDE + numerical modes + linear system)

%% Parameters
% Params for model

a=0.01;
b=1;
d_1 = 0.004;
d_2 = 0.1;
Tf=100;

% Growth params (can vary in time, for simplicity and ease of computation
% just choose to be constants, still get cool behaviour)
rate_int = 0.00001;
rate_L = -0.03;
rate_R = 0.13;

% Solution params
max_mode = 65;
kappa = max_mode+1;
x = linspace(0,1,1001);
t = linspace(0,Tf,1001);

% ICs for linear ODE system
IC1 = 0.001*(rand(2,kappa)-0.5);

%% Solve PDE

function zvec = kenetics(w)
    zvec = [a-w(1) + w(1)^2 * w(2); b - w(1)^2 * w(2)];
end

function u0 = pdeic(xx)
    u0 = [(a+b); b/(a+b)^2];
    for kk=0:max_mode
        u0 = u0 + (IC1(:,kk+1)*sqrt(2)*cos(pi*kk*xx));
    end
end

function [c,f,s] = pdefun(x,t,w,dw)

    c = [1;1];
    f = (1/r(t).^2)*[d_1; d_2].* dw;
    s = kenetics(w) + S_int(t).*w + (x*(rdot_L(t)+rdot_R(t))-rdot_L(t))*dw./r(t); %(x*(rdot_L(t)+rdot_R(t))*intS(t)-rdot_L(t))*dw./r(t)
end

function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
    pl = [0;0];
    pr = [0;0];
    ql = [1;1];
    qr = [1;1];
end

sol = pdepe(0, @pdefun, @pdeic, @pdebc, x, t);
u = sol(:,:,1);
v = sol(:,:,2);

%% Plot solution

figure('Color','white')

[X,T] = meshgrid(x,t);

X = r(t)'.*X;

ax1 = subplot(1,3,1);
pcolor(X,T,u)
xlabel('Position $x$', Interpreter='latex', FontSize=18)
ylabel('Time $t$', Interpreter='latex', FontSize=18)
title('System evolution', Interpreter='latex', FontSize=20)
colormap(ax1,parula)
shading interp

%% Plot modes over time

krange = 1:max_mode;
mag_matrix = zeros(length(t),max_mode);
mag_maxs = zeros(length(t),1);
[K,T] = meshgrid(krange,t);

for tt = 1:length(t)
    u_coef = extract_gen_fourier_coff(u(tt,:), x, max_mode);
    v_coef = extract_gen_fourier_coff(u(tt,:), x, max_mode);
    mags = sqrt(u_coef.^2 + v_coef.^2);
    mag_matrix(tt,:) = mags;
    [p,ind] = max(mags);
    mag_maxs(tt) = krange(ind);
end

ax2 = subplot(1,3,2);
pcolor(K,T,mag_matrix)
xlabel('Wavenumber $k$', Interpreter='latex', FontSize=18)
ylabel('Time $t$', Interpreter='latex', FontSize=18)
title('Numerically extracted modes', Interpreter='latex', FontSize=20)
shading flat

%% Plot predicted mode growth over time (Solve trunc ODE system)

Dif = [d_1, 0; 0, d_2];
mode_ICs = [IC1, [a+b; b/(a+b)^2]];

[time_modes,curr_sol] = ode45(@(t,kvec) derivModes(t,kvec), [0,Tf], mode_ICs);
sol = reshape(curr_sol(:,:),[length(time_modes),2,kappa+1]);
mags = squeeze(sqrt(sol(:,1,1:kappa).^2+sol(:,2,1:kappa).^2));
mags_max = max(mags,[],2);
rel_mags = mags./mags_max;

[K,T] = meshgrid(krange,time_modes);

ax3 = subplot(1,3,3);

hold on
pcolor(K,T,rel_mags(:,2:end));
plot(mag_maxs(10:end),t(10:end),'red',LineStyle='-',LineWidth=2);
hold off

xlabel('Wavenumber $k$', Interpreter='latex', FontSize=18)
ylabel('Time $t$', Interpreter='latex', FontSize=18)
axis([1,max_mode,0,Tf]);
title('Linear mode evolution', Interpreter='latex', FontSize=20)
colormap(ax3,'parula')
shading flat

sgtitle(['$\dot{r}_L = -0.03,\, \dot{r}_R = 0.13,\, S_{\mathrm{int}} = 0$'], Interpreter="latex", FontSize=20)

%% Helper functiosn

% Compute RHS of Eq (6.11) in Thm 6
function D = derivModes(t,kvec)

    old = reshape(kvec,[2,kappa+1]);
    auto_sol = old(:,kappa+1:end);
    old = old(:,1:kappa);

    D = zeros(2, kappa+1);
    
    D(:,kappa+1) = kenetics(auto_sol) + S_int(t).*D(:,kappa+1);

    D(:,1) = (Jacob(auto_sol) + S_int(t)*eye(2))*old(:,1) + old(:,2:end)*coef_zero(t,kappa);

    Ma = coef_nonzero_matrix(t,kappa);
    J = Jacob(auto_sol);
    Aval = A(t);
    S_intval = S_int(t);
    rval2 = r(t).^2;
    loop_over = old(:,2:end);
    parfor l=2:kappa
        D(:,l) = (-((l-1)^2 * pi^2./rval2) * Dif + J...
            +(S_intval+Aval/2)*eye(2))*old(:,l) + loop_over*Ma(:,l-1);
    end

    D = reshape(D, [2*kappa+2,1]);
end

% Jacobian evaluated along the homogenous solution (u_0(t),v_0(t))
function J = Jacob(soln)
    u = soln(1);
    v = soln(2);
    J = [-1+2*u*v, u.^2; -2*u*v, -u.^2];
end

% Series coefficents in Eq (6.11a)
function M = coef_zero(t,kappa)
    M = zeros(kappa-1,1);

    for k=1:(kappa-1)
        M(k) = 2*B(t)+(-1)^k*(sqrt(2)*A(t)-2*B(t));
    end
end

% Series coefficents in Eq (6.11b)
function M = coef_nonzero_matrix(t,kappa)
   M = zeros(kappa-1,kappa-1);

   for l=1:(kappa-1)
        for k=1:(kappa-1)
            if k==l
                M(k,l)= 0;
            else
                M(k,l) = (2*k^2/(k^2-l^2))*((-1)^(k+l)*(A(t)+B(t))-B(t));
            end
        end
   end
end

% Growth functions, have hardcoded the precalculated integrals for speed
function z = rdot_L(t)
    z = rate_L; 
end

function z = rdot_R(t)
    z = rate_R; 
end

function z = S_int(t)
    z = rate_int*ones(size(t));
end

function z = r_int(t)
    z = exp(rate_int*t);
end

function I = denom_int(t)
    I = (rate_L+rate_R)*(1-exp(-rate_int*t))/rate_int;
end

function A = A(t)
    A = rate_int*(rate_L+rate_R)./(exp(rate_int*t)*(rate_int+rate_L+rate_R)-rate_L-rate_R); %(rdot_L(t) + rdot_R(t))./(r_int(t).*(denom_int(t)+1));
end

function B = B(t)
    B = -rate_int*rate_L./(exp(rate_int*t)*(rate_int+rate_L+rate_R)-rate_L-rate_R); %-rdot_L(t)./(r_int(t).*(denom_int(t)+1));
end

function z = r(t)
    z = (exp(rate_int*t)*(rate_int+rate_L+rate_R)-rate_L-rate_R)/rate_int;%r_int(t).*(denom_int(t)+1);
end 

function z = rdot(t)
    z = S_int(t)*r(t)+rdot_L(t)+rdot_R(t);
end

end