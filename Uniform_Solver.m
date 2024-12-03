function Uniform_Solver
% Function produces the plots for the hybrid growths (solve the PDE + numerical modes + Thm 4)

%% Parameters
a=0.01;
b=1;
d_1 = 0.004;
d_2 = 0.1;
L=1;
Tf=170;

%% Growth functions
function z = S(t)
    z = 0.053-0.00064*t';%().*ones(1,length(t));
end

function z = r(t)
    z = exp(0.053*t-0.00032*t.^2);
end

%% Solve PDE
function zvec = kenetics(w)
    zvec = [a-w(1) + w(1)^2 * w(2); b - w(1)^2 * w(2)];
end

function w0 = pdeic(x)
    u_0 = a+b;
    v_0 = b/(a+b)^2;
    w0 = [u_0; v_0]+0.05*(rand(1)-0.5)*[1;1];
end

function [c,f,s] = pdefun(x,t,w,dw)

    c = [1;1];
    f = (1/r(t).^2)*[d_1; d_2].* dw;
    s = kenetics(w) + S(t)*w;
end

function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
    pl = [0;0];
    pr = [0;0];
    ql = [1;1];
    qr = [1;1];
end

x = linspace(0,L,401);
t = linspace(0,Tf,1001);

sol = pdepe(0, @pdefun, @pdeic, @pdebc, x, t);
u = sol(:,:,1);
v = sol(:,:,2);

%% Plot solution

figure('color','white')

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
max_mode = 45;
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

%% Analytic instability (Plot result of Thm 4)

% Solve for homogenous state (u_0(t),v_0(t))
IC = [(a+b); b/(a+b)^2];
[t_base,y_base] = ode45(@(t,y) [derX(t,y(1),y(2),a,b);derY(t,y(1),y(2),a,b)], [0,t(end)], IC);

% Use Thm 4 to predict stability
stability_matrix = uniform_condition_calc(a,b,d_1,d_2,L,...
    max_mode, t_base, y_base, S(t_base)', r(t_base));

% Plot
ax3 = subplot(1,3,3);
[K,T] = meshgrid(krange,t_base);
hold on
pcolor(K,T,1-double(stability_matrix))
plot(mag_maxs(10:end),t(10:end),'red',LineStyle='-',LineWidth=2);
hold off
xlabel('Wavenumber $k$', Interpreter='latex', FontSize=18)
ylabel('Time $t$', Interpreter='latex', FontSize=18)
axis([1,max_mode,0,Tf]);
title('Analytic unstable modes', Interpreter='latex', FontSize=20)
colormap(ax3, gray(2))
shading flat

sgtitle(['$S(t) = 0.053-0.00064t$'], Interpreter="latex", FontSize=20)

%% Helper functions (Appying THM 4 function in uniform_condition_calc.m)

function du = derX(t,u,v,a,b)
    du = a - u + u.^2 .* v + S(t).*u;
end

function dv = derY(t,u,v,a,b)
    dv = b - u.^2 .* v +S(t).*v;
end

end