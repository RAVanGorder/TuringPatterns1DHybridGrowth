function Apical_Solver
% Function produces the plots for the apical growths (solve the PDE + numerical modes + linear system)

%% Parameters
% Params for model
a=0.01;
b=1;
d_1 = 0.004;
d_2 = 0.1;
L=1;

% Params for solution
Tf=100;
max_mode = 60;
kappa = max_mode+1;
x = linspace(0,1,501);
t = linspace(0,Tf,1001);

% Inital conditions for linear ODE system
IC1 = 0.0001*(rand(2,kappa)-0.5);

%% Solve PDE

function zvec = kenetics(w)
    zvec = [a-w(1) + w(1)^2 * w(2); b - w(1)^2 * w(2)];
end

function u0 = pdeic(xx) % Reconstruct appropriate ICs from mode ICs
    u0 = [(a+b); b/(a+b)^2];
    for kk=0:max_mode
        u0 = u0 + (IC1(:,kk+1)*sqrt(2)*cos(pi*kk*xx));
    end
end

function [c,f,s] = pdefun(x,t,w,dw)

    c = [1;1];
    f = (1/r(t).^2)*[d_1; d_2].* dw;
    s = kenetics(w) + x*rdot(t)*dw/r(t);
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

%% Plot numeric modes over time

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

%% Plot predicted mode growth over time (Solve truncated ODE system)

J = [(b-a)/(a+b), (a+b)^2; -2*b/(a+b), -(a+b)^2];
Dif = [d_1, 0; 0, d_2];
[time_modes,curr_sol] = ode45(@(t,kvec) derivModes(t,kvec,J,Dif,kappa,L), [0,Tf], IC1);

% Convert solution into relative magnitudes over time
sol = reshape(curr_sol(:,:),[length(time_modes),2,kappa]);
mags = squeeze(sqrt(sol(:,1,:).^2+sol(:,2,:).^2));
mags_max = max(mags,[],2);
rel_mags = mags./mags_max;


% Plot
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

sgtitle(['$r(t) = \exp(0.04t)$'], Interpreter="latex", FontSize=20)

end

%% Helper functions

function D = derivModes(t,kvec,J,Dif,kappa,L) % RHS of Eqn (5.14) in Thm 5

    old = reshape(kvec,[2,kappa]);
    D = zeros(2,kappa);

    D(:,1) = J*old(:,1) + sqrt(2)*rdot(t)./r(t) * sum((-1).^(1:(kappa-1)) .*old(:,2:end),2);

    coefs = getCoefs(kappa);
    for l=2:kappa
        D(:,l) = (-((l-1)^2 * pi^2./(r(t).^2 * L^2)) * Dif + J...
            + rdot(t)./(2*r(t)) .*[1,0;0,1])*old(:,l) +...
            2*rdot(t)./r(t) * old(:,2:end)*coefs(l-1,:)';
    end

    D = reshape(D, [2*kappa,1]);
end

% Prescribed domain growth (and deriv)
function z = r(t)
    z = 1+0.05*t;
end

function z = rdot(t)
    z = 0.05; 
end

% Coefficents in Eqn (5.14) series
function M = getCoefs(max)
    M = zeros(max-1,max-1);

    for l=1:(max-1)
        for k=1:(max-1)
            if k==l
                M(l,k) = 0;
            else
                M(l,k) = (k^2 *(-1)^(k+l))/(k^2-l^2);
            end
        end
    end
    
end