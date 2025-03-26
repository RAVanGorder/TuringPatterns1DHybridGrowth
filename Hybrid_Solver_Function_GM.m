function Hybrid_Solver_Function_GM(node_initals, s_dots, S_ints, final_time, time_res, space_res, max_k)
% Function obtains solutions and plots for general multi-compartment hybrid evolution
% where Gierer–Meinhardt reaction kenetics are used.

%% Parameters
% Reaction-kinetics and diffusion parameters for reaction-diffusion system
a=1.4;       
b=5;
d_1 = 0.004;
d_2 = 0.4;

function zvec = kenetics(w) %  Gierer–Meinhardt reaction kenetics
    zvec = [1+w(1)^2/w(2)-a*w(1); w(1)^2-b*w(2)];
end

% Jacobian matrix of rection kenetics
function M = J_kenetics(w)
    M = [2*w(1)/w(2)-a, -(w(1)/w(2))^2; 2*w(1), -b];
end

% Homogenous equilibrium of reaction kenetics
IC_nopb = [(b+1)/a; (b+1)^2/(b*a^2)];

% Arrays of growth functions
R_bar_vec = node_initals; % Apical growth node locations at t=0 (starts at R_0)
s_dot_vec = s_dots; % Limiting apical growth (starts at s_0)
S_int_vec = S_ints; % Uniform growth in each compartment (starts at S_int,1)
N = length(S_int_vec); % N 'compartment(s)', N+1 apical growth nodes

% Numerical simulation parameters
Tf=final_time;
max_mode = max_k;
x = linspace(0,1,space_res);
t = linspace(0,Tf,time_res);
t_nodes = linspace(0,Tf,time_res*10);

% Precalculate the location of apical nodes over time
node_locations = calc_node_locations(t_nodes);
max_size = max(node_locations(N+1,:));

%% Solve the reaction-diffusion PDE system

% Set the inital conditions
function u0 = pdeic(xi) 
    u0 = IC_nopb + 0.05*(rand-0.5);
end

% Internal dynamics of PDE
function [c,f,s] = pdefun(xi,tt,w,dw) 
    c = [1;1];
    f = (1/r(tt,N).^2)*[d_1; d_2].* dw;
    s = kenetics(w) - S(xi,tt).*w + chi(xi,tt).*dw;
end

% Boundary Conditions for PDE (Neumann, no flux)
function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,tt) 
    pl = [0;0];
    pr = [0;0];
    ql = [1;1];
    qr = [1;1]; 
end

% Simulation of the reaction-diffusion system using method of linear and FDM in space
sol = pdepe(0, @pdefun, @pdeic, @pdebc, x, t);
u = sol(:,:,1);
v = sol(:,:,2);

%% Plot of the u component from the reaction-diffusion system
figure('Color','white')
[X,T] = meshgrid(x,t);

% Display domain evolution
r_vec = zeros(size(t));
for ii =1:length(t)
    r_vec(ii) = r(t(ii),N);
end
X = r_vec'.*X;

ax1 = subplot(1,4,1);

hold on
pcolor(X,T,u)
colorbar('eastoutside');
clim([0 10]);

% Lines describing locations of apical nodes over time
for i=1:N+1
    plot(node_locations(i,:),t_nodes,'red',LineStyle='-',LineWidth=2);
end
hold off

xlabel('Position $x$', Interpreter='latex', FontSize=18)
ylabel('Time $t$', Interpreter='latex', FontSize=18)
title('System evolution', Interpreter='latex', FontSize=20)
axis([0,max_size,0,Tf]);
colormap(ax1,parula)
shading interp

%% Extract and plot spatial modes from the numerical simulation
% Base state in each compartment
base_each_compartment = {};

% Solve ODE for time varying base state in each compartment
for ii = 1:N
    [~,sol_base] = ode45(@(t,w) kenetics(w)-S_int_vec{ii}(t)*w, t, IC_nopb);
    base_each_compartment{ii} = sol_base;
end

% Stitch together each compartment's base state 
base = zeros(length(t), length(x),2);
for tt_i=1:length(t)
    for xi_i = 1:length(x)
        xx = x(xi_i)*r(t(tt_i),N);
        for ii=1:N 
            if xx <= r(t(tt_i),ii)
                base(tt_i,xi_i,:) = base_each_compartment{ii}(tt_i,:);
                break
            end
        end
    end
end

u_base = base(:,:,1);
v_base = base(:,:,2);

% Determine coefficents of spatial modes of (u-u_base, v-v_base)
krange = 1:max_mode;
mag_matrix = zeros(length(t),max_mode);
mag_maxs = zeros(length(t),1);
[K,T] = meshgrid(krange,t);

for tt = 1:length(t)
    u_coef = extract_gen_fourier_coff(u_base(tt,:)-u(tt,:), x, max_mode);
    v_coef = extract_gen_fourier_coff(v_base(tt,:)-v(tt,:), x, max_mode);
    mags = sqrt(u_coef.^2 + v_coef.^2);
    mag_matrix(tt,:) = mags;
    [~,ind] = max(mags);
    mag_maxs(tt) = krange(ind);
end

% Plot numerically calculated coefficients of spatial modes
ax2 = subplot(1,4,2);
pcolor(K,T,mag_matrix);
xlabel('Wavenumber $k$', Interpreter='latex', FontSize=18)
ylabel('Time $t$', Interpreter='latex', FontSize=18)
title('Numerical modes', Interpreter='latex', FontSize=20)
shading flat


%% Solve linear instability system from Corollary 5.1
D = [d_1,0;0,d_2];
IC_modes = 0.001*(rand(2,max_mode+1)-0.5);
IC_compartments = [IC_nopb(1)*ones(1,N);IC_nopb(2)*ones(1,N)];
IC_colvec = reshape([IC_modes, IC_compartments], [2*(max_mode+1+N),1]);

% Use ode113 for accuracy when considering purely uniform growth, and
% ode15s for all other cases.
[t_modes,modes_sol] = ode15s(@(tt,kvec) derivModes(tt,kvec), [0,Tf], IC_colvec);

[K,T] = meshgrid(krange,t_modes);


%% Plot solution to linear instability system

% Plot mode coefficents, normalised at each time step
sol = reshape(modes_sol(:,:),[length(t_modes),2,max_mode+1+N]);
mags = squeeze(sqrt(sol(:,1,2:max_mode+1).^2+sol(:,2,2:max_mode+1).^2));
mags_max = max(mags,[],2);
rel_mags = mags./mags_max;

ax3 = subplot(1,4,3);

hold on
pcolor(K,T,rel_mags);
plot(mag_maxs(10:end),t(10:end),'red',LineStyle='-',LineWidth=2);
hold off

xlabel('Wavenumber $k$', Interpreter='latex', FontSize=18)
ylabel('Time $t$', Interpreter='latex', FontSize=18)
axis([1,max_mode,0,Tf]);
title('Linear mode evolution', Interpreter='latex', FontSize=20)
colormap(ax3,'parula')
shading flat

% Plot relative growth rates for each mode coefficent
logderiv_mags = zeros(size(mags));
for modek=1:(max_mode)
    logderiv_mags(:,modek) = smoothdata(n_Deriv(mags(:,modek),t_modes)./mags(:,modek),'movmean',5); %smoothdata(n_Deriv(log(mags(:,modek)),t_modes),'movmean',5);
end

logderiv_mags_binary = logderiv_mags>0;
mags_big_enough = rel_mags>10^-16;
logderiv_mags_positive = logderiv_mags_binary.*logderiv_mags.*mags_big_enough;

% Find appropriate scale for colours
max_colour = prctile(logderiv_mags_positive(logderiv_mags_binary),93,"all");

ax4 = subplot(1,4,4);

hold on
pcolor(K,T,logderiv_mags_positive);%pcolor(K,T,rel_mags);%
plot(mag_maxs(10:end),t(10:end),'red',LineStyle='-',LineWidth=2);
clim([0,max_colour])
hold off

colorbar('eastoutside');
xlabel('Wavenumber $k$', Interpreter='latex', FontSize=18)
ylabel('Time $t$', Interpreter='latex', FontSize=18)
axis([1,max_mode,0,Tf]);
title('Linear mode growth', Interpreter='latex', FontSize=20)
colormap(ax4,'parula')
shading flat


%% Functions for useful quantities

% Dynamics relevant for Corollary 5.1
function dkvec = derivModes(tt,kvec)
    old = reshape(kvec, [2,max_mode+1+N]);
    old_modes = old(:,1:max_mode+1);
    old_compartments = old(:,max_mode+2:max_mode+1+N);

    new_modes = zeros(2,max_mode+1);
    new_compartments = zeros(2,N);

    r_vec = arrayfun(@(nn) r(tt,nn), 0:N);
    r_now = r_vec(N+1);
    coefs = get_matrix_coefs(tt,old_compartments,r_vec/r_now);   

    % Deriv of U and V in each compartment base state
    for ii=1:N
        new_compartments(:,ii) = kenetics(old_compartments(:,ii)) - S_int_vec{ii}(tt)*old_compartments(:,ii);
    end

    % Terms relevant for spatial modes
    for k = 0:max_mode
        new_modes(:,k+1) = -(k*pi/r_now)^2*D*old_modes(:,k+1);
        for l = 0:max_mode
            new_modes(:,k+1) = new_modes(:,k+1)+coefs(:,:,k+1,l+1)*old_modes(:,l+1);
        end
    end

    dkvec = reshape([new_modes, new_compartments], [2*(max_mode+1+N),1]);
end

% Compute the matrix coefficents M_{k,l}(tt) in Corollary 5.1
function coefs = get_matrix_coefs(tt, cmps, seg_ends)
    res = 101; % can be increased (to improve accuracy of simulations) or decreased (to improve speed of simulations)

    coefs = zeros(2,2,max_mode+1,max_mode+1);

    xvecs = cell(1,N);
    J11s = cell(1,N);
    J12s = cell(1,N);
    J21s = cell(1,N);
    J22s = cell(1,N);
    Ss   = cell(1,N);
    for jj=1:N
        xvecs{jj} = linspace(seg_ends(jj),seg_ends(jj+1),res);
        A = J_kenetics(cmps(:,jj));
        J11s{jj} = A(1,1)*ones(1,res);
        J12s{jj} = A(1,2)*ones(1,res);
        J21s{jj} = A(2,1)*ones(1,res);
        J22s{jj} = A(2,2)*ones(1,res);
        Ss{jj}   = S_int_vec{jj}(tt)*ones(1,res);
    end
    
    xvec = [xvecs{:}];
    J11 = [J11s{:}];
    J12 = [J12s{:}];
    J21 = [J21s{:}];
    J22 = [J22s{:}];
    S = [Ss{:}];

    chi_vec = arrayfun(@(xx) chi(xx, tt), xvec);

    for k=0:max_mode
        for l=0:max_mode
            if k == 0 && l == 0
                coefs(1,1,k+1,l+1) = trapz(xvec, J11-S);
                coefs(1,2,k+1,l+1) = trapz(xvec, J12);
                coefs(2,1,k+1,l+1) = trapz(xvec, J21);
                coefs(2,2,k+1,l+1) = trapz(xvec, J22-S);
            elseif k >= 1 && l == 0
                coefs(1,1,k+1,l+1) = sqrt(2)*trapz(xvec, (J11-S).*cos(k*pi*xvec));
                coefs(1,2,k+1,l+1) = sqrt(2)*trapz(xvec, (J12).*cos(k*pi*xvec));
                coefs(2,1,k+1,l+1) = sqrt(2)*trapz(xvec, (J21).*cos(k*pi*xvec));
                coefs(2,2,k+1,l+1) = sqrt(2)*trapz(xvec, (J22-S).*cos(k*pi*xvec));
            elseif k == 0 && l >= 1
                coefs(1,1,k+1,l+1) = sqrt(2)*trapz(xvec, (J11-S).*cos(l*pi*xvec)-l*pi*chi_vec.*sin(l*pi*xvec));
                coefs(1,2,k+1,l+1) = sqrt(2)*trapz(xvec, (J12).*cos(l*pi*xvec));
                coefs(2,1,k+1,l+1) = sqrt(2)*trapz(xvec, (J21).*cos(l*pi*xvec));
                coefs(2,2,k+1,l+1) = sqrt(2)*trapz(xvec, (J22-S).*cos(l*pi*xvec)-l*pi*chi_vec.*sin(l*pi*xvec));
            else
                coefs(1,1,k+1,l+1) = sqrt(2)*trapz(xvec, ((J11-S).*cos(l*pi*xvec)-l*pi*chi_vec.*sin(l*pi*xvec)).*cos(k*pi*xvec));
                coefs(1,2,k+1,l+1) = sqrt(2)*trapz(xvec, (J12).*cos(l*pi*xvec).*cos(k*pi*xvec));
                coefs(2,1,k+1,l+1) = sqrt(2)*trapz(xvec, (J21).*cos(l*pi*xvec).*cos(k*pi*xvec));
                coefs(2,2,k+1,l+1) = sqrt(2)*trapz(xvec, ((J22-S).*cos(l*pi*xvec)-l*pi*chi_vec.*sin(l*pi*xvec)).*cos(k*pi*xvec));
            end
        end
    end
end

% Precompute positions of apical nodes, jj=N for r
function z = calc_node_locations(req_times)
    z = zeros(N+1, length(req_times)); 
    for ii=1:length(req_times)
        tvec = req_times(1:ii); % Integration Domain
        if length(tvec)== 1
            tvec = [tvec,tvec];
        end

        % First node location
        int_strain_1 = exp(cumtrapz(tvec, S_int_vec{1}(tvec)));
        integrand = (s_dot_vec{1}(tvec) + s_dot_vec{2}(tvec))./(int_strain_1);
        z(2,ii) = int_strain_1(end) * (R_bar_vec(2) + trapz(tvec, integrand));

        % 2nd through Nth node locations (each uses previous locations)
        for jj = 2:N
            int_strain_jj = exp(cumtrapz(tvec, S_int_vec{jj}(tvec)));
            sum_on_top = zeros(size(tvec));
            for ell=1:jj-1
                sum_on_top = sum_on_top + (2*s_dot_vec{ell+1}(tvec) + S_int_vec{ell}(tvec).*(z(ell+1,1:ii)-z(ell,1:ii)));
            end   
            integrand = (s_dot_vec{1}(tvec) + s_dot_vec{jj+1}(tvec) - S_int_vec{jj}(tvec).*z(jj,1:ii) + sum_on_top)./(int_strain_jj);
            z(jj+1,ii) = int_strain_jj(length(tvec)) * (R_bar_vec(jj+1) + trapz(tvec, integrand));
        end
    end
end

% Interpolate between precalculated values
function z = r(tt, apical_node)
    z = interp1(t_nodes, node_locations(apical_node+1,:),tt);
end

% Piecewise linear coefficents for Chi
function [z1,z2] = script_XY_part(tt,jj)

    big_sum_1 = 0;
    big_sum_2 = 0;
    for ell=1:jj-1
        big_sum_2 = big_sum_2 + (2*s_dot_vec{ell+1}(tt) + S_int_vec{ell}(tt).*(r(tt,ell) - r(tt,ell-1)));
    end 
    for ell=1:N-1
        big_sum_1 = big_sum_1 + (2*s_dot_vec{ell+1}(tt) + S_int_vec{ell}(tt).*(r(tt,ell) - r(tt,ell-1)));
    end

    z1 = (s_dot_vec{1}(tt)+s_dot_vec{N+1}(tt) + S_int_vec{N}(tt).*(r(tt,N) - r(tt,N-1))-S_int_vec{jj}(tt)*r(tt,N)+big_sum_1)/r(tt,N);
    z2 = (s_dot_vec{1}(tt)-S_int_vec{jj}(tt)*r(tt,jj-1)+big_sum_2)/r(tt,N); 
end

% Chi coefficents
function [z1,z2] = script_XY(xi,tttt)
    flag = 0;
    for iiii=1:N-1
        if xi <= r(tttt,iiii)/r(tttt,N)
            [z1,z2] = script_XY_part(tttt,iiii);
            flag = 1;
            break
        end
    end
    if flag == 0
        [z1,z2] = script_XY_part(tttt,N);
    end
end

% Advection term
function z = chi(xi, tt)
    [X_f,Y_f] = script_XY(xi,tt);
    z = xi*X_f - Y_f;
end

% Strain term
function strain = S(xi,tt)
    flag = 0;
    xx = xi*r(tt,N);
    for ii=1:N-1
        if xx <= r(tt,ii)
            strain = S_int_vec{ii}(tt);
            flag=1;
            break
        end
    end
    if flag == 0
        strain = S_int_vec{N}(tt);
    end
end

% Numerical calculation of derivative
function v = n_Deriv(x,t) 
    v = zeros(size(x));
    v(1) = (x(2)-x(1))/(t(2)-t(1)); % Forward difference
    for kk = 1:(length(x)-1)
        v(kk) = (x(kk+1)-x(kk))/(t(kk+1)-t(kk)); % Forward difference
    end
    v(end) = (x(end)-x(end-1))/(t(end)-t(end-1)); % Backwards difference
end
end