function M = uniform_condition_calc(a,b,d_1,d_2,l,kmax, t, y_base, S, r)
% Returns a matrix representing the modes we expect to be unstable at each
% time step, based on linearisation about the steady-state solution y_base

%% Init matricies
M = zeros(length(t), kmax);
LHS_matirx = zeros(length(t), kmax);
RHS_matrixA = zeros(length(t), kmax);
RHS_matrixB = zeros(length(t), kmax);
U=y_base(:,1);
V=y_base(:,2);

%% Calculate LHS of Eq (4.12) in Thm 4
for k=1:kmax
    %% 
    LHS_matirx(:,k) = deter(U,V) + S.*(S+trac(U,V)+(k^2*pi^2./(r.^2*l^2)) * (d_1+d_2))...
        - (k^2*pi^2./(r.^2*l^2)) .* (d_1 * gv(U,V) + d_2 * fu(U,V))...
        + (k^2*pi^2./(r.^2*l^2)).^2 * d_1 * d_2;
end


%% Calculate RHS of Eq (4.12) in Thm 4
d1 = n_Deriv((fu(U,V)+S)./fv(U,V), t);
d2 = n_Deriv(1./(r.^2.*fv(U,V)), t);
d3 = n_Deriv((gv(U,V)+S)./gu(U,V), t);
d4 = n_Deriv(1./(r.^2.*gu(U,V)), t);

RHS_matrixA = (fv(U,V).*d1)*ones(1,kmax) + d_1*(fv(U,V).*d2)*(1:kmax);
RHS_matrixB = (gu(U,V).*d3)*ones(1,kmax) + d_2*(gu(U,V).*d4)*(1:kmax);

RHS_matrix = max(RHS_matrixA,RHS_matrixB);

%% Return True/False matrix for where condition holds
M = (RHS_matrix - LHS_matirx) >0 ;

%% Helper functiosn

    function v = n_Deriv(x,t) % Numerical derivative
        xshift = [0;x(1:end-1)];
        tshift = [-1;t(1:end-1)];
        v = (x-xshift)./(t-tshift);
        v(1) = v(2);
    end

    function d = fu(u,v) % Partial derivs of reaction kinetics 
        d = 2*u.*v -1;
    end

    function d = fv(u,v)
        d = u.^2;
    end

    function d = gu(u,v)
        d = -2*u.*v;
    end

    function d = gv(u,v)
        d = -u.^2;
    end

    function d=deter(u,v) % Det of jacobian
        d= fu(u,v).*gv(u,v) - fv(u,v).*gu(u,v);
    end

    function d=trac(u,v) % Tr of jacobian
        d = fu(u,v) + gv(u,v);
    end
end