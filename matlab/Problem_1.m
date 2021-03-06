function [] = Problem_1()

Set_Default_Plot_Properties();

% Solution domain
Nx = 101;
x0 = 0;
xf = 1;
x = linspace(x0, xf, Nx)';
dx = x(2) - x(1);

% Solution boundary conditions
u0 = 0;
uf = 0;

% Karhunen-Loeve expansion (KLE) options
sigma = 2.0;    % Standard deviation
ell = 2.0;      % Correlation length
a = 1/2;        % Support of eigenproblem
d = 2;          % Number of terms

% Total order of PCE for K_i
pk = 14;
index_pc = nD_polynomial_array(d, pk);

% Calculate the PCE expansion for K_i(x)
[Ki, Pk] = Compute_Ki(pk, sigma, ell, a, d, x);

%%%
% Compute the "true" solution based on sparse sampling
%%%

[m_true, v_true, ~] = Sample_Sparse(dx, u0, uf, Ki, Pk, x, index_pc);

figure()

subplot(2,1,1)
plot(x, [nan, m_true, nan])
xlabel('x');
ylabel('mean');

subplot(2,1,2)
plot(x, [nan, v_true, nan])
xlabel('x');
ylabel('variance');

%%%
% Compute the full PCE solution
%%%

% Figures
hfig1 = figure();
hfig2 = figure();

% Loop over total solution order
p_list = 0:7;
re_m_pt5 = nan(length(p_list),1);
re_v_pt5 = nan(length(p_list),1);
for p_ind = 1:length(p_list)
    
    p = p_list(p_ind);
    
    fprintf('Working on solution of degree p = %i\n', p);
    
    % Total number of terms in solution's PC expansion (minus one)
    P = factorial(p+d) / ( factorial(p) * factorial(d) ) - 1;

    %%%
    % Construct the global matrix equation
    %%%
    
    % Sub-matrix dimension
    Nsub = Nx-2;

    K = zeros((Nsub)*(P+1));
    f = zeros((Nsub)*(P+1), 1);
    
    for i = 1:(P+1)
    for j = 1:(P+1)
        
        % Build the K_ij matrix
        Kij = zeros(Nsub);
        for k = 1:(Pk+1)
        
            % LHS sub-matrix and RHS sub-vector (K_ij and f_i)
            [dag, sub, sup, rhs] = Assemble_u(dx, u0, uf, Ki(k,:));

            % Build dense sub-matrix from sparse storage
            Kk = zeros(Nsub);
            Kk = Kk - diag(dag);
            Kk = Kk - [[zeros(1,length(sub));diag(sub)],zeros(length(sub)+1,1)];
            Kk = Kk - [zeros(length(sup)+1,1),[diag(sup);zeros(1,length(sup))]];
            
            % Add contribution to Kij
            Kij = Kij + Kk * cijk_hermite(i, j, k, index_pc);
            
            % Set f
            if i == 1
                f(1:Nsub,1) = rhs;
            end
            
        end
        
        A = 1 + (i-1)*(Nsub);
        B = 1 + (j-1)*(Nsub);
        Arng = A:(A+(Nsub-1));
        Brng = B:(B+(Nsub-1));
        
        K(Arng,Brng) = Kij;
    
    end
    end
    
    %%%
    % Solve the global matrix equation
    %%%
    
    u = K\f;
    u = reshape(u, Nsub, [])';
    
    %%%
    % Compute statistics
    %%%
    
    m = u(1,:);         % Mean
    if p > 0
        v = sum(u.^2);  % Second moment
    else
        v = u.^2;
    end
    v = v - m.^2; % Use the fact that Var(X) = E[X^2] - (E[X])^2
    
    m_err = abs(m - m_true) ./ m_true;
    v_err = abs(v - v_true) ./ v_true;
    
    re_m_pt5(p_ind) = m_err(floor(Nx/2));
    re_v_pt5(p_ind) = v_err(floor(Nx/2));
    
    m = [nan, m, nan];
    v = [nan, v, nan];
    m_err = [nan, m_err, nan];
    v_err = [nan, v_err, nan];
    
    %%%
    % Plot solution
    %%%
    
    figure(hfig1);
    
    subplot(2,1,1);
    hold on;
    plot(x, m', 'DisplayName', sprintf('p = %i', p));
    xlabel('x');
    ylabel('mean');
    
    subplot(2,1,2);
    hold on;
    plot(x, v', 'DisplayName', sprintf('p = %i', p));
    xlabel('x');
    ylabel('variance');
    
    figure(hfig2);
    
    subplot(2,1,1);
    hold on;
    plot(x, m_err', 'DisplayName', sprintf('p = %i', p));
    xlabel('x');
    ylabel('rel. err. in mean');
    
    subplot(2,1,2);
    hold on;
    plot(x, v_err', 'DisplayName', sprintf('p = %i', p));
    xlabel('x');
    ylabel('rel. err. in variance');
    
end

figure(hfig1);
subplot(2,1,1);
hleg = legend('show');
set(hleg, 'Location', 'eastoutside');
subplot(2,1,2);
hleg = legend('show');
set(hleg, 'Location', 'eastoutside');

figure(hfig2);
subplot(2,1,1);
set(gca, 'YScale', 'log');
hleg = legend('show');
set(hleg, 'Location', 'eastoutside');
subplot(2,1,2);
set(gca, 'YScale', 'log');
hleg = legend('show');
set(hleg, 'Location', 'eastoutside');

figure();
hold on;
plot(p_list, re_m_pt5, 'DisplayName', 'mean');
plot(p_list, re_v_pt5, 'DisplayName', 'variance');
xlabel('p');
ylabel('rel. err.');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

end





















