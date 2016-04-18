function [] = Problem_1()

Set_Default_Plot_Properties();

% Solution domain.
Nx = 11;
x0 = 0;
xf = 1;
x = linspace(x0, xf, Nx)';
dx = x(2) - x(1);

% Solution boundary conditions.
u0 = 0;
uf = 0;

% Karhunen-Loeve expansion (KLE) options
sigma = 2.0;    % Standard deviation
ell = 2.0;      % Correlation length
a = 1/2;        % Support of eigenproblem
d = 2;          % Number of terms

% Polynomial chaos expansion (PCE) options for K_i
pk = 14;         % Total order

% Polynomial chaos expansion (PCE) options for solution uhat
p_list = 0:7;

% Calculate the PCE expansion for K_i(x)
[Ki, Pk] = Compute_Ki(pk, sigma, ell, a, d, x);

% Loop over total solution order (Homework's part 'b')
for p = p_list
    
    % Total number of terms in solution's PC expansion (minus one)
    P = factorial(p+d) / ( factorial(p) * factorial(d) ) - 1

    %%%
    % Construct the global matrix equation.
    %%%

    K = zeros((Nx-2)*(P+1));
    f = zeros((Nx-2)*(P+1), 1);
    
    for i = 1:(P+1)
    for j = 1:i
        
        % Build the K_ij matrix
        Kij = zeros(Nx-2);
        for k = 1:(Pk+1)
        
            % LHS sub-matrix and RHS sub-vector (K_ij and f_i)
            [dag, sub, sup, rhs] = Assemble_u(dx, u0, uf, Ki);

            % Build dense sub-matrix from sparse storage
            Kk = zeros(Nx-2);
            Kk = Kk + diag(dag);
            Kk = Kk + [[zeros(1,length(sub));diag(sub)],zeros(length(sub)+1,1)];
            Kk = Kk + [zeros(length(sup)+1,1),[diag(sup);zeros(1,length(sup))]];
            
            % Add contribution to Kij
            Kij = Kij + Kk * cijk(i,j,k);
    
    end
    end
    
end

end





















