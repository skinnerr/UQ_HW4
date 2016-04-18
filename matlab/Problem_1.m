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

% Total order of PCE for K_i
pk = 14;

% Calculate the PCE expansion for K_i(x)
[Ki, Pk] = Compute_Ki(pk, sigma, ell, a, d, x);

% Polynomial chaos expansion (PCE) options for solution uhat
p_list = 0:7;
p_list = 1;

% Loop over total solution order (Homework's part 'b')
for p = p_list
    
    % Total number of terms in solution's PC expansion (minus one)
    P = factorial(p+d) / ( factorial(p) * factorial(d) ) - 1;
    
    index_pc = nD_polynomial_array(d, pk);

    %%%
    % Construct the global matrix equation
    %%%
    
    % Sub-matrix dimension
    Nsub = Nx-2;

    K = zeros((Nsub)*(P+1));
    f = zeros((Nsub)*(P+1), 1);
    
    for i = 1:(P+1)
    for j = 1:i
        
        % Build the K_ij matrix
        Kij = zeros(Nsub);
        for k = 1:(Pk+1)
        
            % LHS sub-matrix and RHS sub-vector (K_ij and f_i)
            [dag, sub, sup, rhs] = Assemble_u(dx, u0, uf, Ki(k,:));

            % Build dense sub-matrix from sparse storage
            Kk = zeros(Nsub);
            Kk = Kk + diag(dag);
            Kk = Kk + [[zeros(1,length(sub));diag(sub)],zeros(length(sub)+1,1)];
            Kk = Kk + [zeros(length(sup)+1,1),[diag(sup);zeros(1,length(sup))]];
            
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
        if i ~= j
            K(Brng,Arng) = Kij;
        end
    
    end
    end
    
    %%%
    % Solve the global matrix equation
    %%%
    
    u = f \ K;
    u = reshape(u, Nsub, [])';
    plot(u(1,:));
    
end

end





















