function [ u ] = Heat_Solution( h, u0, uf, x1, x2, Ki, Pk, x, index_pc )

%%%
% Returns the solution to the heat equation.
%%%

% Generate realization of K_{Pk}
K = zeros(1,length(x));
Phi_i = piset_hermite([x1, x2], index_pc);
for i = 1:Pk+1
    K = K - Ki(i,:) * Phi_i(i);
end

% Solve heat equation with sampled inputs.
u = Solve_u(h, u0, uf, K);

end

