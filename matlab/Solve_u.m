function [ sol ] = Solve_u( h, u0, uf, K )

    [diag, sub, sup, rhs] = Assemble_u( h, u0, uf, K );
    sol = Thomas(diag, sub, sup, rhs);

end

