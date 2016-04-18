function [dag, sub, sup, rhs] = Assemble_u( h, u0, uf, K )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the theta-system.
    %   diag -- diagonal
    %    sub -- sub-diagonal
    %    sup -- super-diagonal
    %    rhs -- right-hand side vector
    %%%
    
    N = length(K);
    
    dag_range = 2:N-1;
    sub_range = 3:N-1;
    sup_range = 2:N-2;
    
    dag = (-2/h^2) * K(dag_range);
    sub = ( 1/h^2) * ( K(sub_range) - K(sub_range+1)/4 + K(sub_range-1)/4);
    sup = ( 1/h^2) * ( K(sup_range) + K(sup_range+1)/4 - K(sup_range-1)/4);
    rhs = ones(length(dag_range),1);
    
    % Account for boundary conditions, even though thf = 0.
    rhs(1)   = rhs(1)   - (( 1/h^2) * ( K(2)   - K(2+1)  /4 + K(2-1)  /4)) * u0;
    rhs(end) = rhs(end) - (( 1/h^2) * ( K(N-1) + K(N-1+1)/4 - K(N-1-1)/4)) * uf;
    
end