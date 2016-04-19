function [ x ] = Thomas( diag, sub, sup, rhs )
    
    %%%%%%
    % Solves a tri-diagonal matrix system using the Thomas algorithm.
    %   diag -- diagonal
    %    sub -- sub-diagonal
    %    sup -- super-diagonal
    %    rhs -- right-hand side vector
    %    sol -- solution vector
    %
    % Ryan Skinner, October 2015
    %%%
    
    % Eliminate the sub-diagonal.
    for i = 1:length(sup)
        r = sub(i) / diag(i);
        diag(i+1) = diag(i+1) - r * sup(i);
         rhs(i+1) =  rhs(i+1) - r * rhs(i);
    end
    
    % Back-substitute and calculate the solution vector.
    x = zeros(length(diag),1);
    x(end) = rhs(end) / diag(end);
    for i = length(diag)-1:-1:1
        x(i) = (rhs(i) - sup(i) * x(i+1)) / diag(i);
    end

end