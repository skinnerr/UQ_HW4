function [ l, phix ] = Analytical_Eigs( sigma, ell, a, N, x )
%%%
% IN:
%      sigma - standard deviation of correlation function
%      ell   - correlation length parameter
%      a     - support of eigenproblem is [-a,a]
%      N     - # of eigenpairs to compute
%      x     - points to calculate eigenfunctions
% OUT:
%      l     - eigenvalues
%      phix  - eigenfunctions evaluated over x
%%%

% Pre-allocate variables.
omega = nan(N,1);
l     = nan(N,1);
phix  = nan(N,length(x));

% Number of roots of each system to find.
max_even = floor(N/2);
max_odd  = floor((N+1)/2);

% Larger values make MATLAB happier, but sacrifice root finding interval accuracy.
eps_scale = 10000;

%%%%%%%%%%%%%%%%%%
% Compute omegas %
%%%%%%%%%%%%%%%%%%

% Odd.
fun = @(omega) (1/ell) - omega * tan(omega * a);
for i = 1:max_odd
    range = pi/(2*a) * [max(0,(i*2)-3), (i*2)-1] + eps_scale*[eps,-eps];
    omega(1+2*(i-1)) = fzero(fun, range);
end

% Even.
fun = @(omega) omega + (1/ell) * tan(omega * a);
for i = 1:max_even
    range = pi/(2*a) * [i*2-1, i*2+1] + eps_scale*[eps,-eps];
    omega(i*2) = fzero(fun, range);
end

%%%%%%%%%%%%%%%%%%%%%%%
% Compute eigensystem %
%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:N
    
    % Odd eigenfunctions.
    if mod(i,2)
        phix(i,:) = cos(omega(i)*x) ./ sqrt(a + sin(2*a*omega(i)) / (2*omega(i)));
        
    % Even eigenfunctions.
    else
        phix(i,:) = sin(omega(i)*x) ./ sqrt(a - sin(2*a*omega(i)) / (2*omega(i)));
    end
    
    % Eigenvalues use same formula for both odd and even.
    l(i) = 2*ell / (ell^2 * omega(i)^2 + 1);
    
end

% Correct for sigma dependence in eigenvalues.
l = l * sigma^2;

% Transpose so index is (i,x).
phix = phix';

end







