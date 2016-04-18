function [ K, Pk ] = Compute_Ki( pk, sigma, ell, a, d, x )
%%%
% Computes Pk shifted lognormal polynomial chaos expansion basis functions
%%%

% KL expansion eigensystem of G(x, omega)
[ l, phix ] = Analytical_Eigs( sigma, ell, a, d, linspace(-a, a, length(x)));

% Order of Phi_i along the direction j, written as i_j
index_pc = nD_polynomial_array(d, pk);

% Total number of terms in PCE expansion (minus one)
Pk = size(index_pc, 1) - 1;

%%%
% PCE expansion's K_i matrices for the thermal coefficient K = 2 + exp(G)
%%%

% Initialize K (note: index goes from 0..Pk, so size is Pk+1)
K = nan(Pk+1, length(x));

% K_0
tmp = zeros(length(x),1);
for j = 1:d
    tmp = tmp + l(j) * phix(:,j).^2;
end
expsum = exp(1 + 0.5*tmp);
K(1,:) = 2 + expsum;

% Subsequent K_i's
for i = 2:(Pk+1)
    numprod = ones(length(x),1);
    for j = 1:d
        numprod = numprod .* ( sqrt(l(j)) * phix(:,j) ).^index_pc(i,j);
    end
    denomprod = ones;
    for j = 1:d
        denomprod = denomprod * factorial(index_pc(i,j));
    end
    K(i,:) = expsum .* numprod / sqrt(denomprod);
end

end