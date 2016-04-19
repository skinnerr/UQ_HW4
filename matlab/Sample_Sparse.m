function [ m, v, N ] = Sample_Sparse( h, u0, uf, Ki, Pk, x, index_pc )

%%%
% Use tensor-product stochastic collocation (Gauss-Hermite quadrature) to sample sol'n.
%%%

nq = 40;
[pts, w] = mherzo(nq);
N = length(pts)^2;

m = 0;
v = 0;

wsum = 0;

for i = 1:length(pts)
    for j = 1:length(pts)
        x1 = pts(i);
        x2 = pts(j);
        u = Heat_Solution(h, u0, uf, x1, x2, Ki, Pk, x, index_pc);
        m = m + u    * w(i) * w(j);
        v = v + u.^2 * w(i) * w(j);
        wsum = wsum +  w(i) * w(j);
    end
end

m = m' / wsum;
v = v' / wsum;

v = v - m.^2; % Use the fact that Var(X) = E[X^2] - (E[X])^2

end

