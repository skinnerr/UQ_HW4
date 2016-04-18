
% This is a code to calculate <\psi_i psi_j psi_k> for Hermite polynomials

% Assumption: The basis functions \psi_i is normalized such that
% E[\psi_i^2] = 1

function cijk = cijk_hermite(i,j,k,index_pc)

cijk = 1;
idim = 1;
while cijk ~= 0 && idim<= size(index_pc,2)
    s = index_pc(i,idim) + index_pc(j,idim) + index_pc(k,idim);
    if mod(s,2)==0 && max([index_pc(i,idim),index_pc(j,idim),index_pc(k,idim)])<=s/2
        s = s/2;
        cijk_idim = sqrt(factorial(index_pc(i,idim))*factorial(index_pc(j,idim))*factorial(index_pc(k,idim)))/...
                    (factorial(s-index_pc(i,idim))*factorial(s-index_pc(j,idim))*factorial(s-index_pc(k,idim)));
    else
        cijk_idim = 0;
    end
    cijk = cijk * cijk_idim;
    idim = idim + 1;
end

