function [x,w] = mherzo(n)
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ===========================================================
%       Purpose : This program computes the zeros of Hermite
%                 polynomial Ln(x)in the interval[-ì,ì]and the
%                 corresponding weighting coefficients for Gauss-
%                 Hermite integration using subroutine HERZO
%       Input :   n    --- Order of the Hermite polynomial
%                 X(n)--- Zeros of the Hermite polynomial
%                 W(n)--- Corresponding weighting coefficients
%       ===========================================================
% n =[];
x=[];w=[];
x=zeros(1,100);
w=zeros(1,100);
[n,x,w]=herzo(n,x,w);
x = x(1:n);
w = w(1:n);

function [n,x,w]=herzo(n,x,w,varargin);
%       ========================================================
%       Purpose : Compute the zeros of Hermite polynomial Ln(x)
%                 in the interval[-ì,ì], and the corresponding
%                 weighting coefficients for Gauss-Hermite
%                 integration
%       Input :   n    --- Order of the Hermite polynomial
%                 X(n)--- Zeros of the Hermite polynomial
%                 W(n)--- Corresponding weighting coefficients
%       ========================================================
hn=1.0d0./fix(n);
zl=-1.1611d0+1.46d0.*fix(n).^0.5;
for  nr=1:fix(n./2);
if(nr == 1)z=zl; end;
if(nr ~= 1)z=z-hn.*(fix(n./2)+1-nr); end;
it=0;
while (1);
it=it+1;
z0=z;
f0=1.0d0;
f1=2.0d0.*z;
for  k=2:n;
hf=2.0d0.*z.*f1-2.0d0.*(k-1.0d0).*f0;
hd=2.0d0.*k.*f1;
f0=f1;
f1=hf;
end;  k=fix(n)+1;
p=1.0d0;
for  i=1:nr-1;
p=p.*(z-x(i));
end;  i=nr-1+1;
fd=hf./p;
q=0.0d0;
for  i=1:nr-1;
wp=1.0d0;
for  j=1:nr-1;
if(~(j == i))wp=wp.*(z-x(j)); end;
end;  j=nr-1+1;
q=q+wp;
end;  i=nr-1+1;
gd=(hd-q.*fd)./p;
z=z-fd./gd;
if(~(it <= 40&abs((z-z0)./z)> 1.0d-15))break; end;
end;
x(nr)=z;
x(n+1-nr)=-z;
r=1.0d0;
for  k=1:n;
r=2.0d0.*r.*k;
end;  k=fix(n)+1;
w(nr)=3.544907701811d0.*r./(hd.*hd);
w(n+1-nr)=w(nr);
end;
if(n ~= 2.*fix(n./2));
r1=1.0d0;
r2=1.0d0;
for  j=1:n;
r1=2.0d0.*r1.*j;
if(j >= fix((n+1)./2))r2=r2.*j; end;
end;  j=fix(n)+1;
w(fix(n./2)+1)=0.88622692545276d0.*r1./(r2.*r2);
x(fix(n./2)+1)=0.0d0;
end;

% Next normalize to change to probabilistic weights (Modiflied by AD)
x = x * sqrt(2);
w = w / sqrt(pi);

return;
