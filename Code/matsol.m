function[u,iter] = matsol(A,q,tol)
% [u] = matsol(A,q,tol)
% solves the linear system u = A^-1*q
% Uses an SSOR preconditioner and either precg or bicstb.
% Solve using ssor preconditioner
%
% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2006

if norm(q) < 1e-14,
    u = q*0;
    return;
end;

%Generate the preconditioner
% SSOR preconditioner  
Dg = spdiags(diag(A),0,size(A,1),size(A,2));
Dhf = spdiags(1./sqrt(diag(A)),0,size(A,1),size(A,2));
M1 = (Dg + tril(A,-1)) * Dhf;
M2 = Dhf * (Dg + triu(A,1));
 
%Check to see which solver to use (see if the problem is self adjoint)
if norm(A-A','fro') < 1e-14      
     [u, iter] = precg(A, q, M1, M2, 500, tol);
else  
     [u, err1,iter] = bicgstb(A, q, M1, M2, 500, tol);
end