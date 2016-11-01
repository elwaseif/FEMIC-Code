function [x,iter] = precg(A, b, M1, M2, max_it, tol)
% [x,numiter] = precg(A, b, M1, M2, max_it, tol)
% 
% Preconditioned CG solver
%

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2006
  

%initialize the starting variables
  r = b;
  n = length(b);
  x = zeros(n,1);
  rho = r'*r;
  
  if  ( rho == 0.0 ), rho = 1.0; end

  err = norm( r ) / sqrt(rho);
  if ( err < tol ) disp('err < tol'); return, end

  for iter = 1:max_it                       % begin iteration

     % preconditioning step %%%     
     z = M2\(M1\r);  
     
     rho_1 = rho;
     rho = (r'*z);
     if iter == 1, rho0 = norm(z); end;

     if ( iter > 1 ),                       % Calculate the direction vector
        beta = rho / rho_1;
        p = z + beta*p;
     else
        p = z;
     end
     %%%%%%  Matrix times a vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % q = A*p
     q = A*p;
 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     alpha = rho / (p'*q );
     x = x + alpha * p;                    % update approximation vector

     r = r - alpha*q;                      % compute residual
     err = norm( r ) / rho0;               % check convergence
     numiter = iter;
     %fprintf('        PRECG iteration %d, Relative residual = %e\n',iter, err);
     if ( err <= tol )
        break, 
     end 
     res(iter) = err;
     
  end 
  %END precg.m

