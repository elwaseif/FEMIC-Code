% FUNCTION FEMIC_inverse2D.m
%
% This function conmputes the weighted nonlinear least-squares inverse via a
% modified Levenberg-Marquardt scheme with regularized smoothing
% constraints.  Added regularization for smoothing are selected by user to
% produce 2D electrical conducttivity models from frequency-domain EM data.
% The 2D regularization constraint formulation is similar to the that
% developed by Constable et al. for inversion of magnetotellurics data.
%
% The inverse model iteratively call the forward model function
% FEMIC_Jacobian.m until the convergence criteria are met.  The modified LM
% objective function ||d-G(m)||p + muh(R) + muv(R) + (gamma/2)(X) are
% linearized using Taylor-series expansion and lead to:
% A= muv(deltaV'*deltaV)+muh(deltaH'*deltaH)+(gamma/2)*X*'X+(W*J)'*(W*J)
% b=(W*J')*W*(d-G+J*dp)+(gamma/2)*(X'*X)
% where p_trial = A\b (by Cholesky factorization)
% 
% Originated by: Greg Schultz
% Modified from original codes produced in 2004/2005
% Significantly modifield in 2008 to incorporate the log-barrier function
% constraint to enforce positivity and add a number of other features
% Code delivered to the USGS under PO XXXXXXXXX in June 2008
%
%%
% INPUTS:
%   params = model parameters to be optimizaed
%   d = depths (initial) (1 x mlayers)
%   pobs = measurements (expected to be VDM cat HDM both (1 x 2*length(f))
%   sigma = standard deviations on measurements (1 x 2*length(f))
%   f = frequencies
%   r = separation distances between Rx and Tx for bistatic case
%   muh = horizontal regularization coeffecient
%   muv = vertical regularization coeffectient
%   tol_eca = tolerance on changes to conductivity
%   err_tol = convergence criteria for changing errors
%   max_iter = maximum no. of allowable iterations
%   q = [=1,2,3] to designate the data types represented in the pobs input
%       array: 1=Vertical Magnetic Dipole only
%              2=Horizontal Magnetic Dipole only
%              3=Both VMD and HMD data
%%
% OUTPUTS:
%   p_final = the final model array in [P stations x N frequencies (2xN for
%       both VMD and HMD data)] form
%   muh_final = the final horizontal regulatization coeffecient
%   rms_error = the history of the Lp norm rms errors between forward model
%       results (G(m)) and data (d)
%
% EXAMPLE USAGE:
% [g,mmu,of,mdata]=inverse_2d_guiCHC_SIMU(params,d,pobs,sigma,f,r,muh,muv,tol_eca,err_tol,max_iter,q);

function [p_final, muh_final, rms_error2, G]=FEMIC_tost1(MM,S,el,pobs,sigma,muv,err_tol,max_iter,q,sz,wta)

%% Initialize array size and constant parameters
porder = 2;             % norm order (p-value) for error calc (p=2 -> L2 norm)
szp=size(pobs);         % number of frequencies (x2) by number of stations
P=szp(2);               % number of data (freq x 2)
params=MM.con;
d=MM.thk;
f=S.freq;
r=S.r;
M = length(d);%size(params,1);     % number of model parameters
Md = length(d);         % number of layer thicknesses
Ms = M-Md;              % should be the number of conducitivities
N=length(f);            % number of frequencies
if q==3,
    NN=2*N;             % twice the frequenices if both VDM and HDM used
else
    NN=N;               % total number of data points (same as N if only one orientation of data is used)
end
NP=NN*P;                % size of the data set
MP=M*P;  
% size of the model output set
xc=params;
params=repmat(params,1,P);

obs=pobs(:);%reshape(pobs,NP,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Initialize the regularization matrices
W=diag(1./sigma(:));MTX.W=W;
dp=1;
tic;
mref=xc;
%% Outer loop over maximum number of iterations (breaks if alternate
%% convergence criteria are met
itc = 0; 
misfit = []; gc = 1; normg0 = 1;
while(norm(gc)/normg0 > err_tol & itc < max_iter & norm(gc)>1e-20)
   % iteration count
   itc = itc+1;
     for ii=1:P,          % loop over all measurement positions (B1)
       fprintf('          Computing profile %3.0f\n',ii)
       pin(:,ii)=xc(:,ii);  % model input parameters for this position (ii)
       %[GG(:,ii),JJ]=jacnmodelMdl_mat(f,r,pin,d);
       MM.con=pin(:,ii);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Compute the model kernel Gout, and Jacobian matrix J
       [Gout,Jout] =fdem1dfwd(S,MM,el,1); Gout=imag(Gout); Jout=imag(Jout);%[Gout,Jout] = FEMIC_JacobianMdl(f,r,pin,d,N,M,Md,q);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       dlmwrite('Jout.dat',Jout);
       %[t,x,y] = sim('jacobianMdl3');
       
       GG(:,ii) = Gout;%(1:NN);       % place current Gout vector in the global GG array
       JJ = Jout;%(1:NN,1:6);         % organize the global Jacobian matrix
       J((ii-1)*NN+(1:NN),(ii-1)*(M)+(1:M))=JJ(:,1:M); % reorder Jacobian
     end  % end (B1) loop over measurement positions
   % Reshape G into a vector
   G=GG(:);%reshape(GG,NP,1);
   elap_time=toc/60; 
   fd = 0.5*(G-obs)'*(W)*(G-obs);
   if itc == 1  & isempty(muv);
          %%Temporarily assign a value, to be corrected later
          para.BETA = 0;
   end;   
  dx=1*ones(P,1);dz=ones(length(d),1);
  nx=length(dx);nz=length(dz);V = spdiags((wta(:)), 0, nx*nz, nx*nz); 
  dlmwrite('wta.dat',wta);dlmwrite('dx.dat',dx);dlmwrite('dz.dat',dz);
  [Gx, Gz]=grad(dx,dz);sx=1;als=7;%V=0;
   Gs = [sx*Gx;sz*Gz];
   Wt = spdiags(wta(:),0,nx*nz,nx*nz);
   MW = Wt' * ( Gs' * Gs + als * V ) * Wt;pin=pin(:);%reshape(pin,length(pin),1);
   mref=mref(:);%reshape(mref,length(mref),1);
   fm = muv*((xc(:)-mref)'*MW*(xc(:)-mref));
    %%Add them for the total Objective function
   fc =fd+fm;
   grad_fm = muv*MW*(xc-mref);
   grad_fd = J'*(G-obs);dlmwrite('gf.dat',grad_fd);
   % Combine the gradients
   gc = grad_fd + grad_fm;
   %%%% Store some quantities for later use
  
   misfit = sqrt((G-obs)'*W*(G-obs))/sqrt(obs'*W*obs);
   rms_error = misfit;
   if itc == 1, normg0 = norm(gc); f0 = fc; mis0 = misfit; end;
   MTX.mc = xc;
   para.intol = 0;     %  tol for inexact newton solver (ipcg)
   para.inintol = 1e-9;     %  tol for the forward and adjoint problems  
   para.ininintol = 1e-6;   %  tol for the inner solution in the ipcg
   para.init = 4;           %  number of ipcg iterations
   para.ilutol = 0;     %ilu preconditioner tolerance, reduce to 1e-3 if you run into memory issues
   para.alp=1e-4;          % Parameter for line search
   MTX.WTW=MW;dlmwrite('J.dat',J);
   s = ipcg(MTX, muv, -gc, para.intol, para.ininintol, para.init,J,W); 
  
   if max(abs(s)) < 1e-3,
      fprintf('    max_s = %e,  norm(g) = %e\n', max(abs(s)), norm(gc));
      fprintf('STEP size too small CONVERGE  '); % return; 
   end;
pin=pin(:);
   % Try the step 
   mu_LS = 1; 
   iarm = 0;     
   while 1,
      
       
      xt = xc(:) + mu_LS*s;
      
     xt=reshape(xt,length(d),length(xt)/length(d));xt=xt;
     %pin=reshape(pin,length(d),length(xt)/length(d)); 
     for ii=1:P,%% Loop over all station positions (B2)
         %pin=reshape(pin,6,10);
         dlmwrite('pin.dat',pin);dlmwrite('xt.dat',xt);
           pin=[xt(:,ii)];
           MM.con=xt(:,ii); dlmwrite('ccon.dat',xt(:,ii));
           %[G_trial(:,ii),J_trial]=jacnmodelMdl_mat(f,r,pin,d);
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % Compute the model kernel Gout, and Jacobian matrix J
           [Gout,Jout] =fdem1dfwd(S,MM,el,1); Gout=imag(Gout); Jout=imag(Jout);%[Gout,Jout] = FEMIC_JacobianMdl(f,r,pin,d,N,M,Md,q);
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           G_trial(:,ii) = Gout(1:NN);
           J_trial((ii-1)*NN+(1:NN),(ii-1)*(M)+(1:M))=Jout(1:NN,1:M);
      end  % end (B2) loop over station positions
       G_try=G_trial(:);%reshape(G_trial,NP,1);
       
   
      fd = 0.5*(G_try-obs)'*W*(G_try-obs);
      
      %automatically determine a beta guess
      if itc == 1 & muv ==0;
          para.BETA = 0.5*(fd./( (xt-mref)'*muv*(xt-mref))) 
      end;
      fm = 0.5*muv*( (xt(:)-mref)'*muv*(xt(:)-mref));       
      ft = fd+fm;
      fgoal = fc - para.alp*mu_LS*(s'*gc);
      
      if ft < fgoal, 
        break,
      else
   	    % break,
      iarm = iarm+1;
        mu_LS = mu_LS/2;    
      end;  
      
      if(iarm > 5)
       %    disp(' Line search FAIL EXIT(0)');     
        %   return;  
        break,
	  end
     fgoal = fc - para.alp*mu_LS*(s'*gc);
   end  % end line search
   % Update model
   
   % L-Curve
       
   
   xc = xt; 
   misfitnew = misfit;
   misfitold = misfitnew;
   rms_error2(:,itc)=rms_error;
   
  
end  % end loop (A1) over iterations 
p_final=xt;
muh_final=muv;
  
