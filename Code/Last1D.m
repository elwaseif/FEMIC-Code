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
% where p_trial = A\b (by Cholesky factorization, Discrepancy, or Max Entropy)
%

% See Refrence: Schultz and Ruppel, 2005, Geophysics
% Also see FEMIC code technical note and manual, 2008
% 
% Originated by: Greg Schultz
% Modified from original codes produced in 2004/2008
% Significantly modifield in 2008 to incorporate the log-barrier function
% constraint to enforce positivity and add a number of other features
% Code delivered to the USGS-Storrs in June 2008
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
%// Inputs are generally specificied by the Matlab GUI code FEMIC_InvGUI.m
%%
% OUTPUTS:
%   p_final = the final model array in [P stations x N frequencies (2xN for
%       both VMD and HMD data)] form
%   muh_final = the final horizontal regulatization coeffecient
%   rms_error = the history of the Lp norm rms errors between forward model
%       results (G(m)) and data (d)
%
% EXAMPLE USAGE:
% [p,mu,errRMS,g]=FEMIC_inverse2D(params,d,pobs,sigma,f,r,muh,muv,tol_eca,err_tol,max_iter,q);

function [p_final,  rms_error, G, best_index]=Last1D(Ma,S,pobs,el,sigma,muh,muv,...
    tol_eca,err_tol,max_iter,q,pmin,pmax,barrier,invType)
LCURVEf = 0;                % Initialize LCURVE method to DEFAULT (=not used)
initmodel_cond=Ma.con;
init_lyrthick=Ma.thk;
f=S.freq;
r=S.r;
params=initmodel_cond;
d=init_lyrthick;
porder = 2;             % norm order (p-value) for error calc (p=2 -> L2 norm)
szp=size(pobs);         % number of frequencies (x2) by number of stations
P=szp(2);               % number of data (freq x 2)
M = length(params);     % number of model parameters
Md = length(d);         % number of layer thicknesses
Ms = M-Md;              % should be the number of conducitivities
N=length(f);            % number of frequencies
if q==3,
    NN=2*N;             % twice the frequenices if both VDM and HDM used
else
    NN=N;               % total number of data points (same as N if only one orientation of data is used)
end
NP=NN*P;                % size of the data set
MP=M*P;                 % size of the model output set
params=(repmat(params,1,P));
obs=pobs(:);%(reshape(pobs,NP,1));
deltav = -diag(ones((M)*P,1))+diag(ones((M)*P-P,1),P);
deltah = -diag(ones((M)*P,1))+diag(ones((M)*P-1,1),1);
sigma=sigma(:);%reshape(sigma, NP, 1);
po=(params(1:M,:));
W=diag(1./sigma);
dp=1;
tic;
switch invType
    case(1),
        inversionTypechar = 'Occams Inversion (Fixed Reg. Coeff.)\n';
    case(2),
        inversionTypechar = 'Truncated SVD (Discrepancy Principle)\n';
    case(3)
        inversionTypechar = 'Maximum Entropy\n';
    case(4),
        inversionTypechar = 'Occams Inversion (L-curve)\n';
    case (5)
         inversionTypechar = 'Biconjugate gradient stabilizing method\n';
    otherwise
        inversionTypechar = 'Unkwown Inversion Method!!!\n';
end
fprintf('******  STARTING FEMIC INVERSION  ******\n');
fprintf(['  Inversion Type: ',inversionTypechar]);
fprintf('  Maximum Iterations: %i; Error Tolerance: %f; Model Change Tolerance: %f\n',max_iter,err_tol,tol_eca);

for i=1:max_iter,       % loop over iterations (A1)
   for ii=1:P,          % loop over all measurement positions (B1)
       fprintf('          Computing profile %3.0f\n',ii)
       pin(:,i)=[po(:,ii)];  % model input parameters for this position (ii)
        Ma.con=((pin(:,i)));
       [Gout,Jout] =fdem1dfwd(S,Ma,el,1); Gout=(imag(Gout));%./real(Gout); 
       Jout=(imag(Jout));%./real(Jout);%dlmwrite('J.dat',Gout);%sqrt((imag(Jout)).^2 + (imag(Jout)).^2);%FEMIC_JacobianMdl(f,r,pin,d,N,M,Md,q);
       GG(:,ii) = Gout(1:NN);       % place current Gout vector in the global GG array
       JJ = Jout(1:NN,1:M);         % organize the global Jacobian matrix
       J((ii-1)*NN+(1:NN),(ii-1)*(M)+(1:M))=JJ(:,1:M); % reorder Jacobian
   end  % end (B1) loop over measurement positions
   G(:,i)=GG(:);
   rms_error(i) = norm(W*G(:,i)-W*obs,porder)/sqrt(length(G(:,i)));
   elap_time=toc/60; 
   fprintf('   Iteration Number: %2.0f rms error: %3.3f Elasped time (mins): %3.2f\n',i,(rms_error(i)),elap_time);
   
   hj=1;vj=1;
   if length(muv)>=length(muh), mu=muv; else mu=muh; end
   for j=1:length(mu),     % Loop over all regularization coeffecients (C1)
       ppo=reshape(po,MP,1);% reorder the the model parameters into a vector
        X = diag( (1./(ppo-pmin)) + (1./(pmax-ppo)),0 );dlmwrite('x.dat',X);
       
       A = muv(j)*deltav'*deltav + muh(j)*deltah'*deltah + (W*J)'*(W*J);%+ barrier.*(X);
       
       b = (W*J)'*W*(obs - G(:,i) + J*ppo);
    switch invType
           case(1),
               p_trial(:,j)=(A\b);
           case(2),
               [U,s,V] = FEMIC_svd(A);
               [p_trial(:,j),mul]=FEMIC_discrep(U,s,V,b,err_tol/i,ppo);
               mu = mul;
           case(3),
               [p_trial(:,j)] = FEMIC_maxentropy(A,b,mu(j),1,ppo);
           case(4),
               p_trial(:,j)=(A\b);
               LCURVEf=1;
           case (5)
               p = symrcm(A);
               [iii,up] = sort(p);
               [PL,PU] = luinc(sparse(A), 1e-4);
               p_trial(:,j) = bicgstb(A, b, PL, PU, 500, 1e-9);
           otherwise
               fprintf('ERROR: No proper inversion method selected\n');
               return
       end
        lbf = 2*barrier*sum( (log(ppo-pmin)+log(pmax-ppo)) );
       ptemp=(reshape(p_trial(:,j),M,P));  % vector format for temporary trial model parameters
       dlmwrite('ptemp.dat',ptemp);
       
       for ii=1:P,%% Loop over all station positions (B2)
           pin2=(ptemp(:,ii));
          Ma.con=((pin2));
           [Gout,Jout] = fdem1dfwd(S,Ma,el,1); Gout=(imag(Gout));%./real(Gout);
           Jout=(imag(Jout));%./real(Jout);%sqrt((imag(Jout)).^2 + (imag(Jout)).^2); %FEMIC_JacobianMdl(f,r,pin,d,N,M,Md,q);
           
           G_trial(:,ii) = Gout(1:NN);
           J_trial((ii-1)*NN+(1:NN),(ii-1)*(M)+(1:M))=Jout(1:NN,1:M);
       end  % end (B2) loop over station positions
       G_try(:,j)=reshape(G_trial,NP,1);
       rms_trial(j) = norm(W*G_try(:,j)-W*obs, porder)/sqrt(length(G_try(:,j)));
       if LCURVEf & length(mu)>3,  % needs at least 4 points to find L-curve
           residual_norm(j) = rms_trial(j)*sqrt(length(G_try(:,j)));
           solution_norm(j) = mu(j);
       end
       if length(muv)>=length(muh), vj=vj+1; else hj=hj+1; end
   end  % end loop (C1) over all regularization coeffecients (C1)
   barrier = barrier./5;%((rms_error(i)-rms_trial(j))/rms_error(i)) )
   % Evaluate error tolerance in rms errors
   aq1=sum(G);aq2=sum(G_try);aq3=find(aq1<-0.01);aq4=find(aq2<-0.01);  
   rms_error(isnan(rms_error))=100000000;
   rms_trial(isnan(rms_trial))=100000000;
   rms_error(aq3)=100000000;
   rms_trial(aq4)=100000000;

   if min(rms_trial) > err_tol,
       [best_rms, best_index] = min(rms_trial);

   else

       index = find(rms_trial <= err_tol);
       if isempty(index)
            %rms_trial(isnan(rms_trial))=100000000000;
           [best_rms, best_index] = min(rms_trial);
           p_final = p_trial(:,best_index);
           rms_error(i+1)  = best_rms;
           G(:,i+1) = G_try(:,best_index);
           mu_final = mu(best_index);
           return;
     
       end
   end
   pnew = p_trial(:,best_index);
  rms_vs_chg = norm((pnew' - po')./po')/sqrt(length(po'));
   if i==1
   if rms_vs_chg < tol_eca | (best_rms <= err_tol) | (i==max_iter)
     [best_rms, best_index] = min(rms_error);
       p_final =(pin(:,best_index));
      rms_error(i+1)  = best_rms;
      G(:,i+1) = G(:,best_index);
      fprintf('!!*******      Model Sequence Completed Successfully     **********!!\n')
      return;
   end
   else
       chg=(abs(rms_error(i-1)-rms_error(i)))/rms_error(i-1);
       if rms_vs_chg < tol_eca | (best_rms <= err_tol) | (i==max_iter) | chg<0.1 %| rms_error(i)>rms_error(i-1)
     [best_rms, best_index] = min(rms_error);
       p_final = (pin(:,best_index));
      rms_error(i+1)  = best_rms;
      G(:,i+1) = G(:,best_index);
      fprintf('!!*******      Model Sequence Completed Successfully     **********!!\n')
     
      return;
       end
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  L-curve Determination of Regularization coeffecient
    if LCURVEf, 
        if length(mu)>3,
            % Use cubic spline interpolant to creat L-curve at 351 points
            [nresid_norm, nsoln_norm] = cubicSpline([residual_norm' solution_norm'],32);
            % Use Hansen et al., 2007 pruning algorithm to find the "corner" of
            % the L-curve
            [k_corner,info] = FEMIC_corner(nresid_norm, nsoln_norm);
            new_mu = nsoln_norm(k_corner);
            if isempty(new_mu), new_mu=mu; end
           mu = [new_mu/20;new_mu/8;new_mu/4;4*new_mu;8*new_mu;20*new_mu];
            if length(muv)>length(muh), muv=mu; else muh=mu; end
        end
    end


   po=((pnew)); %Ma.con=(ptemp);
end  % end loop (A1) over iterations 
  
