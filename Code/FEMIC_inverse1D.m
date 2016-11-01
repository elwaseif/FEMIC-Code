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
% priori = user defined priori information 
% sx=smoothness in the x-direction
%sz= smoothness in the z-direction
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

function [p_final, mu_final, rms_error, G]=FEMIC_inverse1D(initmodel_cond, init_lyrthick,pobs,sigma,coords,f,r,muh,muv,...
    tol_eca,err_tol,max_iter,q,pmin,pmax,barrier,invType,priori,sz, plotsType, statusUpdate,el)
%% Pre-initialization Flags
LCURVEf = 0;                % Initialize LCURVE method to DEFAULT (=not used)
%% Initialize array size and constant parameters
wprof = 1;   % selection of the station number for interim plots
boo=init_lyrthick;
sigma=sigma(:);
for i=1:length(initmodel_cond)
    aaaa(i)=length(initmodel_cond{i});
    boo{i}=zeros(length(boo{i}),1);
end
for i=1:length(initmodel_cond)
    
    if length(initmodel_cond{i})<max(aaaa);
        dd=initmodel_cond{i};
        da=init_lyrthick{i};
        fa=boo{i};
        dd(length(initmodel_cond{i}):max(aaaa))=dd(length(initmodel_cond{i}))*(ones(length(length(initmodel_cond{i}):max(aaaa)),1));
        da(length(init_lyrthick{i}):max(aaaa))=da(length(init_lyrthick{i}))*(ones(length(length(init_lyrthick{i}):max(aaaa)),1));
        fa(length(boo{i}):max(aaaa))=zeros(length(length(boo{i}):max(aaaa)),1);
    else
       dd=initmodel_cond{i};
       da=init_lyrthick{i};
       fa=boo{i};
    end
    params(:,i)=dd;
    d(:,i)=da;
    wta(:,i)=fa;
end

porder = 2;             % norm order (p-value) for error calc (p=2 -> L2 norm)
P=size(pobs,2);         % number of frequencies (x2) by number of stations
PP=size(pobs,1);               % number of data (freq x 2)
%M=length(d);     
M = size(params,1);     % number of model parameters
Md = size(d,1);         % number of layer thicknesses
Ms = M-Md;              % should be the number of conducitivities
%N=szp(1);
N=length(f);            % number of frequencies
if q==3,
    NN=2*N;             % twice the frequenices if both VDM and HDM used
else
    NN=N;               % total number of data points (same as N if only one orientation of data is used)
end
NP=NN*P;                % size of the data set
MP=M*P;                 % size of the model output set
wta=zeros(Md,P);
dx=ones(P,1);dz=ones(Md,1);nx=length(dx);nz=length(dz);
% import the user defined priori information in the starting model
[faa fo]=size(priori);

        for lo=1:faa
        ll=priori(lo,1);
        qq=priori(lo,2);
        params(qq,ll)=priori(lo,3);
        wta(qq,ll)=priori(lo,4);
        end
        obs=pobs(:);%(reshape(pobs,NP,1));
        [Gx Gz]=grad(dx,dz);sx=1;
        Gs = [sx*Gx;sz*Gz];
        wta=1-wta;
        Wt = spdiags(wta(:),0,nx*nz,nx*nz);
        MW = Wt' * ( Gs' * Gs) * Wt;
       %sigma=repmat(sigma, P, 1);
       po=params(1:M,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

W=diag(1./sigma);
dp=1;
tic;
MM.con=po;
MM.k=length(po);
MM.thk=d;
MM.chie=zeros(length(d),1);
MM.chim=zeros(length(d),1);
S.freq=f;
for i=1:length(f);S.tor{i}='z';end;S.tor=reshape(S.tor,length(f),1);
S.tmom=ones(length(f),1);
S.tx=zeros(length(f),1);
S.ty=zeros(length(f),1);
S.tzoff=zeros(length(f),1);
S.ror=S.tor;
S.rmom=ones(length(f),1);
S.rx=1.66*ones(length(f),1);
S.ry=zeros(length(f),1);
S.rzoff=zeros(length(f),1);
S.nf=length(f);
S.r=S.rx;
%el=-1;
deltav = -diag(ones((M)*P,1))+diag(ones((M)*P-P,1),P);
deltah = -diag(ones((M)*P,1))+diag(ones((M)*P-1,1),1);
% case 5 added by Mehrez
switch invType
    case(1),
        inversionTypechar = 'Occams Inversion (Fixed Reg. Coeff.)';
    case(2),
        inversionTypechar = 'Truncated SVD (Discrepancy Principle)';
    case(3)
        inversionTypechar = 'Maximum Entropy';
    case(4),
        inversionTypechar = 'Occams Inversion (L-curve)';
    case (5)
         inversionTypechar = 'Biconjugate gradient stabilizing method';
    case (6)
         inversionTypechar = 'Thiknov inversion';
    otherwise
        inversionTypechar = 'Unknown Inversion Method!!!';
end

fprintf('******  STARTING FEMIC INVERSION  ******\n');

str='******  STARTING FEMIC INVERSION  ******';
statusUpdate(str);

fprintf(['  Inversion Type: ',inversionTypechar]);
fprintf('\n');

str=strcat('Inversion Type: ', inversionTypechar);
statusUpdate(str);

fprintf('  Maximum Iterations: %i; Error Tolerance: %f; Model Change Tolerance: %f\n',max_iter,err_tol,tol_eca);
str=strcat('Maximum Iterations: ', num2str(max_iter), ' Error tolerance: ', num2str(err_tol), ' Model Change tolerance: ', num2str(tol_eca));
statusUpdate(str);
%% Outer loop over maximum number of iterations (breaks if alternate
%  convergence criteria are met)
for i=1:max_iter,       % loop over iterations (A1)
   
   if length(muv)>length(muh), mu=muv; else mu=muh; end
   totLoops=length(mu)+P;
   
   str=strcat('ITERATION: ', num2str(i));
   progressbar(str)
   
   str='Looping through measurement positions';
   statusUpdate(str);
    
    for ii=1:P,          % loop over all measurement positions (B1)
       fprintf('          Computing profile %3.0f\n',ii)
              
       status_perc=ii/totLoops;
       progressbar(status_perc);
       
       str=strcat('Computing profile: ', num2str(ii));
       statusUpdate(str);
        
       
        pin(:,i)=po(:,ii);  % model input parameters for this position (ii)
        MM.con=((pin(:,i)));
        [Gout,Jout] =fdem1dfwd(S,MM,el,1); Gout=(imag(Gout));%./real(Gout); 
        Jout=(imag(Jout));
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       GG(:,ii) = Gout(1:NN);       % place current Gout vector in the global GG array
       JJ = Jout(1:NN,1:M);         % organize the global Jacobian matrix
       J((ii-1)*NN+(1:NN),(ii-1)*(M)+(1:M))=JJ(:,1:M); % reorder Jacobian
   end  % end (B1) loop over measurement positions
   G(:,i)=GG(:);
   % Compute the Lp (p=porder) rms error
   rms_error(i) = norm(W*G(:,i)-W*obs,porder)/sqrt(length(G(:,i)));
   elap_time=toc/60; 
   
   fprintf('   Iteration Number: %2.0f rms error: %3.3f Elasped time (mins): %3.2f\n',i,(rms_error(i)),elap_time);

   str=strcat('Iteration Number: ', num2str(i), '  rms error: ', num2str(rms_error(i)), ' Elasped time (mins): ', num2str(elap_time));
   statusUpdate(str);

   hj=1;vj=1;
   if length(muv)>length(muh), mu=muv; else mu=muh; end
   for j=1:length(mu),     % Loop over all regularization coeffecients (C1)
       status_perc=(P+j)/totLoops;
       progressbar(status_perc);
       
       ppo=po(:);%reshape(po,MP,1);% reorder the the model parameters into a vector
       
        X = diag( (1./(ppo-pmin)) + (1./(pmax-ppo)),0 );
       
       % Compute the RHS matrix system for inversion
       A = muv(vj)*MW'*MW + (W*J)'*(W*J);% + barrier.*X;
      
       switch invType
           case(2),
               p_trial(:,j)=(A\b);
           case(3),
               [U,s,V] = FEMIC_svd(A);
               [p_trial(:,j),mul]=FEMIC_discrep(U,s,V,b,err_tol/i,ppo);
               mu = mul;
           case(4),
               [p_trial(:,j)] = FEMIC_maxentropy(A,b,mu(j),1,ppo);
           case(5),
               p_trial(:,j)=(A\b);
               LCURVEf=1;
           case (6)
               p = symrcm(A);
               [iii,up] = sort(p);
               setup.type = 'crout';
               setup.milu = 'col';
               setup.droptol = 1e-4;
               [PL,PU] = ilu(sparse(A), setup);
               p_trial(:,j) = bicgstb(A, b, PL, PU, 500, 1e-9);
           case (1)           
               muv=muv(vj);
              [p_final, mu_final, rms_error, G]=FEMIC_tost1(MM,S,el,pobs,sigma,muv,err_tol,max_iter,q,sz,wta);
              [hplot] = FEMIC_plotINVinterim(G(N*(wprof-1)+1:N*wprof,1),obs(N*(wprof-1)+1:N*wprof),...
            f,i,max_iter,p_final(:,wprof),d,rms_error(length(rms_error)),params,pmin,pmax);
              return;
           otherwise
               fprintf('ERROR: No proper inversion method selected\n');
               return
       end
       ptemp=(reshape(p_trial(:,j),M,P));     
       ccc=zeros(max(aaaa),length(init_lyrthick));
       for cnt=1:length(init_lyrthick)
            if length(init_lyrthick{cnt})>=max(aaaa)
                bbb=init_lyrthick{cnt};
            end
            ccc(:,cnt)=bbb;
       end
       [m n]=size(ccc);
       ccc(m,:)=ccc(m-1,:);
       ea=cumsum(ccc);zz=max(ea')';
       x=coords(:,2);
       [X,Y]=meshgrid(x,zz);
     
       lbf = 2*barrier*sum( (log(ppo-pmin)+log(pmax-ppo)) );
        % vector format for temporary trial model parameters
       
       for ii=1:P,%% Loop over all station positions (B2)
           pin2=ptemp(:,ii);
           MM.con=((pin2));
           [Gout,Jout] = fdem1dfwd(S,MM,el,1); Gout=(imag(Gout));%./real(Gout);
           Jout=(imag(Jout));
          G_trial(:,ii) = Gout(1:NN);
          J_trial((ii-1)*NN+(1:NN),(ii-1)*(M)+(1:M))=Jout(1:NN,1:M);
       end  % end (B2) loop over station positions
       G_try(:,j)=G_trial(:);%reshape(G_trial,NP,1);
       rms_trial(j) = norm(W*G_try(:,j)-W*obs, porder)/sqrt(length(G_try(:,j)));
       %keyboard
       if LCURVEf & length(mu)>3,  % needs at least 4 points to find L-curve
           residual_norm(j) = rms_trial(j)*sqrt(length(G_try(:,j)));
           solution_norm(j) = mu(j);
       end
       if length(muv)>length(muh), vj=vj+1; else hj=hj+1; end
   end  % end loop (C1) over all regularization coeffecients (C1)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
  
       
      
%% Compute convergence criteria

   str=strcat('Computing convergence criteria');
   statusUpdate(str);
    % reduce barrier function by order of magnitude/2 = 5
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
           mu_final = muv;% mu(best_index);
           return;
       end
   end
   pnew = p_trial(:,best_index);
  
    %check for convergence
   rms_vs_chg = norm((pnew' - po')./po')/sqrt(length(po'));
   if i==1
   if rms_vs_chg < tol_eca | (best_rms <= err_tol) | (i==max_iter)
       [best_rms, best_index] = min(rms_error);
       pnew= pin(:,best_index);
       baa= shiftdim(pnew(find(wta)));
       da = params;
       da(find(wta)) = baa;%pnew;
       pnew = da;
       p_final = pnew;

       rms_error(i+1)  = best_rms;
       G(:,i+1) = G(:,best_index);
       mu_final = muv;% mu(best_index);
       fprintf('!!*******      Model Sequence Completed Successfully     **********!!\n')
      
       
      str=strcat('!!*******      Model Sequence Completed Successfully     **********!!');
      statusUpdate(str);
      progressbar(1);
      %if plotsType(1),
        [hplot] = FEMIC_plotINVinterim(G_try(N*(wprof-1)+1:N*wprof,best_index),obs(N*(wprof-1)+1:N*wprof),...
            f,i,max_iter,p_final(:,wprof),d,best_rms,params,pmin,pmax);
      %end
      return;
   end
   else
        chg=(abs(rms_error(i-1)-rms_error(i)))/rms_error(i-1);
       if rms_vs_chg < tol_eca | (best_rms <= err_tol) | (i==max_iter) | chg<0.1 %| rms_error(i)>rms_error(i-1)
        %  rms_error(isnan(rms_error))=10000000000;
     [best_rms, best_index] = min(rms_error);
       p_final = (pin(:,best_index));
      rms_error(i+1)  = best_rms;
      G(:,i+1) = G(:,best_index);
      mu_final = muv;%mu(best_index);
      fprintf('!!*******      Model Sequence Completed Successfully     **********!!\n')
      [hplot] = FEMIC_plotINVinterim(G_try(N*(wprof-1)+1:N*wprof,best_index),obs(N*(wprof-1)+1:N*wprof),...
            f,i,max_iter,p_final(:,wprof),d,best_rms,params,pmin,pmax);
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
            
            if plotsType(2),
                %figure(11);
                set(gcf,'Position',[1000 500 600 600]);
                plot(log10(nresid_norm),log10(1./nsoln_norm) ,'k.-','LineWidth',i);
                hold on; plot(log10(residual_norm),log10(1./solution_norm),'go','MarkerSize',12);
                hold on; plot(log10(nresid_norm(k_corner)),log10(1./nsoln_norm(k_corner)),'r*','MarkerSize',14);
                hold off;
                xlabel('log[residual norm]'); ylabel('log[solution seminorm]');
                set(gca,'FontSize',12)
                shg
            end
          new_mu = nsoln_norm(k_corner);
            if isempty(new_mu), new_mu=mu; end
            mu = [new_mu/20;new_mu/8;new_mu/4;4*new_mu;8*new_mu;20*new_mu];
            if length(muv)>length(muh), muv=mu; else muh=mu; end
        end
    end

%%   Plotting of intermediate results
    %if plotsType(1),
        [hplot] = FEMIC_plotINVinterim(G_try(N*(wprof-1)+1:N*wprof,best_index),obs(N*(wprof-1)+1:N*wprof),...
            f,i,max_iter,pnew(:,wprof),d,best_rms,params,pmin,pmax);
    %end
%%
   po=(pnew);
end  % end loop (A1) over iterations 
  
