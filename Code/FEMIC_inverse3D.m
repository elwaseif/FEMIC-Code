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

function [p_final, mu_final, rms_error, G, x, zz]=FEMIC_inverse3D(initmodel_cond, init_lyrthick,pobs,sigma,f,r,coords,muh,muv,...
    tol_eca,err_tol,max_iter,q,pmin,pmax,barrier,invType,priori,sx,sy,sz, statusUpdate,sens,vall,plotdoi,perc,el)

%% Pre-initialization Flags
%XXXxGSVDf = 0;             % Flag for Generalized Singular Value Decomposition method for regularization (never used!!- gms)
%XXXxCONDf = 0;             % Flag to determin condition number  (never %used!!- gms)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotTypes: Plot types 
% NONE = [0 0];
% PLOTINTER = [1 0];        % Flag for plotting simple intermediate results
% LCURVEf = [1 1];          % Flag for automatic calcuation of the regularization parameter by max curvature of L-curve
LCURVEf = 0;                % Initialize LCURVE method to DEFAULT (=not used)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% invTypes: Inversion Method types (Examples)
% 1 = CHOLESKYf;            % Flag for Regular Cholesky facotrization
% 2 = DISCREPf;             % Flag for Discrepancy Principle (effecient Newton step; set reg. coeff.
% 3 = MAXENTROPYf;          % Flag for Maximum Entropy Inversion (nonlinear conjugate gradient algorithm with inexact
%                             % line search and a step-length control which ensures that all elements of the iteration
%                             % vector are positive.)
%% Initialize array size and constant parameters

parpool('local',vall)

boo=init_lyrthick;
sigma=sigma(:);

if LCURVEf, 
    muv=muv;
else
    muv=muh;
end
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
PPP=size(pobs,1);
P=length(unique(coords(:,2)));
PP=length(unique(coords(:,3)));
pobs=reshape(pobs,PPP,P*PP);
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
NP=NN*P*PP;                % size of the data set
MP=M*P*PP;                 % size of the model output set
obs=(reshape(pobs,NP,1));
dx=ones(P,1);dz=ones(M,1);dy=ones(PP,1);
als=1;
dx=1*ones(PPP,1);dz=ones(PPP,1);dy=ones(PP,1);
  nx=length(dx);nz=length(dz);ny=length(dy);
  [Gx, Gy,Gz]=grad3D(dx,dy,dz);als=1;%V=0;
   Gs = [sx*Gx;sy*Gy;sz*Gz];
   Wt = spdiags(wta(:),0,nx*ny*nz,nx*ny*nz);
   MW = Wt' * ( Gs' * Gs) * Wt;
po=params;%(1:M,:);
% import the user defined priori information in the starting model
[faa fo]=size(priori);
        for lo=1:faa
        ll=priori(lo,1);
        qq=priori(lo,2);
        params(qq,ll)=priori(lo,3);
        wta(qq,ll)=priori(lo,4);
        end
        obs=(reshape(pobs,NP,1));
        po=params(1:M,:);

W=diag(1./sigma(:));
dp=1;
tic;
MM.k=length(M);
MM.thk=d;
MM.chie=zeros(MM.k,1);
MM.chim=zeros(MM.k,1);
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

switch invType
    case(1),
        inversionTypechar = 'Occams Inversion (Fixed Reg. Coeff.)';
    case(1),
        inversionTypechar = 'Truncated SVD (Discrepancy Principle)';
    case(1)
        inversionTypechar = 'Maximum Entropy';
    case(1),
        inversionTypechar = 'Occams Inversion (L-curve)';
    case (1)
         inversionTypechar = 'Biconjugate gradient stabilizing method';
    case (1)
         inversionTypechar = 'Thiknov inversion';
    otherwise
        inversionTypechar = 'Unkwown Inversion Method!!!';
end
fprintf('******  STARTING FEMIC INVERSION  ******\n');
str='******  STARTING FEMIC INVERSION  ******';
statusUpdate(str);

fprintf(['  Inversion Type: ',inversionTypechar]);
fprintf('\n');
str=strcat('Inversion Type: ', inversionTypechar);
statusUpdate(str)

fprintf('  Maximum Iterations: %i; Error Tolerance: %f; Model Change Tolerance: %f\n',max_iter,err_tol,tol_eca);
str=strcat('Maximum Iterations: ', num2str(max_iter), ' Error tolerance: ', num2str(err_tol), ' Model Change tolerance: ', num2str(tol_eca));
statusUpdate(str);

first=true;
%% Outer loop over maximum number of iterations (breaks if alternate
%  convergence criteria are met)
for i=1:max_iter,       % loop over iterations (A1)

   if plotdoi
       if i<2
       [doi]=FEMIC_DOI(d(:,1),pobs,sigma(1:length(f)),f,r,muv(1),err_tol,max_iter,q,vall,perc);
       end
   end
    if length(muv)>length(muh), mu=muv; else mu=muh; end
    totLoops=(length(mu)*P)+P + length(mu);
    
    str=strcat('ITERATION: ', num2str(i));
    progressbar(str)
   
    str='Looping through measurement positions';
    statusUpdate(str);
    
    for ii=1:P*PP,          % loop over all measurement positions (B1)
       fprintf('          Computing profile %3.0f\n',ii)
       
       status_perc=ii/totLoops;
       progressbar(status_perc);
       
       str=strcat('Computing profile: ', num2str(ii));
       statusUpdate(str);
       pin(:,:,i)=po;  % model input parameters for this position (ii)
       pin2=po(:,ii);
       %[GG(:,ii),JJ]=jacnmodelMdl_mat(f,r,pin,d);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Compute the model kernel Gout, and Jacobian matrix J
       MM.con=pin2;MM.thk=init_lyrthick{ii};
       %save MM;
       [Gout,Jout] =fdem1dfwd(S,MM,el,1); Gout=(imag(Gout)); Jout=(imag(Jout));
       GG(:,ii) = Gout(1:NN);       % place current Gout vector in the global GG array
       JJ = Jout(1:NN,1:M);         % organize the global Jacobian matrix
       J((ii-1)*NN+(1:NN),(ii-1)*(M)+(1:M))=JJ(:,1:M);
   end  % end (B1) loop over measurement positions
  
   %keyboard;
   % Reshape G into a vector
   G(:,i)=GG(:);
   % Compute the Lp (p=porder) rms error
   rms_error(i) = norm(W*G(:,i)-W*obs,porder)/sqrt(length(G(:,i)));
   elap_time=toc/60; 
   fprintf('   Iteration Number: %2.0f rms error: %3.3f Elasped time (mins): %3.2f\n',i,(rms_error(i)),elap_time);
   
   str=strcat('Iteration Number: ', num2str(i), '  rms error: ', num2str(rms_error(i)), ' Elasped time (mins): ', num2str(elap_time));
   statusUpdate(str);
   
   deltav = -diag(ones((M)*P,1))+diag(ones((M)*P-P,1),P);

   hj=1;vj=1;
   if length(muv)>length(muh), mu=muv; else mu=muh; end
   for j=1:length(mu),     % Loop over all regularization coeffecients (C1)

       status_perc=(P+j+(P*(j-1)))/totLoops;
       progressbar(status_perc);

       ppo=reshape(po,MP,1);%reshape(po,MP,1);% reorder the the model parameters into a vector
       
       % Compute the log-barrier function matrix for the inequality
       % constraint
       X = diag( (1./(ppo-pmin)) + (1./(pmax-ppo)),0 );
       
       % Compute the RHS matrix system for inversion
       %A = muv(vj)*MW'*MW + (W*J)'*(W*J) + barrier.*X;
       A =  muv(j)*MW' + (W*J)'*(W*J);% + barrier.*X;

      % A = WTW'*WTW+ (W*J)'*(W*J) + barrier.*X;

       % Compute the LHS matrix system
       
       b = (W*J)'*W*(obs - G(:,i) + J*ppo);
%        % Compute the trial model parameters by Cholesky facotrization
%        p_trial(:,j)=(A\b);LCURVEf=1;
ccc=zeros(max(aaaa),length(init_lyrthick));
       for cnt=1:length(init_lyrthick)
           if length(init_lyrthick{cnt})>=max(aaaa)
             bbb=init_lyrthick{cnt};
             ccc(:,cnt)=bbb;
           end
       end
       [m n]=size(ccc);
       ccc(m,:)=ccc(m-1,:);
       ea=cumsum(ccc);zz=max(ea')';
       %x=coords(:,2);
       %[X,Y]=meshgrid(x,zz);
       %keyboard
       switch invType
           case(2),
               p_trial(:,j)=(A\b);
           case(3),
               [U,s,V] = FEMIC_svd(A);
               [p_trial(:,j),mul]=FEMIC_discrep(U,s,V,b,err_tol/i,ppo);
               %max_iter=3;
               mu = mul;
           case(4),
               [p_trial(:,j)] = FEMIC_maxentropy(A,b,mu(j),1,ppo);
               %max_iter=1;
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
                dlmwrite('params.dat',params);
               muv=max(muv);MM.con=params;
               save MM
               pobs=reshape(pobs,PPP,P,PP);
               [p_final, muh_final, rms_error, G]=FEMIC_tost3(MM,S,el,pobs,sigma,muv,err_tol,max_iter,q,sx,sy,sz,wta,pmin,pmax)
               mu_final=muv(1);rms_error=min(rms_error);
               [X Y Z]=meshgrid(1:length(dx),1:length(dy),cumsum(MM.thk));
               xslice=[1 7];yslice=[7];zslice=[6];
               figure;
                slice(X,Y,Z,reshape(p_final,PPP,P,PP),xslice,yslice,zslice)
               return;
           otherwise
               fprintf('ERROR: No proper inversion method selected\n');
               return
       end
       
       
              
       % A priori barrier term
       lbf = 2*barrier*sum( (log(ppo-pmin)+log(pmax-ppo)) );
       ptemp=(reshape(p_trial(:,j),M,P*PP));  % vector format for temporary trial model parameters
       
       for ii=1:P*PP,%% Loop over all station positions (B2)
           
           
           status_perc=(P+j+(P*(j-1))+ii)/totLoops;
           progressbar(status_perc);
           
           pin2=ptemp(:,ii);
           %[G_trial(:,ii),J_trial]=jacnmodelMdl_mat(f,r,pin,d);
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           MM.con=(pin2);MM.thk=init_lyrthick{ii};
           [Gout,Jout] = fdem1dfwd(S,MM,el,1); Gout=(imag(Gout)); Jout=(imag(Jout));
           % Compute the model kernel Gout, and Jacobian matrix J
           %[Gout,Jout(:,:,ii)] = FEMIC_JacobianMdl(f,r,pin,d(:,ii),N,M,Md,q);
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %[G,J] = jacobianMdl_mat_opt_depths(freq,rspa,params,d,N,M,nlayers)
           %[t,x,y] = sim('jacobianMdl3');
           G_trial(:,ii) = Gout(1:NN);%Jout(1:NN,1:M,ii)=Jout;
           J_trial((ii-1)*NN+(1:NN),(ii-1)*(M)+(1:M))=Jout(1:NN,1:M);
          
       end  % end (B2) loop over station positions
       
       G_try(:,j)=reshape(G_trial,NP,1);
       rms_trial(j) = norm(W*G_try(:,j)-W*obs, porder)/sqrt(length(G_try(:,j)));
       %keyboard
       if LCURVEf & length(mu)>3,  % needs at least 4 points to find L-curve
           residual_norm(j) = rms_trial(j)*sqrt(length(G_try(:,j)));
           solution_norm(j) = mu(j);
       end
       if length(muv)>length(muh), vj=vj+1; else hj=hj+1; end
   end  % end loop (C1) over all regularization coeffecients (C1)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [m n]=size(d);
   d(m,:)=d(m-1,:);
   ea=cumsum(d);zz=max(ea')';
   for iii=0:length(coords(:,1))+1
      xy{iii+1} = find(coords(:,2) == iii);
      yx{iii+1} = find(coords(:,3) == iii);
      x = xy(~cellfun('isempty',xy));
      y = yx(~cellfun('isempty',yx));
   end
   x=cell2mat(x);x=coords(x(1,:),2);
   y=cell2mat(y);y=coords(y(1,:),3);
   dlmwrite('x.txt',x);
   dlmwrite('y.txt',y);
   [X Y Z]=meshgrid(x,y,zz);
   dlmwrite('zz.txt',zz);
   xslice = max(y); yslice = max(y); zslice = [min(zz),zz(floor(length(zz)/3)),zz(floor(2*length(zz)/3)),max(zz)]; 
   
   str=strcat('FEMIC 3D INVERSION   Iteration: ', num2str(i));
   if first
        figure('Name', str);
        first=false;
   else
        close();
        figure('Name', str);
   end
   sense=diag(J'*(W'*W)*J);sense=reshape(sense,M,PP*P);sense=sense';
   sense=reshape(sense,length(y),length(x), M);
   p_trail=reshape(p_trial(:,j),length(y),length(x), M);
   dlmwrite('ptrial.txt',p_trail)
   if sens || plotdoi
       subplot(1,2,1);
       slice(X,Y,Z,sense,xslice,yslice,zslice);xlabel('Distance[m]');ylabel('Distance [m]');colorbar;title(colorbar,'Conductivity [mS/m]')
        title('Sensitivity plot');subplot(1,2,2);slice(X,Y,Z,p_trail,xslice,yslice,zslice);hold on;plot(doi,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10); 
     
        xlabel('Distance[m]');ylabel('Distance [m]');colorbar;title(colorbar,'Conductivity [mS/m]'); title('Inverted Model (S/m)');
        drawnow
       
   else
       slice(X,Y,Z,sense,xslice,yslice,zslice);xlabel('Distance[m]');ylabel('Distance [m]');colorbar;title(colorbar,'Conductivity [mS/m]')
       title('Inverted Model (S/m)');
       drawnow        
   end
%% Compute convergence criteria
    % reduce barrier function by order of magnitude/2 = 5
   barrier = barrier./5;%((rms_error(i)-rms_trial(j))/rms_error(i)) )
   % Evaluate error tolerance in rms errors
    aq1=sum(G);aq2=sum(G_try);aq3=find(aq1<-0.01);aq4=find(aq2<-0.01);  
   rms_error(isnan(rms_error))=100000000;
   rms_trial(isnan(rms_trial))=100000000;
   rms_error(aq3)=100000000;
   rms_trial(aq4)=100000000;
   
   if min(rms_trial) > err_tol
   %    index = find(min(rms_trial));
    %   best_rms =min(rms_trial);
     %  best_index = index;
     [best_rms, best_index] = min(rms_trial);
     else
       %rms_trial(isnan(rms_rial))=100000000;
       index = find(rms_trial <= err_tol);
        pnew = p_trial(:,index(1));
        if isempty(index)
            %rms_trial(isnan(rms_trial))=100000000000;
           [best_rms, best_index] = min(rms_trial);
           p_final = p_trial(:,best_index);
           rms_error(i+1)  = best_rms;
           G(:,i+1) = G_try(:,best_index);
           mu_final = mu(best_index);
           return;
       %else
       %best_index = index(1);
       %best_rms = rms_trial(best_index);
      %[best_rms, best_index] = min(rms_trial);

       end
       
   end
   
     pnew = p_trial(:,best_index);
     rms_vs_chg = norm((pnew' - po(:)')./po(:)')/sqrt(length(po(:)'));
 
  % p_final = pnew;
       %check for convergence
  % rms_vs_chg = norm((pnew' - po')./po')/sqrt(length(po'));
 if i==1
   if rms_vs_chg < tol_eca | (best_rms <= err_tol) | (i==max_iter)
       [best_rms, best_index] = min(rms_error);
       p_final= pin(:,:,best_index);
       %baa= shiftdim(pnew(find(wta)));
       %da = params;
       %da(find(wta)) = baa;%pnew;
       %pnew = da;
      % p_final = pnew;

       rms_error(i+1)  = best_rms;
       G(:,i+1) = G(:,best_index);
       %mu_final = muv;% mu(best_index);
      dlmwrite('inverted_models.dat',p_final);
      dlmwrite('predicted_data.dat',G);
      dlmwrite('rms_errors.dat',rms_error);
       fprintf('!!*******      Model Sequence Completed Successfully     **********!!\n')

      str=strcat('!!*******      Model Sequence Completed Successfully     **********!!');
      statusUpdate(str);
      progressbar(1);
       
     return;
   end
 else
      chg=(abs(rms_error(i-1)-rms_error(i)))/rms_error(i-1);
       if rms_vs_chg < tol_eca | (best_rms <= err_tol) | (i==max_iter) | chg<0.01 %| rms_error(i)>rms_error(i-1)
        %  rms_error(isnan(rms_error))=10000000000;
     [best_rms, best_index] = min(rms_error);
       pnew= (pin(:,:,best_index));
       baa= shiftdim(pnew(find(wta)));
       da = params;
       da(find(wta)) = baa;%pnew;
       pnew = da;
       p_final = pnew;
      rms_error(i+1)  = best_rms;
      G(:,i+1) = G(:,best_index);
     % mu_final = mu(best_index);
      mu_final = muv;% mu(best_index);
       
      dlmwrite('inverted_models.dat',p_final);
      dlmwrite('predicted_data.dat',G);
      dlmwrite('rms_errors.dat',rms_error);
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
            % k-corner is the index of the muh values (seminorm
            % equivalent=nsoln_norm)
            %new_mu = cornerval;
            new_mu = nsoln_norm(k_corner);
            if isempty(new_mu), new_mu=mu; end
            mu = [new_mu/20;new_mu/8;new_mu/4;4*new_mu;8*new_mu;20*new_mu];
            if length(muv)>length(muh), muv=mu; else muh=mu; end
        end
    end


%%
   po=reshape(pnew,M,P*PP);
end  % end loop (A1) over iterations 
  
