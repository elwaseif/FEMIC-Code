% function [out]=FEMIC_forward2D(params,d,f,r,q,mdltype)
%
% This MATLAB code calls the FEMIC forward modeling routines to compute
% observed conuctivity values for FDEM instruments.
%
% G. Schultz 2008
%%
% INPUTS:
%   params = model parameters to be optimizaed
%   d = depths (initial) (1 x mlayers)
%   pobs = measurements (expected to be VDM cat HDM both (1 x 2*length(f))
%   sigma = standard deviations on measurements (1 x 2*length(f))
%   f = frequencies
%   r = separation distances between Rx and Tx for bistatic case
%   q = [=1,2,3] to designate the data types represented in the pobs input
%       array: 1=Vertical Magnetic Dipole only
%              2=Horizontal Magnetic Dipole only
%              3=Both VMD and HMD data
%%
% OUTPUTS:
%   p_final = the final model array in [P stations x N frequencies (2xN for
%
% EXAMPLE USAGE:
% [g]=FEMIC_forward2D(params,d,f,r,q,mdltype);

function [out]=FEMIC_forward2D(params,d,vall)

%% Initialize array size and constant parameters
szp=size(params);       % number of stations
P=szp(2);               
%M=length(d);     
M = szp(1);             % number of model parameters
Md = length(d);         % number of layer thicknesses
Ms = M-Md;              % should be the number of conducitivities
%N=szp(1);
load SS
N=length(S.freq);            % number of frequencies 
%if q==3,
 %   NN=2*N;             % twice the frequenices if both VDM and HDM used
%else                    % number of data (freq x 2)
 %   NN=N;               % total number of data points (same as N if only one orientation of data is used)
%end
%NP=NN*P;                % size of the data set
%MP=M*P;                 % size of the model output set
%dlmwrite('q.txt',q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Initialize the regularization matrices
po=params(1:M,:);
% Create a NPxNP diagonal weighting matrix whose elements are the inverse of the
% estimated (or measured) stanard deviations of the measurements
dp=1;
tic;
isOpen = matlabpool('size') > 0;%plotdoi=0;
if isOpen==0
matlabpool('local',vall)
end
%dlmwrite('d.dat',d);
M.k = length(d); % number of layers
M.thk = d'; % layer thicknesses (m), enter zero for half-space layer
%M.con = [.01;.001;.1]; % layer conductivities (S/m)
M.chie = zeros(length(S.freq),1); % layer dielectric susceptibilities (SI) eps = eps_0 * (1 + chie)
M.chim = zeros(length(S.freq),1); % layer magnetic susceptibilities (SI)   mu = mu_0 * (1 + chim)
el = -1;
%% Outer loop over maximum number of iterations (breaks if alternate
%% convergence criteria are met
warning off
   for ii=1:P,          % loop over all measurement positions (B1)
      % fprintf(' Calibration and filtering the raw data are in progress.. %3.0f\n')
       pin=[po(:,ii)]; %da=d(:,ii); % model input parameters for this position (ii)
       M.con=(1./pin);
       [Gout] = fdem1dfwd(S,M,el,0);  
       %[GG(:,ii),JJ]=jacnmodelMdl_mat(f,r,pin,d);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Compute the model kernel Gout, and Jacobian matrix J
 %      [Gout] =  FEMIC_JacobianMdl(f,r,pin,d,N,M,Md,2);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %[t,x,y] = sim('jacobianMdl3');
       mdl(:,ii) = Gout;       % place current Gout vector in the global GG array
   end  % end (B1) loop over measurement positions
   
   elap_time=toc/60;
   mdlIdx = mdl(:,1)>0;
   j=0;
   for i=1:length(mdlIdx),
       if mdlIdx(i)~=0,
           j=j+1;
           out(j,:)=mdl(i,:);
       end
   end
  % fo=length(out);
   %if q==2
    %   out(length(f)+1:fo)=0;
   %end