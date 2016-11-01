function [doi]=FEMIC_DOI(d,pobs,sigma,f,r,muv,err_tol,max_iter,q,vall,perc)
fprintf('******  STARTING MAXIMUM DOI ESTIMATION  ******\n');
[m n]=size(pobs);sigma=sigma(:,1);
parfor i=1:n-1
po(i)=abs(sum(pobs(:,i))-sum(pobs(:,i+1)))\min(min([sum(pobs(:,i)) sum(pobs(:,i+1))]));
end
po(po==inf)=0;po=100.*po;
ee=find(po>perc);dff=n-length(ee);
if dff>1 && min(ee)~=1 && max(ee)~=n
    ff=[1 min(ee)-1 ee max(ee)+1 n];%ff(length(ff)+1:length(ee)+dff)=[ee max(ee)+1 n];
else
    ff=[ee];
end
ff=unique(ff);
if isempty(ff) || length(ff)<=3
    redpobs=pobs(:,[1 floor(n/3) floor(n/2) n]);
    ff =[1  floor(n/3) floor(n/2) n];
else
    redpobs=pobs(:,ff);
end
params=40*ones(length(f),1);%round((min(min(redpobs)) + max(max(redpobs)))/2)*ones(length(f),1);
   [sa ta]=size(redpobs);
   %model1=zeros(sa,ta);
  parfor u=1:ta
  [model1, muh_final, rms_error, G]=FEMIC_tost22(params,d,redpobs(:,u),sigma,f,r,muv,...
    err_tol,max_iter,q);
  model11(:,u)=model1;
  end
    params2=60*ones(length(f),1);%abs(params(1)-round(max(max(redpobs))/2))*ones(length(f),1); %model2=zeros(sa,ta);
  parfor u=1:ta
    [model2, muh_final, rms_error, G]=FEMIC_tost22(params2,d,redpobs(:,u),sigma,f,r,muv,...
    err_tol,max_iter,q);
model22(:,u)=model2;
  end
R=abs(model11-model22)./min(min(abs(params-params2)));
[u y]=size(R);dd=cumsum(d);
for i=1:y
    
    bl=find(R(1:u,i)>0.2);
    to(i)=dd(bl(1));
end
ap=1:n;%rq=length(ff);ff(rq+1)=n;to(rq+1)=to(1);
doi=interp1(ff,to,ap);
fprintf('******  END OF DOI ESTIMATION  ******\n');
fprintf('******  STARTING THE INVERSION  ******\n');