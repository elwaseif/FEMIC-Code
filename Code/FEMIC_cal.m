function [rawGem, calGem, filtGem]=FEMIC_cal(dcale,dobsy,f,r,q,vall)
%dcal = prdGEM; % data you are calibrating to (predicted response to the DC model)
%dobs = obsGEM; % data you want to calibrate (the observed GEM)
isOpen = matlabpool('size') > 0;%plotdoi=0;
fprintf('******  STARTING THE RAW DATA FILTERING AND CALIBRATION PROCESS  ******\n');
if isOpen==0
matlabpool('local',vall)
end
m=mean(dcale(:,4));
parfor ii=1:length(f)
	depth(ii)=503/sqrt(f(ii)*(m*1000.)); %estimate depth of penetration for frequency
	d2(:,ii)=abs(dcale(:,3)-depth(ii)); %find difference between depth of resistivity and depth of EM
	ff=min(d2(:,ii)) %which resisitivity depth is closest to the EM depth
	fdepth(:,ii)=find(d2(:,ii)==ff); %find all resistivity measurements that correspond to that depth
end
fdepth=fdepth';fdepth=fdepth(:);
dcale2(:,1)=dcale(fdepth,1);
dcale2(:,2)=dcale(fdepth,2);
dcale2(:,3)=dcale(fdepth,3);
dcale2(:,4)=dcale(fdepth,4);dcale=dcale2;
xx=(dcale(:,1));yy=(dcale(:,2));dd=dcale(:,3);
dcall=dcale(:,4);
dcali=reshape(dcall,length(f),length(dd)/length(f));
warning off
dobss=dobsy(:,8);
xx2=dobsy(:,1);
yy2=dobsy(:,2);
rawGem=dobss;
coor1=[xx,yy];
coor2=[xx2,yy2];
[m n]=size(coor1);
[ma na]=size(coor2);

%finding distance
for i=1:m
	parfor j=1:ma
		bo(j,i)=sum(abs(coor1(i,:)-coor2(j,:)));
	end
end

[mm nn]=size(bo);
for ii=1:nn
	ba=bo(:,ii);
	boo(ii)=min(find(ba==min(ba)));
end
dcali=dcali(:);
dcali2=reshape(dcali,length(f),length(dcali)/length(f));

[dcal]=FEMIC_forward2D(dcali2,depth,vall);
dcal=dcal';
dobsa=dobss(boo);
dobs=reshape(dobsa,length(f),length(dobsa)/length(f));
dobs=dobs';
ns = size(dobs,1); % soundings
nf = size(dobs,2); % frequencies
% initial values
G = ones(1,nf); % gain
B = zeros(1,nf); % bias
dprd = (dcal + repmat(B,ns,1))*diag(G);

wd = ones(1,nf); % enter 1/sigma values
wm = [1 1]; % G, B

% build block matrices - kind of klugey for now...
for i=1:nf   
    A(i).J = [(dcal(:,i) + B(i)) G(i)*ones(ns,1)];
    A(i).D = dobs(:,i) - dprd(:,i);
    A(i).Wd = wd(i)*eye(ns);
    A(i).Wm = diag(wm);
    A(i).m = [G(i);B(i)];
end
J = blkdiag(A.J);
D = blkdiag(A.D);
Wd = blkdiag(A.Wd);
Wm = blkdiag(A.Wm);
m = blkdiag(A.m);
m0 = blkdiag(A.m);


%% solve
for n = 1:10   
  
    H = J'*Wd'*Wd*J + Wm'*Wm;
    r = J'*Wd'*D + Wm'*Wm*(m-m0);
    
    % calculate update
    dm = H\r;
    m = m + dm;
    G = G + diag(dm(1:2:end,:)).';
    B = B + diag(dm(2:2:end,:)).';
    
    % update predicted & calibrated data
    dprd = (dcal + repmat(B,ns,1))*diag(G);
    dcor = dobs*diag(1./G) - repmat(B,ns,1);
	calGem=dcor;
    
    % update block matrices - kind of klugey for now...
    for i=1:nf
        A(i).J = [(dcal(:,i) + B(i)) G(i)*ones(ns,1)];
        A(i).D = dobs(:,i) - dprd(:,i);
    end
    J = blkdiag(A.J);
    D = blkdiag(A.D);
    
    % misfit
    err(n,:) = diag(D'*(Wd'*Wd)*D);
    
end
dcor=dcor';
dcor2=dcor(:);
dobss(boo)=dcor2;
%dlmwrite('dob.dat',dobss);


%% PCA filter
%dobss(ffo)=dcor2;
dcorr=reshape(dobss,length(f),length(dobss)/length(f));%length(find(xx2==xx2(1))),length(dobss)/length(find(xx2==xx2(1))));dcorr=dcorr';
m = mean(dcorr.',2);ns=length(f);
M = repmat(m,1,ns);
[u,s,v] = svd(dcorr.' - M);
lp = u(:,1)*u(:,1)'*(dcorr.'-M) + M; filtGem=lp;% 1st component only
Filtered_data=dobsy;Filtered_data(:,8)=lp(:);
%dlmwrite('Filtered_data.dat',Filtered_data,'precision','%12.5f');
dlmwrite('Filtered_data.dat',Filtered_data);
fprintf('******  THE RAW DATA HAS BEEN FILTERED AND CALIBRATED SUCCESSFULLY  ******\n');
