%input should be # cal points x # frequencies
load GEMCalWorkspace

x=calX; % x coordinate for co-located GEM & DC
dcal = prdGEM; % data you are calibrating to (predicted response to the DC model)
dobs = obsGEM; % data you want to calibrate (the observed GEM)

ns = size(dobs,1); % soundings
nf = size(dobs,2); % frequencies

% initial values
G = ones(1,nf); % gain
B = zeros(1,nf); % bias
dprd = (dcal + repmat(B,ns,1))*diag(G);

% weights  - leave at ones for now... it works, but should be fixed
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

%% PCA filter
m = mean(dcor.',2);
M = repmat(m,1,ns);
[u,s,v] = svd(dcor.' - M);
lp = u(:,1)*u(:,1)'*(dcor.'-M) + M; % 1st component only

%% Display
for i = 1:nf
    figure
    subplot(2,3,[1 2]);hold on;
    ylabel('ppm');title(['I - ' num2str(f(i))]);
    plot(x,real(dobs(:,i)),'b.-')
    %plot(flipud(gemXfilt),real(filtGEMobs(:,i)),'g.-')
    plot(x,real(dcal(:,i)),'r.-')
    plot(x,real(dcor(:,i)),'m.-')
    plot(x,real(lp(i,:)),'c.-')
    legend({'observed GEM','prd GEM from DC','calibrated GEM','calibrated & filtered GEM'})
    xlim([0 max(x)])
    subplot(2,3,[4 5]);hold on;
    xlabel('distance (m)');ylabel('ppm');title(['Q - ' num2str(f(i))]);
    plot(x,imag(dobs(:,i)),'b.-')
    plot(x,imag(dcal(:,i)),'r.-')
    plot(x,imag(dcor(:,i)),'m.-')
    plot(x,imag(lp(i,:)),'c.-')
    xlim([0 max(x)])
    subplot(2,3,3);hold on
    plot(angle(dobs(:,i))*180/pi,abs(dobs(:,i)),'b.')
   % plot(angle(filtGEMobs(:,i))*180/pi,abs(filtGEMobs(:,i)),'g.')
    plot(angle(dcal(:,i))*180/pi,abs(dcal(:,i)),'r.')
    plot(angle(dcor(:,i))*180/pi,abs(dcor(:,i)),'m.')
    plot(angle(lp(i,:))*180/pi,abs(lp(i,:)),'c.')
    xlabel('phase');ylabel('amplitude')
    title(['G = ' num2str(abs(G(i))) ', \phi = ' num2str(angle(G(i))*180/pi) ', B = ' num2str(B(i))]);
    subplot(2,3,6);hold on
    plot(dobs(:,i),'b.')
    %plot(filtGEMobs(:,i),'g.')
    plot(dcal(:,i),'r.')
    plot(dcor(:,i),'m.')
    plot(lp(i,:),'c.')
    xlabel('in phase');ylabel('quadrature')
    %title(['G_p_c_a = ' num2str(abs(G(i))) ', \phi_p_c_a = ' num2str(angle(G(i))*180/pi) ', B_p_c_a = ' num2str(B(i))])
end
    

