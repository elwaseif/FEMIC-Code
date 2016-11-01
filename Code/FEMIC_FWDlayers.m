% FUNCTION FEMIC_Jacobian.m
%
% This function computes the forward layered-Earth conductivity model for a
% given set of bistatic frequency-domain EMI data.  The forward model is
% coded to follow the formulation for the classic layered-Earth recursion
% given in texts and papers by Keller and Frischnecht, Koefed, Wait, and
% others.  The recursion relationship uses the trigonemetric form and
% relies on computation of the Hankel transform.  This bessel-kernel
% integral is solved by an effecient gauss-quadrature numerical integration
% function over a variable length of integration points - the length is set
% by a functional formula that is mostly dependent on the induction number
% and number of layers to be solved for.  
% 
% Originated by: Greg Schultz
% Modified from original codes produced in 2004/2005
% Significantly modifield in 2008
% Code delivered to the USGS under PO XXXXXXXXX in June 2008
%
%%
% INPUTS:
%   freq = frequencies
%   rspa = separation distances between Rx and Tx for bistatic case
%   params = model parameters (layer conductivities and layer thicknesses
%   d = depths (initial) (1 x mlayers)
%   N = number of data per station
%   M = number of model parameters
%   nlayers = number of layers to solve for
%   q = [=1,2,3] to designate the data types represented in the pobs input
%       array: 1=Vertical Magnetic Dipole only
%              2=Horizontal Magnetic Dipole only
%              3=Both VMD and HMD data
%%
% OUTPUTS:
%   G = solution to the model system [P stations x N frequencies (2xN for
%       both VMD and HMD data)] form
%
% EXAMPLE USAGE:
% [G] = FEMIC_FWDlayers(freq,rspa,params,d,N,M,nlayers,q)

function [G] = FEMIC_FWDlayers(freq,rspa,params,d,N,M,nlayers,q)
%tic;
%% Define constants
aint1 = 0.0;
aint2 = 1.3799;
gnpts = 100;
muo = 1.25663706143591e-6;
pi = 3.14159265358979;
%% Pre-initialize Parameters
% p-locations,
% m-layers
% n-frequencies or spacings
%params = zeros(M,1);% conductivities
% d = zeros(M,1);     % thicknesses
% G = zeros(2*N,1);   % model matrix
% J=zeros(2*N,M);     % Jacobian
% B= zeros(N,1);      % induction no.
% psigma = zeros(M,1);% modeled conductivities
% sigma = zeros(M,1); % modeled conductivities initial

G=zeros(30,1);
B=zeros(30,1);
psigma = zeros(30,1);
sigma = zeros(30,1);
% use S/m instead of input mS/m (Assumes input in mS)
for ii=1:nlayers,
    sigma(ii) = params(ii)./1000.0;
    lyr(ii) = d(ii);
end
% Intialize layer thicknesses
if M>nlayers,
    for ii=1:M-nlayers
        lyr(ii) = params(nlayers+ii);
    end
end
spac = 1.0;         % gain factor (never used - always 1 - debug mode only)


%% Loop over all frequencies to compute the rows of G and J
for k=1:N,
    % radial frequency
    w = 2*pi*freq(k);
    r2 = rspa(k)*spac;     % spacing
    skin=sqrt(2.0/(max(sigma)*muo*w));% skin depth
    % Auto-selection of Hankel transform integration limit
	B(k)=r2/skin;
    gnpts = 100+20.*ceil(((nlayers.^1.5).*sqrt(B(k).^3)));
    if gnpts>2500; gnpts = 2500; end
    % Compute the T-values 
    if q~=2,    % if NOT HMD only data - thus, if VMD or both VMD and HMD
        [tval0] = gaussinterp(r2,sigma,0,w,lyr,N,nlayers,aint1,aint2,gnpts);
        % complex magnetic field ratios for horizontal coplanar and coaxial
        % coil configurations (VMD) ^
        s1 = 1-tval0(1)*(r2^3.0);
        s2 = -tval0(2)*(r2^3.0);
    end
    
%time1 = toc; fprintf('Gauss 1 time = %f\n',time1);tic;
    if q~=1,
        [tval2] = gaussinterp(r2,sigma,2,w,lyr,N,nlayers,aint1,aint2,gnpts);
        % complex magnetic field ratios for vertical coplanar and coaxial
        % coil configurations (HMD) >>
        s3 = 1.0-tval2(1)*(r2^2.0);
        s4 = -tval2(2)*(r2^2.0);
    end
%time1 = toc; fprintf('Gauss 2 time = %f\n',time1);tic;
    % 2*i*B^2 --> complex
    % twoiB2=2.0*1i*B(k)^2;
    %% compute the G system matrix
    if q~=2,
        G(k) = s2*(4000.0/(muo*w*(r2^2.0)));
    end
    if q~=1,
        G(k+N) = s4*(4000.0/(muo*w*(r2^2.0)));
    end
    if q==2,
        G(k) = s4*(4000.0/(muo*w*(r2^2.0)));
    end

end   % end loop over frequencies
%return
%time1 = toc; fprintf('Total time = %f\n',time1);
%% cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function [val] = gaussinterp(r2,sigma,ntype,w,d,N,M,aint1,aint2,gnpts)
% function to integrate the hankel transform via gauss-quadrature

% intialize variables
xabsc=zeros(gnpts,1);
gwts = zeros(gnpts,1);
ffr = zeros(gnpts,1);
ffi = zeros(gnpts,1);
val = zeros(2,1);

% set local variables
epsd = 3.0e-14;
x1 = 1.0;
x2  = -1.0;
mm = (gnpts+1)/2;
xm = 0.5*(x2+x1);
xl = 0.5*(x2-x1);
pVec = zeros(1002,1);
%tic; 
% Loop over all integration points
for ii=1:mm,
    z=cos(pi*(ii-0.250)/(gnpts+0.50));
    z1 = 0;
%     p1=1.0;
%     p2=0.0;
%     for jj=1:gnpts,
%         p3 = p2;
%         p2=p1;
%         p1 = ((2.0*jj-1.00)*z*p2-(jj-1.0)*p3)/jj;
%     end
%     pp = gnpts*(z*p1-p2)/(z*z-1);
%     z1=z;
%     z=z1-p1/pp;

%     pVec(1) = 0;
%     pVec(2) = 1;
%     for jj = 1:gnpts
%         pVec(jj+2) = ((2.0*jj-1.00)*z*pVec(jj+1)-(jj-1.0)*pVec(jj))/jj;
%     end
%     pp = gnpts*(z*pVec(gnpts+2)-pVec(gnpts+1))/(z*z-1);
%     z1=z;
%     z=z1-pVec(gnpts+2)/pp;
    
%     while (abs(z-z1)>epsd),
%         p1=1.0;
%         p2=0.0;
%         for jj=1:gnpts,
%             p3 = p2;
%             p2=p1;
%             p1 = ((2.0*jj-1.00)*z*p2-(jj-1.0)*p3)/jj;
%         end
%         pp = gnpts*(z*p1-p2)/(z*z-1);
%         z1=z;
%         z=z1-p1/pp;
%     end
    pVec(1) = 0;
    pVec(2) = 1;
    while (abs(z-z1)>epsd),
        for jj = 1:gnpts
             pVec(jj+2) = ((2.0*jj-1.00)*z*pVec(jj+1)-(jj-1.0)*pVec(jj))/jj;
        end
        pp = gnpts*(z*pVec(gnpts+2)-pVec(gnpts+1))/(z*z-1);
        z1=z;
        z=z1-pVec(gnpts+2)/pp;
    end
    xabsc(ii)=-(xm+x1*z);
    xabsc(gnpts+1-ii)=xm+x1*z;
    gwts(ii) =2.0*x1/((1.0-z*z)*pp*pp);
    gwts(gnpts+1-ii) = 2.0*x1/((1.0-z*z)*pp*pp);
end
%time1 = toc; fprintf('---time = %f\n',time1);tic; 
for ii=1:gnpts,
    x=(aint1+aint2)/2+xabsc(ii)*(aint2-aint1)/2.0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% EValuate the recursion relationship
    valu = evlfun(x,r2,sigma,ntype,w,d,N,M);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ffr(ii) = real(valu)*(aint2-aint1)/2.0;
    ffi(ii) = imag(valu)*(aint2-aint1)/2.0;
end

%time1 = toc; fprintf('---time = %f\n',time1);
%[size(gwts) size(ffr)]
val(1) = dot(gwts,ffr);
val(2) = dot(gwts,ffi);
%return

%% cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function [valu] = evlfun(x,r2,sigma,ntype,w,d,N,M)
    
% this function evaluates the recursion relationship to generate the valu
% of Ro given a set of layer conducitivities and thicknesses

% MAgnetic permeability of free space
muo = 1.25663706143591e-6;
% intialize the squared integration parameter
l2 = 0.0;
% nu2 = complex(0.*[1:20],1.*[1:20]);
% v = complex(0.*[1:20],1.*[1:20]);
% rtop = complex(0.*[1:20],1.*[1:20]);
% rbot = complex(0.*[1:20],1.*[1:20]);
% ro = complex(0.*[1:20],1.*[1:20]);
% vik = complex(0.*[1:20],1.*[1:20]);
% initial the complex intermediate recursion variables
nu2 = 0.*[1:20]'+1.*[1:20]'*1i;
v = 0.*[1:20]'+1.*[1:20]'*1i;
rtop = 0.*[1:20]'+1.*[1:20]'*1i;
rbot = 0.*[1:20]'+1.*[1:20]'*1i;
ro = 0.*[1:20]'+1.*[1:20]'*1i;
vik = 0.*[1:20]'+1.*[1:20]'*1i;
M=M+1;
% Add a an air layer of 0 conducitivity
sigma = [0; sigma];
% Add an air layer (boundary) with zero thickness to ensure convergence of
% the solution - IMPORTANT
d=[0 d];
% initialize the wavenumber, nu
for ij=1:M,
    nu2(ij)=(0+1i)*w*muo*sigma(ij);
end

% Compute the real and imaginary propagation terms (exponents)
l2=x^2;
for ij=1:M,
    vi=imag(sqrt(l2+nu2(ij)));
    vr=real(sqrt(l2+nu2(ij)));
    v(ij)=vr+1i*vi;
end

% there was a loop of initializiation of the complex variables here in the
% old FORTRAN code - not needed in matlab b/c auto-initialization
ro(M) = 0;

% Recurison relationship - starts at depth and peels back to surface
% for ij=M:-1:2,
%     vik(ij) = (v(ij-1)-v(ij))/(v(ij-1)+v(ij));
%     rtop(ij-1)=vik(ij)+ro(ij)*exp(-2.0*d(ij)*v(ij));
%     rbot(ij-1)=1.0+vik(ij)*ro(ij)*exp(-2.0*d(ij)*v(ij));
%     ro(ij-1)=rtop(ij-1)/rbot(ij-1);
% end

% Recurison relationship - starts at depth and peels back to surface
for ij=M:-1:2,
    va(ij) = (v(ij-1)-v(ij));
    vb(ij) = (v(ij-1)+v(ij));
    vc(ij) = ro(ij)*exp(-2.0*d(ij)*v(ij));
    ro(ij-1) = (va(ij)+vb(ij)*vc(ij))/(vb(ij)+va(ij)*vc(ij));
end

% this defines the sommerfeld integral order
% ntype = 0: need for horizonatal coplanar loops (vert. dipoles) and
%            vertical coaxial loops (with #2)
% ntype = 1: need for perpendicular loops
% ntype = 2: need for vertical coplanar loops (horiz. dipoles)
argu=r2*x;
[bj0, bj1] = jbessel(argu);
if ntype == 0
    valur=12*real(ro(1))*bj0;
    valui=l2*imag(ro(1))*bj0;
elseif ntype == 1
    valur=12*real(ro(1))*bj1;
    valui=l2*imag(ro(1))*bj1; 
elseif ntype==2,
	valur=x*real(ro(1))*bj1;
    valui=x*imag(ro(1))*bj1;
end
valu=valur+1i*valui;
%return


%% cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function [bj0,bj1]=jbessel(x)
        %       =======================================================
        %       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
        %                Y1(x), and their derivatives
        %       Input :  x   --- Argument of Jn(x)& Yn(x,x ò 0)
        %       Output:  BJ0 --- J0(x)
        %                DJ0 --- J0'(x)
        %                BJ1 --- J1(x)
        %                DJ1 --- J1'(x)
        %                BY0 --- Y0(x)
        %                DY0 --- Y0'(x)
        %                BY1 --- Y1(x)
        %                DY1 --- Y1'(x)
        %       =======================================================
        a=zeros(1,12);
        b=zeros(1,12);
        a1=zeros(1,12);
        b1=zeros(1,12);
        pi=3.141592653589793d0;
        rp2=0.63661977236758d0;
        x2=x.*x;
        % if(x == 0.0d0);
        %     bj0=1.0d0;
        %     bj1=0.0d0;
        %     dj0=0.0d0;
        %     dj1=0.5d0;
        %     by0=-1.0d+300;
        %     by1=-1.0d+300;
        %     dy0=1.0d+300;
        %     dy1=1.0d+300;
        %     return;
        % end;
        if(x <= 12.0d0);
            bj0=1.0d0;
            r=1.0d0;
            for  k=1:30;
                r=-0.25d0.*r.*x2./(k.*k);
                bj0=bj0+r;
                if(abs(r)< abs(bj0).*1.0d-15)break; end;
            end;
            bj1=1.0d0;
            r=1.0d0;
            for  k=1:30;
                r=-0.25d0.*r.*x2./(k.*(k+1.0d0));
                bj1=bj1+r;
                if(abs(r)< abs(bj1).*1.0d-15)break; end;
            end;
            bj1=0.5d0.*x.*bj1;
            %ec=log(x./2.0d0)+0.5772156649015329d0;
            cs0=0.0d0;
            w0=0.0d0;
            r0=1.0d0;
            for  k=1:30;
                w0=w0+1.0d0./k;
                r0=-0.25d0.*r0./(k.*k).*x2;
                r=r0.*w0;
                cs0=cs0+r;
                if(abs(r)< abs(cs0).*1.0d-15)
                    break;
                end;
            end;
            %by0=rp2.*(ec.*bj0-cs0);
            cs1=1.0d0;
            w1=0.0d0;
            r1=1.0d0;
            for  k=1:30;
                w1=w1+1.0d0./k;
                r1=-0.25d0.*r1./(k.*(k+1)).*x2;
                r=r1.*(2.0d0.*w1+1.0d0./(k+1.0d0));
                cs1=cs1+r;
                if(abs(r)< abs(cs1).*1.0d-15)
                    break;
                end;
            end;
            %by1=rp2.*(ec.*bj1-1.0d0./x-0.25d0.*x.*cs1);
        else
            a(:)=[-.7031250000000000d-01,.1121520996093750d+00,-.5725014209747314d+00,.6074042001273483d+01,-.1100171402692467d+03,.3038090510922384d+04,-.1188384262567832d+06,.6252951493434797d+07,-.4259392165047669d+09,.3646840080706556d+11,-.3833534661393944d+13,.4854014686852901d+15];
            b(:)=[.7324218750000000d-01,-.2271080017089844d+00,.1727727502584457d+01,-.2438052969955606d+02,.5513358961220206d+03,-.1825775547429318d+05,.8328593040162893d+06,-.5006958953198893d+08,.3836255180230433d+10,-.3649010818849833d+12,.4218971570284096d+14,-.5827244631566907d+16];
            a1(:)=[.1171875000000000d+00,-.1441955566406250d+00,.6765925884246826d+00,-.6883914268109947d+01,.1215978918765359d+03,-.3302272294480852d+04,.1276412726461746d+06,-.6656367718817688d+07,.4502786003050393d+09,-.3833857520742790d+11,.4011838599133198d+13,-.5060568503314727d+15];
            b1(:)=[-.1025390625000000d+00,.2775764465332031d+00,-.1993531733751297d+01,.2724882731126854d+02,-.6038440767050702d+03,.1971837591223663d+05,-.8902978767070678d+06,.5310411010968522d+08,-.4043620325107754d+10,.3827011346598605d+12,-.4406481417852278d+14,.6065091351222699d+16];

            t1=x-0.25d0.*pi;
            p0=1.0d0;
            q0=-0.125d0./x;
            k0=12;
                if(x >= 35.0)k0=10; end;
                if(x >= 50.0)k0=8; end;            
                for  k=1:k0;
                    p0=p0+a(k).*x.^(-2.*k);
                    q0=q0+b(k).*x.^(-2.*k-1);
                end;  %k=k0+1;
            cu=sqrt(rp2./x);
            bj0=cu.*(p0.*cos(t1)-q0.*sin(t1));
            %by0=cu.*(p0.*sin(t1)+q0.*cos(t1));
            t2=x-0.75d0.*pi;
            p1=1.0d0;
            q1=0.375d0./x;
            for  k=1:k0;
                p1=p1+a1(k).*x.^(-2.*k);
                q1=q1+b1(k).*x.^(-2.*k-1);
            end;  %k=k0+1;
            cu=sqrt(rp2./x);
            bj1=cu.*(p1.*cos(t2)-q1.*sin(t2));
            %by1=cu.*(p1.*sin(t2)+q1.*cos(t2));
        end
        % dj0=-bj1;
        % dj1=bj0-bj1./x;
        % dy0=-by1;
        % dy1=by0-by1./x;
        %return;