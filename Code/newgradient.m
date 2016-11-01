function [Gx, Gz]=newgradient(x,z)

%mesh.num_param - (x-1)*(z-1)
%tmp_x # of x parameters -1 (center values)
%tmp_z # of z parameters -1 (center values)


mesh.num_param=length(x)*length(z);
for i=1:length(z)
mesh.param_x(:,i)=x;
end
mesh.param_x=mesh.param_x(:);
for i=1:length(z)
  mesh.param_y(:,i)=z(i)*ones(length(x),1);
end
mesh.param_y=mesh.param_y(:);
c=zeros(mesh.num_param,mesh.num_param);
cx=zeros(mesh.num_param,mesh.num_param);
cy=zeros(mesh.num_param,mesh.num_param);

% tmp_x=union(mesh.tmp_param(:,1),mesh.tmp_param(:,1));
% tmp_y=union(mesh.tmp_param(:,2),mesh.tmp_param(:,2));


tmp_x=unique(mesh.param_x);
tmp_y=unique(mesh.param_y);



for i=1:mesh.num_param
    
    current_x=mesh.param_x(i);
    current_y=mesh.param_y(i);
    ind=find(tmp_x==current_x);
    % search all other parameters that have the same y and the x=ind+1
    for j=1:mesh.num_param
        if ind~=length(tmp_x)
            if mesh.param_y(j)==current_y && mesh.param_x(j)==tmp_x(ind+1) 
               cx(i,j)=1;  
%                cx(i,j)=sqrt(      (mesh.tmp_param(j,6)-mesh.tmp_param(j,5))/ ( mesh.tmp_param(j,1)-mesh.tmp_param(i,1)));
            end
        end
    end
end

for i=1:mesh.num_param
   cx(i,i)=-sum(cx(i,:));    
end


ctc1=cx'*cx;
        

for i=1:mesh.num_param
    
    current_x=mesh.param_x(i);
    current_y=mesh.param_y(i);
    ind=find(tmp_y==current_y);
    % search all other parameters that have the same y and the x=ind+1
    for j=1:mesh.num_param
        if ind~=length(tmp_y)
            if mesh.param_y(j)==tmp_y(ind+1) && mesh.param_x(j)==current_x 
               cy(i,j)=1;  
%                cy(i,j)=sqrt(      (mesh.tmp_param(j,4)-mesh.tmp_param(j,3))/ ( mesh.param_y(j)-mesh.param_y(i)) );
            end
        end
    end
end

for i=1:mesh.num_param
   cy(i,i)=-sum(cy(i,:));    
end
Gx=cx;
Gz=cy;