function [P,Q]=MI_num_dem_poly_zver_multi_x(v,Px,Pxz,Pz,hn,nz,nx,pen)
%fn is the number state variables
%hn is the number meas variables in one set of measurements
% nz is the number of measurement sets.
% v=sym('v',[nx,nz]); where nx is the number of sensors
% Pz is cell for all sensors that ca
% Px is cell for each object x
% Similarly Pxz is a cell for each object
% mpol('x',1,nz)
%if pen ==0 then use the usual way where all the variables are distrinuted
%in Pz
% else use pen to weight the R matrix

mpol('q',nx,1);
mpol('p',nx,1);

if pen==0
   
    for nxi=1:1:nx
   PPz=mpol(Pz{nxi});
   PPxz=mpol(Pxz{nxi});

    for i=1:1:nz
        PPxz(:,(i-1)*hn+1:i*hn)=v(nxi,i)*PPxz(:,(i-1)*hn+1:i*hn);
        for j=1:1:nz
%             if j~=i  % no j~=i gives the correct trace to the TRUE system
                PPz((i-1)*hn+1:i*hn,(j-1)*hn+1:j*hn)=v(nxi,i)*PPz((i-1)*hn+1:i*hn,(j-1)*hn+1:j*hn);
                PPz((j-1)*hn+1:j*hn,(i-1)*hn+1:i*hn)=v(nxi,i)*PPz((j-1)*hn+1:j*hn,(i-1)*hn+1:i*hn);
%             end
        end
    end
    q(nxi,1)=trace(PPz-PPxz'*inv(Px{nxi})*PPxz);
    p(nxi,1)=trace(PPz);
    end
else
    
    R=[];
    for i=1:1:nz
        PPxz(:,(i-1)*hn+1:i*hn)=v(i)*Pxz(:,(i-1)*hn+1:i*hn);
        V=v(i)*eye(hn);
           R=blkdiag(R,V); 
    end
    PPz=Pz+pen*(eye(nz*hn)-R);
end

% p=trace(PPz);
% q=trace(PPz-PPxz'*inv(Px)*PPxz);

% for i=1:1:nz*hn
%     for j=1:1:nz*hn
%       Cpz(i,j)=Cofactor(PPz,size(Pz,1),i,j);         
%     end
% end

% p=trace(q*PPx-PPxz*Cpz*Pxz');

% mpol Q;
% mpol P;
% mpol SS;
% P=1;
% Q=0;
% for i=1:1:nx
%     P=P*p(i,1);
%     SS=1;
%     for j=1:1:nx
%         if j~=i
%         SS=SS*p(j,1);
%         end
%     end
%     Q=Q+q(i,1)*SS;
% end


mpol Q;
mpol P;
P=1;
Q=0;
for i=1:1:nx
    Q=Q+(p(i,1)-q(i,1))^1;
end
