function [p,q]=MI_num_dem_poly(v,Px,Pxz,Pz,hn,nz,pen)
%fn is the number state variables
%hn is the number meas variables in one set of measurements
% nz is the number of measurement sets.
% v=sym('v',[1,nz]);
% PPxz=sym('PPxz',size(Pxz));
% PPz=sym('PPz',size(Pz));
% mpol('x',1,nz)
%if pen ==0 then use the usual way where all the variables are distrinuted
%in Pz
% else use pen to weight the R matrix

PPz=mpol(Pz);
PPxz=mpol(Pxz);
PPx=mpol(Px);
Cpz=mpol(Pz);
if pen==0
    for i=1:1:nz
        PPxz(:,(i-1)*hn+1:i*hn)=v(i)*Pxz(:,(i-1)*hn+1:i*hn);
        for j=1:1:nz
            if j~=i
                PPz((i-1)*hn+1:i*hn,(j-1)*hn+1:j*hn)=v(i)*PPz((i-1)*hn+1:i*hn,(j-1)*hn+1:j*hn);
                PPz((j-1)*hn+1:j*hn,(i-1)*hn+1:i*hn)=v(i)*PPz((j-1)*hn+1:j*hn,(i-1)*hn+1:i*hn);
            end
        end
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

q=det(PPz);

for i=1:1:nz*hn
    for j=1:1:nz*hn
      Cpz(i,j)=Cofactor(PPz,size(Pz,1),i,j);         
    end
end

p=det(q*PPx-PPxz*Cpz*PPxz');
q=q^max(size(PPx));
