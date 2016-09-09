function [pnum,qden]=MI_num_dem_poly_det_zver(v,Px,Pxz,Pz,hn,nz,vindfun)
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

ntar=length(Px);

pnum=1;
qden=1;
for nt=1:1:ntar
PPz=mpol(Pz{nt});
PPxz=mpol(Pxz{nt});


    for i=1:1:nz
        PPxz(:,(i-1)*hn+1:i*hn)=v(vindfun(i,nt))*PPxz(:,(i-1)*hn+1:i*hn);
        for j=1:1:nz
            if j~=i
                PPz((i-1)*hn+1:i*hn,(j-1)*hn+1:j*hn)=v(vindfun(i,nt))*PPz((i-1)*hn+1:i*hn,(j-1)*hn+1:j*hn);
                PPz((j-1)*hn+1:j*hn,(i-1)*hn+1:i*hn)=v(vindfun(i,nt))*PPz((j-1)*hn+1:j*hn,(i-1)*hn+1:i*hn);
            end
        end
    end

pnum=pnum*det(PPz);
qden=qden*det(PPz-PPxz'*inv(Px{nt})*PPxz);

end


