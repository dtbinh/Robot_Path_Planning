function [hx,GG]=comb_meas_props(x,h,g,R,SS)
% combine all the measurement properties into one big ass measurement
% function with big ass R and g for Kalman filtering
%SS has to be defined in a specif manner

nz=size(SS,1);
hx=[];
GG=[];
% keyboard
for i=1:1:nz
    hx=vertcat(hx,h{SS(i,1)}(x,SS(i,2:3)));
    GG=blkdiag(GG,g{SS(i,1)}(x,SS(i,2:3))*R{SS(i,1)}*g{SS(i,1)}(x,SS(i,2:3))');
end


end










