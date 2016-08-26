% grid movement of UAV way points.
function ni=gridmove(Nn,pi,u)

% control u=1-right, 2-up, 3-left, 4-down
if u==1
    ni=pi+1;
elseif u==3
    ni=pi-1;
elseif u==2
    ni=pi+Nn;
elseif u==4
    ni=pi-Nn;
end

