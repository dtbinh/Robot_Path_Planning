%% script to optimize the grid length
Vmin=10;
Vmax=15;
umin=-1;
umax=1;


for p1=-1:2:1
    for p2=-1:2:1
        L=dubminpaths(r0,rf,X0,Xf,p1,p2,dir)
    end
end