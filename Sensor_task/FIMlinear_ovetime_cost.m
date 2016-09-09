function g=FIMlinear_ovetime_cost(mus,FIM0,invPhis,H,invR,nt,nrads,nsats)
% FIM0 is a cell of the initial FIMs of all satellites FIM(sat)
% invPhis is the cell for STMs of all satellites for each time step Phis(sat,k)
%H is the cell of jacobian of measurement H(sat,rad,k)
%R is cell of covariance of radars R(rad)
f=0;

ss=1;
for i=1:1:nsats
    for k=1:1:nt
        P=zeros(6,6);
        for j=1:1:nrads
            P=P+mus(ss)*H{i,j,k}'*invR{j}*H{i,j,k};
            ss=ss+1;
        end
        FIM0{i}=invPhis{i,k}'*FIM0{i}*invPhis{i,k}+P;
    end
    f=f+FIM0{i};
end

g=trace(f);

end



