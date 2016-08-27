function xk1=KIRB_UM_eg_dyn_disc_manv_3(xk,T,XY,dx)
% global T;
% T=para(1);
% T=1;
% keyboard
% weight the direction by the measurmeents

%% nudge velocity of this sigma point

v=xk(3:4);
v=v(:);
vn=v/norm(v);


V=get_side_ind2(xk,XY,75);


Vorig=V;
ind=find(V(:,end)>=0);
V=V(ind,:);

try
    mmu=V(:,3)-Vorig(1,3);
catch
    keyboard
end
mmu=sum(mmu)/length(mmu);


alpha1 =sign(mmu);
V(:,3)=V(:,3)/sum(V(:,3));
thmu=sum([V(:,3).*V(:,1),V(:,3).*V(:,2)],1);
thmu=thmu/norm(thmu);
alpha=0.5;
if alpha1==0 || sum(isnan(thmu))>0
    Newdir=vn;
else
    Newdir=((1-alpha)*vn+alpha* thmu');
end
Newdir=Newdir/norm(Newdir);

xk(3:4)=norm(v)*Newdir;



%%
xk1=[1,0,T,0;
    0,1,0,T;
    0,0,1,0;
    0,0,0,1]*xk;
end
