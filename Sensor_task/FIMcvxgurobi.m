%% cvx learning
f=randn(150000,1).^2;
binint
cvx_begin 
    variable mmu(10) binary
    maximize(customcvxfunc(mmu))
    subject to
        sum(mmu)==1
%         mmu>=0
%         mmu<=1
cvx_end

N=30;
opt=optimset('Algorithm','interior-point');
mmu=fmincon(@(x)-customcvxfunc(x),[1,zeros(1,N-1)],[],[],ones(1,N),1,zeros(1,N),ones(1,N),[],opt)

[m1,m2]=meshgrid(0:0.1:1);
p=zeros(size(m1));
for i=1:1:size(m1,1)
    for j=1:1:size(m1,2)
        p(i,j)=customcvxfunc([m1(i,j),m2(i,j)]);
    end
end

contour(m1,m2,p)
xlabel('m1')
ylabel('m2')

A=randn(5);
 P=magic(5)+1*(A*A').^2
 MU=[1,0,0;
     0,1,0;
     0,0,1];
 P=[eye(2),zeros(2,3);zeros(3,2),MU]*P*[eye(2),zeros(2,3);zeros(3,2),MU]'
 Px=P(1:2,1:2)
 Pz=P(3:5,3:5)
 Pxz=P(1:2,3:5)
 C1=Px-Pxz*inv(Pz)*Pxz'
 C2=Pz-Pxz'*inv(Px)*Pxz
 C2inv=inv(C2+0.001*eye(3))
 
 -trace(inv(Px)*Pxz*C2inv*Pxz')+trace(C2inv*Pz)
 
 -trace(inv(Px)*Pxz([1,2],[1,3])*inv(C2([1,3],[1,3]))*Pxz([1,2],[1,3])')+trace(inv(C2([1,3],[1,3]))*Pz([1,3],[1,3]))
 
 
 P=magic(5)+1*(A*A').^2
 Px=P(1:2,1:2)
 Pz=P(3:5,3:5)
 Pxz=P(1:2,3:5)
 C2=Pz-Pxz'*inv(Px)*Pxz
 C2invt=MU*inv(C2)*MU'
 Pxzt=Pxz*MU'
 Pzt=MU*Pz*MU'
 -trace(inv(Px)*Pxzt*C2invt*Pxzt')+trace(C2invt*Pzt)
 
 
 B=zeros(3);
 Px-Pxz*MU*inv(MU*Pz*MU+B)*MU*Pxz'
 det(Px-Pxz*MU*inv(MU*Pz*MU+B)*MU*Pxz')
 trace(Px-Pxz*MU*inv(MU*Pz*MU+B)*MU*Pxz')
 
 
 tf=48*60*60; % 48 hours 
dt=1*60; % 10 mins time step
t0=0;
Tvec=t0:dt:tf;
Nsat=size(Xsat0,1);
tic
 [r,v] = keplerUniversal(Xsat0(1,1:3)',Xsat0(1,4:6)',Tvec(end)-Tvec(1),398601.2);
toc
 opt = odeset('reltol',1e-12,'abstol',1e-12);
ytruth=cell(Nsat,1);
yplottruth=cell(Nsat,1);
for i=1:1
tic
    [tt,xx]=ode45(@twoBody,Tvec,Xsat0(i,:)',opt);
toc
    ytruth{i}=xx;
end
 

tic
[t,y]=twoBodyKeplerProp(Tvec,Xsat0(20,:));
toc

tic
    [tt,xx]=ode45(@twoBody,Tvec([1,150]),Xsat0(20,:)',opt);
toc