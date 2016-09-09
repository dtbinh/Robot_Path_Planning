 clc
clear
 close all
% P=magic(5)*magic(5)'/100;
Px=[2,-1;-1,4.5];
mx=[3,5];
% Z1=@(x)[norm([x(1)-2,x(2)-2]);atan2(x(2)-2,x(1)-2)];
% Z2=@(x)[norm([x(1)+2,x(2)+2]);atan2(x(2)+2,x(1)+2)];
% Z3=@(x)[norm([x(1)-10,x(2)]);atan2(x(2),x(1)-10)];

Z1=@(x)[x(1);x(2)^2];
Z2=@(x)[x(1)+x(2);x(1)-x(2)];
Z3=@(x)[x(1)^3;3*x(1)+x(2)];

hjack=@(x,a,b)[(x(1)-a)/norm([x(1)-a,x(2)-b]),(x(2)-b)/norm([x(1)-a,x(2)-b]);(b-x(2))/norm([x(1)-a,x(2)-b])^2,(-a+x(1))/norm([x(1)-a,x(2)-b])^2];
% R1=diag([1^2,(1*pi/180)^2]);
% R2=diag([5^2,(0.1*pi/180)^2]);
% R3=diag([0.5^2,(10*pi/180)^2]);

R3=diag([0.1^2,(0.1)^2]);
R2=diag([200^2,(200)^2]);
R1=diag([3^2,(3)^2]);

hn=2;
nz=3;
H=[hjack(mx,2,2);hjack(mx,-2,-2);hjack(mx,10,0)];
H1=hjack(mx,2,2);
H2=hjack(mx,-2,-2);
H3=hjack(mx,10,0);



% [xut,wut]=UT_sigmapoints(mx',Px,2);
% z1ut=zeros(length(wut),2);
% for i=1:1:length(wut)
%     z1ut(i,:)=Z1(xut(i,:));
% end
% [x,w]=UT_sigmapoints(mx',Px,2);
[x,w]=GH_points(mx',Px,15);
z1=zeros(length(w),2);
z2=zeros(length(w),2);
z3=zeros(length(w),2);
for i=1:1:length(w)
    z1(i,:)=Z1(x(i,:));
    z2(i,:)=Z2(x(i,:));
    z3(i,:)=Z3(x(i,:));
end

Z=[z1,z2,z3];
Hlin=StatLinearize_sigmapts(Z,x,w);
H1lin=Hlin(1:2,:);
H2lin=Hlin(3:4,:);
H3lin=Hlin(5:6,:);

[mz,Pz]=MeanCov(Z,w);
Pz=Pz+blkdiag(R1,R2,R3);
Pxz=CrossCov(x,mx,Z,mz,w);


%%
Info_total=0.5*log(det(Px)/det(Px-Pxz*inv(Pz)*Pxz'))

% sequential info
for i=1:1:length(w)
    z1(i,:)=Z1(x(i,:));
end
[mz1,Pz1]=MeanCov(z1,w);
Pz1=Pz1+R1;
Pxz1=CrossCov(x,mx,z1,mz1,w);

Px1=Px-Pxz1*inv(Pz1)*Pxz1'

[x1,w1]=GH_points(mx',Px1,15);
for i=1:1:length(w1)
    z2(i,:)=Z2(x1(i,:));
end
[mz2,Pz2]=MeanCov(z2,w1);
Pz2=Pz2+R2;
Pxz2=CrossCov(x1,mx,z2,mz2,w1);

Px2=Px1-Pxz2*inv(Pz2)*Pxz2'

[x2,w2]=GH_points(mx',Px2,15);
for i=1:1:length(w2)
    z3(i,:)=Z3(x1(i,:));
end
[mz3,Pz3]=MeanCov(z3,w2);
Pz3=Pz3+R3;
Pxz3=CrossCov(x2,mx,z3,mz3,w2);

Px3=Px2-Pxz3*inv(Pz3)*Pxz3'

0.5*log(det(Px)/det(Px1))+0.5*log(det(Px1)/det(Px2))+0.5*log(det(Px2)/det(Px3))




%% Only adding variable in Pxz vs adding variables in Pxz and Pz

% only Pxz make the entries to 0
% removinf z1
PPz=Pz;
PPz([3,4],:)=PPz([3,4],:)*0;
PPz(:,[3,4])=PPz(:,[3,4])*0;
PPz([3,4],[3,4])=Pz([3,4],[3,4])

PPxz=Pxz;
PPxz(:,[3,4])=PPxz(:,[3,4])*0;
Info=0.5*log(det(Px)/det(Px-PPxz*inv(PPz)*PPxz'))
Info=0.5*log(det(Pz)/det(Pz-PPxz'*inv(Px)*PPxz))

PPz=Pz;
PPxz=Pxz;
PPz([3,4],:)=[];
PPz(:,[3,4])=[];
PPxz(:,[3,4])=[];
InfoT=0.5*log(det(Px)/det(Px-PPxz*inv(PPz)*PPxz'))
%%
V=[];
C=cell(9,19);
% disp('log det ratio, det-ratio, det-diff, tr-ratio , tr-diff , norm-ratio , norm-ratio-zver ')
C{1,1}='v1'; C{1,2}='v2'; C{1,3}='v3'; C{1,4}='log(|Px||Py|/|P|)'; C{1,5}='log(|Px|/|Px-PxyPyPxy|)'; C{1,6}='log(|Py|/|Py-PxyPxPxy|)'; C{1,7}='tr(Px)tr(Py)/tr(P)'; C{1,8}='tr(Px)/tr(Px-PxyPyPxy)'; 
C{1,9}='tr(Py)/tr(Py-PxyPxPxy)'; C{1,10}='|Px||Py|-|P|'; C{1,11}='|Px|-|Px-PxyPyPxy|'; C{1,12}='|Py|-|Py-PxyPxPxy|'; C{1,13}='tr(Px)tr(Py)-tr(P)'; C{1,14}='tr(Px)-tr(Px-PxyPyPxy)'; 
C{1,15}='tr(Py)-tr(Py-PxyPxPxy)'; C{1,16}='frob(Px)frob(Py)/frob(P)'; C{1,17}='frob(Px)/frob(Px-PxyPyPxy)';  C{1,18}='frob(Py)/frob(Py-PxyPxPxy)'; C{1,19}='frob(Px)frob(Py)-frob(P)'; 
disp('1')
kk=2;
for v1=0:1
    for v2=0:1
        for v3=0:1
            v=zeros(3,1);
            v(1)=v1;
            v(2)=v2;
            v(3)=v3;
            
            
            % True infor
            if v(1)==1 && v(2)==1 && v(3)==1
                InfoT=Px-Pxz*inv(Pz)*Pxz';
                InfoTzver=Pz-Pxz'*inv(Px)*Pxz;
                PPz=Pz;
                PPxz=Pxz;
                P=[Px,Pxz;Pxz',Pz];
            elseif v(1)==0 && v(2)==0 && v(3)==0
                InfoT=Px;
                InfoTzver=0*Pz;
                PPz=0*Pz;
                PPxz=0*Pxz;
                P=Px;
            else
                
                R=cell(3,1);
                R{1}=R1;R{2}=R2;R{3}=R3;
                RR=[];
                Zind=[];
                for i=1:1:3
                    if v(i)==1
                        Zind=horzcat(Zind,Z(:,(i-1)*hn+1:i*hn));
                        RR=blkdiag(RR,R{i});
                        
                    end
                    
                end
                %     keyboard
                [mmz,PPz]=MeanCov(Zind,w);
                PPz=PPz+RR;
                PPxz=CrossCov(x,mx,Zind,mmz,w);
                InfoT=Px-PPxz*inv(PPz)*PPxz';
                InfoTzver=PPz-PPxz'*inv(Px)*PPxz;
                P=[Px,PPxz;PPxz',PPz];
                %     InfoT=0.5*log(det(Px)/det(Px-PPxz*inv(PPz)*PPxz'));
            end
            % %FIM optimization
            % F(1)=trace(H1lin'*inv(R1)*H1lin);
            % F(2)=trace(H2lin'*inv(R2)*H2lin);
            % F(3)=trace(H3lin'*inv(R3)*H3lin);
            % vfim=bintprog(-F',[],[],[1 1 1],2)
            %
            % F(1)=trace(H1'*inv(R1)*H1);
            % F(2)=trace(H2'*inv(R2)*H2);
            % F(3)=trace(H3'*inv(R3)*H3);
            % vfim=bintprog(-F',[],[],[1 1 1],2)
            disp('new')
%             V=vertcat(V,[v1,v2,v3]);
%             C=vertcat(C,[log(det(Px)/det(InfoT));det(Px)/det(InfoT);det(Px)-det(InfoT);trace(Px)/trace(InfoT);trace(Px)-trace(InfoT);norm(Px,'fro')/norm(InfoT,'fro');norm(PPz,'fro')/norm(PPz-PPxz'*inv(Px)*PPxz,'fro')]');
            C{kk,1}=v1; C{kk,2}=v2; C{kk,3}=v3; C{kk,4}=log(det(Px)*det(PPz)/det(P)); C{kk,5}=log(det(Px)/det(InfoT)); C{kk,6}=log(det(PPz)/det(InfoTzver)); 
            C{kk,7}=trace(Px)*trace(PPz)/trace(P)'; C{kk,8}=trace(Px)/trace(InfoT); C{kk,9}=trace(PPz)/trace(InfoTzver); C{kk,10}=det(Px)*det(PPz)-det(P); C{kk,11}=det(Px)-det(InfoT); 
            C{kk,12}=det(PPz)-det(InfoTzver); C{kk,13}=trace(Px)*trace(PPz)-trace(P); C{kk,14}=trace(Px)-trace(InfoT); 
            C{kk,15}=trace(PPz)-trace(InfoTzver); C{kk,16}=norm(Px,'fro')*norm(PPz,'fro')/norm(P,'fro'); C{kk,17}=norm(Px,'fro')/norm(InfoT,'fro');  C{kk,18}=norm(PPz,'fro')/norm(InfoTzver,'fro'); 
             C{kk,19}=norm(Px,'fro')*norm(PPz,'fro')-norm(P,'fro'); 
            kk=kk+1;
        end
    end
end
 fid = fopen('test.csv', 'w') ;
 fprintf(fid, '%s,', C{1,1:end-1}) ;
 fprintf(fid, '%s\n', C{1,end}) ;
 fclose(fid) ;

 dlmwrite('test.csv', C(2:end,:), '-append') ;
C
%%
[0     0     1]
 3.2248
    4.6023
    5.0877
    4.3957

[ 0     1     0]
3.1614
    3.1037
    4.4058
    2.6607

[011]
4.6616
1.1183

[111]
6.6743
0.67805

[100]
2.3782
1.7881


[010]
3.1626
2.0921



[001]
3.3052
1.3608

%%



% Then we define the optimization problem
mpol('v',1,3)
[pd,qd]=MI_num_dem_poly_det_zver(v,Px,Pxz,Pz,hn,nz);
K=[v(1)-v(1)^2==0;v(2)-v(2)^2==0;v(3)-v(3)^2==0];
K=[K;sum(v)==2];
Prob = msdp(max(pd),mom(qd)==1,K);
[status,obj] = msol(Prob) ;
vsolz = round(double(v))


% Then we define the optimization problem
mpol('v',1,3)
[pt,qt]=MI_num_dem_poly_trace_zver(v,Px,Pxz,Pz,hn,nz);
K=[v(1)-v(1)^2==0;v(2)-v(2)^2==0;v(3)-v(3)^2==0];
K=[K;sum(v)==2];
Prob = msdp(max(pt),mom(qt)==1,K);
[status,obj] = msol(Prob) ;
vsolz = round(double(v))

% Then we define the optimization problem
mpol('v',1,3)
[pt,qt]=MI_num_dem_poly_trace_zver(v,Px,Pxz,Pz,hn,nz);
K=[v(1)-v(1)^2==0;v(2)-v(2)^2==0;v(3)-v(3)^2==0];
K=[K;sum(v)==2];
Prob = msdp(max(pt-qt),K);
[status,obj] = msol(Prob) ;
vsolz = round(double(v))


cvx_begin sdp
variable x(3) nonnegative
maximize(log_det(inv(Px)+x(1)*H1lin'*inv(R1)*H1lin+x(2)*H2lin'*inv(R2)*H2lin+x(3)*H3lin'*inv(R3)*H3lin))
subject to
sum(x)==2
x(1)<=1
x(2)<=1
x(3)<=1
cvx_end
x

cvx_begin sdp
variable y(3)
variable x(2)
maximize(log_det(inv(Px)+x(1)*H1'*inv(R1)*H1+x(2)*H2'*inv(R2)*H2+x(3)*H3'*inv(R3)*H3))
subject to
[x(1),x(2),1;1,2,1;-3,5,1]*y==[]
cvx_end
x


y=coef(qd);
[ys,ind]=sort(y);
cvx_begin
variable x(n) nonnegative
minimize(norm(ys-x,1)+norm(A*x,4))
subject to
ys>=x
cvx_end
for i=1:1:7
   yr(i)=x(ind==i) ;
end
[y,yr']

%% checking if det polys are same as info
pnum=@(v)3.8481-3.0792*v(1)^2*v(2)^2-2.5483*v(1)^2*v(3)^2-3.3745*v(2)^2*v(3)^2+4.5567*v(1)^2*v(2)^2*v(3)^2+0.10998*v(1)^4*v(2)^4+0.45131*v(1)^4*v(2)^2*v(3)^2+0.23587*v(1)^4*v(3)^4+0.15488*v(1)^2*v(2)^4*v(3)^2+0.47874*v(1)^2*v(2)^2*v(3)^4+0.02345*v(2)^4*v(3)^4-0.19914*v(1)^4*v(2)^4*v(3)^2-0.57053*v(1)^4*v(2)^2*v(3)^4-0.073337*v(1)^2*v(2)^4*v(3)^4+0.050292*v(1)^4*v(2)^4*v(3)^4;
qden=@(v)3.8481-3.9391*v(1)^2-3.9673*v(2)^2-3.6422*v(3)^2+0.75396*v(1)^4+3.998*v(1)^2*v(2)^2+3.4444*v(1)^2*v(3)^2+0.26456*v(2)^4+3.7225*v(2)^2*v(3)^2+0.074413*v(3)^4-0.73156*v(1)^4*v(2)^2-0.52837*v(1)^4*v(3)^2-0.25929*v(1)^2*v(2)^4-3.4405*v(1)^2*v(2)^2*v(3)^2+0.1771*v(1)^2*v(3)^4-0.24074*v(2)^4*v(3)^2-0.0471*v(2)^2*v(3)^4+0.043181*v(1)^4*v(2)^4+0.4823*v(1)^4*v(2)^2*v(3)^2-0.15023*v(1)^4*v(3)^4+0.21137*v(1)^2*v(2)^4*v(3)^2-0.22987*v(1)^2*v(2)^2*v(3)^4-0.0037846*v(2)^4*v(3)^4-0.022855*v(1)^4*v(2)^4*v(3)^2+0.17299*v(1)^4*v(2)^2*v(3)^4+0.025717*v(1)^2*v(2)^4*v(3)^4-0.015316*v(1)^4*v(2)^4*v(3)^4;


v=[1,1,0];
log(pnum(v)/qden(v))
4.4224

v=[1,0,1];
log(pnum(v)/qden(v))
3.6972


v=[0,1,1];
log(pnum(v)/qden(v))
4.0747



v=[1,1,1];
log(pnum(v)/qden(v))
5.1436


v=[1,0,0];
log(pnum(v)/qden(v))
1.7586


v=[0,1,0];
log(pnum(v)/qden(v))
3.2761


v=[0,0,1];
log(pnum(v)/qden(v))
2.6194


%% checking if tracve polys are same as info
pnum=@(v)4.2802*v(1)+27.6705*v(2)+4.2835*v(3); %c=[4.2802,27.6705,4.2835]
qden=@(v)1.7013*v(1)+25.102*v(2)+0.69163*v(3); %e=[1.7013,25.102,0.69163]

%cvx linear frac opt



%%
% plotting
[xx,yy]=meshgrid(-30:1:30);
pp=zeros(size(xx));

for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        pp(i,j)=mvnpdf([xx(i,j),yy(i,j)],mx,Px);
    end
end
figure
contour(xx,yy,pp,30);
hold on
plot(mx(1),mx(2),'ro',2,2,'ks',-2,-2,'ks',10,0,'ks')
plot(x(:,1),x(:,2),'b+')


%%
MCx=mvnrnd(mx,Px,10000);
MCz=zeros(size(MCx,1),2);
for i=1:1:size(MCx,1)
    MCz(i,:)=Z1(MCx(i,:));
end
mean(MCz,1)


%%
%    Multiple targets
Px1=[20,-10;-10,45];
mx1=[3,5];
Px2=[30,15;15,10];
mx2=[-3,5];
a1=2;
b1=2;
a2=-2;
b2=-2;
a3=10;
b3=0;
Z1=@(x)[norm([x(1)-a1,x(2)-b1]);atan2(x(2)-b1,x(1)-a1)];
Z2=@(x)[norm([x(1)-a2,x(2)-b2]);atan2(x(2)-b2,x(1)-a2)];
Z3=@(x)[norm([x(1)-a3,x(2)-b3]);atan2(x(2)-b3,x(1)-a3)];
hjack=@(x,a,b)[(x(1)-a)/norm([x(1)-a,x(2)-b]),(x(2)-b)/norm([x(1)-a,x(2)-b]);(b-x(2))/norm([x(1)-a,x(2)-b])^2,(-a+x(1))/norm([x(1)-a,x(2)-b])^2];
R1=diag([1^2,(1*pi/180)^2]);
R2=diag([5^2,(0.1*pi/180)^2]);
R3=diag([0.5^2,(10*pi/180)^2]);
hn=2;
fn=2;
nz=3;
nx=2;

[x1,w1]=GH_points(mx1',Px1,5);
[x2,w2]=GH_points(mx2',Px2,5);
z1_x1=zeros(length(w1),2);
z1_x2=zeros(length(w1),2);
z2_x1=zeros(length(w1),2);
z2_x2=zeros(length(w1),2);
z3_x1=zeros(length(w1),2);
z3_x2=zeros(length(w1),2);

for i=1:1:length(w1)
    z1_x1(i,:)=Z1(x1(i,:));
    z1_x2(i,:)=Z1(x2(i,:));
    
    z2_x1(i,:)=Z2(x1(i,:));
    z2_x2(i,:)=Z2(x2(i,:));
    
    z3_x1(i,:)=Z3(x1(i,:));
    z3_x2(i,:)=Z3(x2(i,:));
end

Z1=[z1_x1,z2_x1,z3_x1];
Z2=[z1_x2,z2_x2,z3_x2];

Hlin_x1=StatLinearize_sigmapts(Z1,x1,w1);
Hlin_x2=StatLinearize_sigmapts(Z2,x2,w2);

H1lin_x1=Hlin_x1(1:2,:);
H2lin_x1=Hlin_x1(3:4,:);
H3lin_x1=Hlin_x1(5:6,:);

H1lin_x2=Hlin_x2(1:2,:);
H2lin_x2=Hlin_x2(3:4,:);
H3lin_x2=Hlin_x2(5:6,:);

[mz1,Pz1]=MeanCo*v(Z1,w1);
Pz1=Pz1+blkdiag(R1,R2,R3);
Pxz1=CrossCov(x1,mx1,Z1,mz1,w1);

[mz2,Pz2]=MeanCov(Z2,w2);
Pz2=Pz2+blkdiag(R1,R2,R3);
Pxz2=CrossCov(x2,mx2,Z2,mz2,w2);

%% FIM optimization
H1_x1=hjack(mx1,a1,b1);
H2_x1=hjack(mx1,a2,b2);
H3_x1=hjack(mx1,a3,b3);

H1_x2=hjack(mx2,a1,b1);
H2_x2=hjack(mx2,a2,b2);
H3_x2=hjack(mx2,a3,b3);

% F(1,1)=trace(H1lin_x1'*inv(R1)*H1lin_x1);
% F(1,2)=trace(H2lin_x1'*inv(R2)*H2lin_x1);
% F(1,3)=trace(H3lin_x1'*inv(R3)*H3lin_x1);
%
% F(2,1)=trace(H1lin_x2'*inv(R1)*H1lin_x2);
% F(2,2)=trace(H2lin_x2'*inv(R2)*H2lin_x2);
% F(2,3)=trace(H3lin_x2'*inv(R3)*H3lin_x2);

F(1,1)=trace(H1_x1'*inv(R1)*H1_x1);
F(1,2)=trace(H2_x1'*inv(R2)*H2_x1);
F(1,3)=trace(H3_x1'*inv(R3)*H3_x1);

F(2,1)=trace(H1_x2'*inv(R1)*H1_x2);
F(2,2)=trace(H2_x2'*inv(R2)*H2_x2);
F(2,3)=trace(H3_x2'*inv(R3)*H3_x2);

vv=bintprog(-[F(1,:)';F(2,:)'],[],[],[eye(3),eye(3)],ones(3,1));
vfim=[vv(1:3)';vv(4:6)']

cvx_begin
variable mmu(2) binary
variable u(2,3) binary
maximize(mmu(1)*trace(inv(Px1))+u(1,1)*F(1,1)+u(1,2)*F(1,2)+u(1,3)*F(1,3)...
        +mmu(2)*trace(inv(Px2))+u(2,1)*F(2,1)+u(2,2)*F(2,2)+u(2,3)*F(2,3) )
subject to
 u(1,1)+u(2,1)==1
 u(1,2)+u(2,2)==1
 u(1,3)+u(2,3)==1
 1/3*(u(1,1)+u(1,2)+u(1,3))<=mmu(1)<=(u(1,1)+u(1,2)+u(1,3))
 1/3*(u(2,1)+u(2,2)+u(2,3))<=mmu(2)<=(u(2,1)+u(2,2)+u(2,3))
cvx_end

%%  Assume MI is sum of individual i.e. maximixe the bound
% eg I(x1:z1,z2,z3)=I(x1:z1)+I(x1:z2)+I(x1:z3)
tic
MI=zeros(2,3);

ZZ1=[z1_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R1);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT11=Px1-PPxz1*inv(PPz1)*PPxz1';
MI(1,1)=log(det(Px1)/det(InfoT11));

ZZ1=[z2_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R2);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT12=Px1-PPxz1*inv(PPz1)*PPxz1';
MI(1,2)=log(det(Px1)/det(InfoT12));

ZZ1=[z3_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R3);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT13=Px1-PPxz1*inv(PPz1)*PPxz1';
MI(1,3)=log(det(Px1)/det(InfoT13));

ZZ2=[z1_x2]
[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+blkdiag(R1);
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);
InfoT21=Px2-PPxz2*inv(PPz2)*PPxz2';
MI(2,1)=log(det(Px2)/det(InfoT21));

ZZ2=[z2_x2]
[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+blkdiag(R2);
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);
InfoT22=Px2-PPxz2*inv(PPz2)*PPxz2';
MI(2,2)=log(det(Px2)/det(InfoT22));

ZZ2=[z3_x2]
[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+blkdiag(R3);
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);
InfoT23=Px2-PPxz2*inv(PPz2)*PPxz2';
MI(2,3)=log(det(Px2)/det(InfoT23));

vv=bintprog(-[MI(1,:)';MI(2,:)'],[],[],[eye(3),eye(3)],ones(3,1));
vMI=[vv(1:3)';vv(4:6)']
toc
%% x1 z1 z2 z3
ZZ1=[z1_x1,z2_x1,z3_x1];

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R1,R2,R3);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2;
clc
log(det(Px1)/det(InfoT1))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)
  3.7373
   39.0318
   [trace(InfoT1),trace(Px2)]
      [det(InfoT1),det(Px2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [       19.055       11.536          800           65
%           800           65          800           65]

%% x2 z1 z2 z3
ZZ2=[z1_x2,z2_x2,z3_x2];

[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+blkdiag(R1,R2,R3);
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1;
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';
clc
log(det(Px2)/det(InfoT2))
trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
3.3676
21.4029
   [trace(Px1),trace(InfoT2)]
      [det(Px1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [800           65       19.975       13.592
%  800           65          800           65]

%% x1 z1 z2 , x2 z3
ZZ1=[z1_x1,z2_x1];
ZZ2=[z3_x2];

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R1,R2);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);


[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+blkdiag(R3);
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';
clc
log(det(Px1)/det(InfoT1))+log(det(Px2)/det(InfoT2))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)+trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
4.2316
31.5397
   [trace(InfoT1),trace(InfoT2)]
      [det(InfoT1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [120.04       22.388       45.496       25.119
%  800           65          800           65]
%% x1 z1 z3 , x2 z2
ZZ1=[z1_x1,z3_x1];
ZZ2=[z2_x2];
RR1=blkdiag(R1,R3);
RR2=blkdiag(R2);

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+RR1;
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);

[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+RR2;
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';
clc
log(det(Px1)/det(InfoT1))+log(det(Px2)/det(InfoT2))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)+trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
 3.7728
 28.0290
    [trace(InfoT1),trace(InfoT2)]
      [det(InfoT1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [28.637       13.143       261.65       33.369
%  800           65          800           65]
%% x1 z2 z3 , x2 z1
ZZ1=[z2_x1,z3_x1];
ZZ2=[z1_x2];
RR1=blkdiag(R2,R3);
RR2=blkdiag(R1);

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+RR1;
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);

[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+RR2;
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';
clc
log(det(Px1)/det(InfoT1))+log(det(Px2)/det(InfoT2))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)+trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
3.6991
39.8977
    [trace(InfoT1),trace(InfoT2)]
      [det(InfoT1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [       44.436       16.546       132.02       30.473
%           800           65          800           65]

%% x1 z1 , x2 z2 Z3
ZZ1=[z1_x1];
ZZ2=[z2_x2,z3_x2];
RR1=blkdiag(R1);
RR2=blkdiag(R2,R3);

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+RR1;
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);

[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+RR2;
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';
clc
log(det(Px1)/det(InfoT1))+log(det(Px2)/det(InfoT2))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)+trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
4.887
20.5370
    [trace(InfoT1),trace(InfoT2)]
      [det(InfoT1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [       221.67       29.947       33.223       19.328
%           800           65          800           65]

%% x1 z2 , x2 z1 Z3
ZZ1=[z2_x1];
ZZ2=[z1_x2,z3_x2];
RR1=blkdiag(R2);
RR2=blkdiag(R1,R3);

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+RR1;
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);

[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+RR2;
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';
clc
log(det(Px1)/det(InfoT1))+log(det(Px2)/det(InfoT2))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)+trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
4.8962
32.4056
    [trace(InfoT1),trace(InfoT2)]
      [det(InfoT1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [       154.03       25.026       21.825       14.565
%           800           65          800           65]

%% x1 z3 , x2 z1 Z2

ZZ1=[z3_x1];
ZZ2=[z1_x2,z2_x2];
RR1=blkdiag(R3);
RR2=blkdiag(R1,R2);

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+RR1;
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);

[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+RR2;
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';

clc
log(det(Px1)/det(InfoT1))+log(det(Px2)/det(InfoT2))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)+trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
3.5786
28.8949
    [trace(InfoT1),trace(InfoT2)]
      [det(InfoT1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [       82.686       23.516        101.5       25.149
%           800           65          800           65]


%%

Px=cell(2,1);
Pxz=cell(2,1);
Px{1}=Px1;  Px{2}=Px2;
Pxz{1}=Pxz1;  Pxz{2}=Pxz2;
Pz=cell(2,1);
Pz{1}=Pz1;  Pz{2}=Pz2;
mpol('v',nx,nz)
[p,q]=MI_num_dem_poly_zver_multi_x(v,Px,Pxz,Pz,hn,nz,nx,0)

K=(ones(nx,nz)-v).*v;
K=reshape(K,nx*nz,1);
mpol('S',nz,1);
for i=1:1:nz
    S(i)=sum(v(:,i));
end
K=[K==0;S==1];

Prob = msdp(max(q),K);
tic

[status,obj] = msol(Prob)
toc
vsolz = round(double(v))

%% is MI a sum of individual MI

% ZZ1=[z1_x1,z2_x1,z3_x1];3.7373
ZZ1=[z1_x1,z2_x1,z3_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R1,R2,R3);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';

ZZ1=[z1_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R1);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT11=Px1-PPxz1*inv(PPz1)*PPxz1';

ZZ1=[z2_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);

PPz1=PPz1+blkdiag(R2);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT12=Px1-PPxz1*inv(PPz1)*PPxz1';

ZZ1=[z3_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R3);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT13=Px1-PPxz1*inv(PPz1)*PPxz1';


log(det(Px1)/det(InfoT1))
log(det(Px1)/det(InfoT11))+log(det(Px1)/det(InfoT12))+log(det(Px1)/det(InfoT13))
%% matching the z-ver trace and det
% v=sym('v',[nx,nz]);
PPxz=Pxz2;
PPz=Pz2;
Px=Px2;
v=[1 0 0];

for i=1:1:nz
    PPxz(:,(i-1)*hn+1:i*hn)=v(i)*PPxz(:,(i-1)*hn+1:i*hn);
    for j=1:1:nz
%         if j~=i
            PPz((i-1)*hn+1:i*hn,(j-1)*hn+1:j*hn)=v(i)*PPz((i-1)*hn+1:i*hn,(j-1)*hn+1:j*hn);
            PPz((j-1)*hn+1:j*hn,(i-1)*hn+1:i*hn)=v(i)*PPz((j-1)*hn+1:j*hn,(i-1)*hn+1:i*hn);
%         end
    end
end

q=trace(PPz-PPxz'*inv(Px)*PPxz);
p=trace(PPz);

PPxz=Pxz2;
PPz=Pz2;
Px=Px2;
PPxz_red=[];
PPz_red=[];
for i=1:1:nz
    if v(i)==1
        PPxz_red=[PPxz_red,PPxz(:,(i-1)*hn+1:i*hn)];
    end
    PP=[];
    for j=1:1:nz
            if v(i)==1 && v(j)==1
                PP=[PP,PPz((i-1)*hn+1:i*hn,(j-1)*hn+1:j*hn)];
            end
    end
    PPz_red=[PPz_red;PP];
end
q_red=trace(PPz_red-PPxz_red'*inv(Px)*PPxz_red);
p_red=trace(PPz_red);

[q_red,q;
    p_red,p]

%% Plotting the targets
[xx,yy]=meshgrid(-20:0.1:20);
pp1=zeros(size(xx));
pp2=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        pp1(i,j)=mvnpdf([xx(i,j),yy(i,j)],mx1,Px1);
        pp2(i,j)=mvnpdf([xx(i,j),yy(i,j)],mx2,Px2);
     end
end
contour(xx,yy,pp1)
hold on
contour(xx,yy,pp2)
plot(2,2,'ks',-2,-2,'ks',10,0,'ks','linewidth',4,'MarkerSize',8)

%% Plotting the targets
Px11_MI4= [5.7102    4.1758
    4.1758   10.8358];
Px22_MI4= [ 21.1347   11.0866
   11.0866    7.3967];

Px11_MUB=[11.4930    2.2765
    2.2765   10.8952];
Px22_MUB=[ 4.5393    3.7835
    3.7835    4.7533];


[xx,yy]=meshgrid(-20:0.1:20);
pp1MUB=zeros(size(xx));
pp2MUB=zeros(size(xx));
pp1MI4=zeros(size(xx));
pp2MI4=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        pp1MUB(i,j)=mvnpdf([xx(i,j),yy(i,j)],[0,0],Px11_MUB);
        pp2MUB(i,j)=mvnpdf([xx(i,j),yy(i,j)],[0,0],Px22_MUB);
        pp1MI4(i,j)=mvnpdf([xx(i,j),yy(i,j)],[0,0],Px11_MI4);
        pp2MI4(i,j)=mvnpdf([xx(i,j),yy(i,j)],[0,0],Px22_MI4);
     end
end
figure(1)
contour(xx,yy,pp1MUB,1e-3,'r','linewidth',2)
hold on
contour(xx,yy,pp1MI4,1e-3,'b.-','linewidth',2)
legend('MI_{UB}','MI_4')
xlabel('x')
ylabel('y')
plot_prop_paper
grid


figure(2)
contour(xx,yy,pp2MUB,1e-3,'r','linewidth',2)
hold on
contour(xx,yy,pp2MI4,1e-3,'b.-','linewidth',2)
legend('MI_{UB}','MI_4')
xlabel('x')
ylabel('y')
plot_prop_paper
grid

