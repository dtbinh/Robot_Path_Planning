clc
clear
close all
% P=magic(5)*magic(5)'/100;
Px1=[3^2,5;5,3^2];
mx1=[3,5];
Px2=[5^2,-10;-10,3^2];
mx2=[7,5];

[xx,yy]=meshgrid(-20:0.2:20);
pp1_prr=zeros(size(xx));
pp1_pst=zeros(size(xx));

pp2_prr=zeros(size(xx));
pp2_pst=zeros(size(xx));


s1=[0,0];
s2=[5,0];
s3=[10,0];
Z1=@(x)[norm([x(1)-s1(1),x(2)-s1(2)]);atan((x(2)-s1(2))/(x(1)-s1(1)))];
Z2=@(x)[norm([x(1)-s2(1),x(2)-s2(2)]);atan((x(2)-s2(2))/(x(1)-s2(1)))];
Z3=@(x)[norm([x(1)-s3(1),x(2)-s3(2)]);atan((x(2)-s3(2))/(x(1)-s3(1)))];
hjack=@(x,s)[(x(1)-s(1))/norm([x(1)-s(1),x(2)-s(2)]),(x(2)-s(2))/norm([x(1)-s(1),x(2)-s(2)]);...
    (s(2)-x(2))/norm([x(1)-s(1),x(2)-s(2)])^2,(-s(1)+x(1))/norm([x(1)-s(1),x(2)-s(2)])^2];

R1=diag([10^2,(10*pi/180)^2]);
R2=diag([5^2,(5*pi/180)^2]);
R3=diag([10^2,(10*pi/180)^2]);
hn=2;
nz=3;
H_x1=[hjack(mx1,s1);hjack(mx1,s2);hjack(mx1,s3)];
H_x2=[hjack(mx2,s1);hjack(mx2,s2);hjack(mx2,s3)];




% [x,w]=UT_sigmapoints(mx',Px,2);
[x1,w_x1]=GH_points(mx1',Px1,5);
[x2,w_x2]=GH_points(mx2',Px2,5);
z1_x1=zeros(length(w_x1),hn);
z2_x1=zeros(length(w_x1),hn);
z3_x1=zeros(length(w_x1),hn);

z1_x2=zeros(length(w_x2),hn);
z2_x2=zeros(length(w_x2),hn);
z3_x2=zeros(length(w_x2),hn);

for i=1:1:length(w_x1)
    z1_x1(i,:)=Z1(x1(i,:));
    z2_x1(i,:)=Z2(x1(i,:));
    z3_x1(i,:)=Z3(x1(i,:));
    
    z1_x2(i,:)=Z1(x2(i,:));
    z2_x2(i,:)=Z2(x2(i,:));
    z3_x2(i,:)=Z3(x2(i,:));
end

Z_x1=[z1_x1,z2_x1,z3_x1];
Z_x2=[z1_x2,z2_x2,z3_x2];



[mz_x1,Pz_x1]=MeanCov(Z_x1,w_x1);
[mz_x2,Pz_x2]=MeanCov(Z_x2,w_x2);

Pz_x1=Pz_x1+blkdiag(R1,R2,R3);
Pz_x2=Pz_x2+blkdiag(R1,R2,R3);

Pxz_x1=CrossCov(x1,mx1,Z_x1,mz_x1,w_x1);
Pxz_x2=CrossCov(x2,mx2,Z_x2,mz_x2,w_x2);




%%
V=[];
C=cell(9,27);
% disp('log det ratio, det-ratio, det-diff, tr-ratio , tr-diff , norm-ratio , norm-ratio-zver ')
cols={'vz1_x1','vz2_x1','vz3_x1','vz1_x2','vz2_x2','vz3_x2','log(|Px||Py|/|P|)','log(|Px|/|Px-PxyPyPxy|)','log(|Py|/|Py-PxyPxPxy|)','tr(Px)tr(Py)/tr(P)','tr(Px)/tr(Px-PxyPyPxy)','tr(Py)/tr(Py-PxyPxPxy)',...
    '|Px||Py|-|P|','|Px|-|Px-PxyPyPxy|','|Py|-|Py-PxyPxPxy|','tr(Px)tr(Py)-tr(P)','tr(Px)-tr(Px-PxyPyPxy)','tr(Py)-tr(Py-PxyPxPxy)','frob(Px)frob(Py)/frob(P)',...
    'frob(Px)/frob(Px-PxyPyPxy)','frob(Py)/frob(Py-PxyPxPxy)','frob(Px)frob(Py)-frob(P)','sumTrPx*sumTrPy/sumTrP','(trPx1*trPz1+trPx2*trPz2)/(trP1+trP2)' ,'sumTrPy/sumTrPypost','sumfroPx*sumfroPx/sumfroP',...
    'sum(froPx*froPy)/sum(froP)','sum(froPx)/sum(froPxpost)','sum(froPy)/sum(froPypost)','detFIM','trFIM','MIUB'};
for i=1:1:length(cols)
    C{1,i}=cols{i};
end
kk=2;
for vz1_x1=0:1
    for vz1_x2=0:1
        for vz2_x1=0:1
            for vz2_x2=0:1
                for vz3_x1=0:1
                    for vz3_x2=0:1
                        if (vz1_x1+vz1_x2)~=1 || (vz2_x1+vz2_x2)~=1 || (vz3_x1+vz3_x2)~=1
                            continue
                        end
                        %FIM
                        FIM=inv(blkdiag(Px1,Px2))+[H_x1,zeros(6,2);zeros(6,2),H_x2]'*blkdiag(vz1_x1*inv(R1),vz2_x1*inv(R2),vz3_x1*inv(R3),vz1_x2*inv(R1),vz2_x2*inv(R2),vz3_x2*inv(R3))*[H_x1,zeros(6,2);zeros(6,2),H_x2];
                        PP=inv(FIM);
                        figure(kk)
                        clf
                        subplot(1,2,1)
                        
                        plot_1sigellip([0,0],Px1,'r',2.5);
                        hold on
                        plot_1sigellip([0,0],PP(1:2,1:2),'b--',2.5);
                        title(strcat('Target 1 : ',num2str([vz1_x1,vz2_x1,vz3_x1])))
                        hold off
                        subplot(1,2,2)
                        
                        plot_1sigellip([0,0],Px2,'r',2.5);
                        hold on
                        plot_1sigellip([0,0],PP(3:4,3:4),'b--',2.5);
                        title(strcat('Target 2 : ',num2str([vz1_x2,vz2_x2,vz3_x2])))
                        hold off
                        pause(0.2)
                        ss=num2str([vz1_x1,vz2_x1,vz3_x1,vz1_x2,vz2_x2,vz3_x2]);
                        ss(find(ss==' '))='_';
                        saveas(gcf, strcat('COVellip_EKF_',ss), 'png')
                        pause(0.3)
                        
                        % True infor
                        if vz1_x1==0 && vz2_x1==0 && vz3_x1==0
                            
                            
                            PPz_x1=Pz_x1;
                            PPxz_x1=0*Pxz_x1;
                            
                            InfoT_x1=Px1-PPxz_x1*inv(PPz_x1)*PPxz_x1';
                            InfoTzver_x1=PPz_x1-PPxz_x1'*inv(Px1)*PPxz_x1;
                            P1=[Px1,PPxz_x1;PPxz_x1',PPz_x1];
                            
                            
                        else
                            
                            
                            R=cell(3,1);
                            R{1}=R1;R{2}=R2;R{3}=R3;
                            RR=[];
                            Zind=[];
                            
                            if vz1_x1==1
                                Zind=horzcat(Zind,z1_x1);
                                RR=blkdiag(RR,R1);
                                
                            end
                            if vz2_x1==1
                                Zind=horzcat(Zind,z2_x1);
                                RR=blkdiag(RR,R2);
                                
                            end
                            if vz3_x1==1
                                Zind=horzcat(Zind,z3_x1);
                                RR=blkdiag(RR,R3);
                                
                            end
                            %     keyboard
                            [mmz_x1,PPz_x1]=MeanCov(Zind,w_x1);
                            PPz_x1=PPz_x1+RR;
                            PPxz_x1=CrossCov(x1,mx1,Zind,mmz_x1,w_x1);
                            InfoT_x1=Px1-PPxz_x1*inv(PPz_x1)*PPxz_x1';
                            InfoTzver_x1=PPz_x1-PPxz_x1'*inv(Px1)*PPxz_x1;
                            P1=[Px1,PPxz_x1;PPxz_x1',PPz_x1];
                            %     InfoT=0.5*log(det(Px)/det(Px-PPxz*inv(PPz)*PPxz'));
                        end
                        
                        % Now for target 2
                        if vz1_x2==0 && vz2_x2==0 && vz3_x2==0
                            
                            
                            
                            PPz_x2=Pz_x2;
                            PPxz_x2=0*Pxz_x2;
                            
                            
                            InfoT_x2=Px2-PPxz_x2*inv(PPz_x2)*PPxz_x2';
                            InfoTzver_x2=PPz_x2-PPxz_x2'*inv(Px2)*PPxz_x2;
                            P2=[Px2,PPxz_x2;PPxz_x2',PPz_x2];
                            
                            
                        else
                            
                            
                            R=cell(3,1);
                            R{1}=R1;R{2}=R2;R{3}=R3;
                            RR=[];
                            Zind=[];
                            
                            if vz1_x2==1
                                Zind=horzcat(Zind,z1_x2);
                                RR=blkdiag(RR,R1);
                                
                            end
                            if vz2_x2==1
                                Zind=horzcat(Zind,z2_x2);
                                RR=blkdiag(RR,R2);
                                
                            end
                            if vz3_x2==1
                                Zind=horzcat(Zind,z3_x2);
                                RR=blkdiag(RR,R3);
                                
                            end
                            %     keyboard
                            [mmz_x2,PPz_x2]=MeanCov(Zind,w_x2);
                            PPz_x2=PPz_x2+RR;
                            PPxz_x2=CrossCov(x2,mx2,Zind,mmz_x2,w_x2);
                            InfoT_x2=Px2-PPxz_x2*inv(PPz_x2)*PPxz_x2';
                            InfoTzver_x2=PPz_x2-PPxz_x2'*inv(Px2)*PPxz_x2;
                            P2=[Px2,PPxz_x2;PPxz_x2',PPz_x2];
                            %     InfoT=0.5*log(det(Px)/det(Px-PPxz*inv(PPz)*PPxz'));
                        end
                        
                        
                        C{kk,1}=vz1_x1; C{kk,2}=vz2_x1; C{kk,3}=vz3_x1; C{kk,4}=vz1_x2; C{kk,5}=vz2_x2; C{kk,6}=vz3_x2;
                        C{kk,7}=log(det(Px1)*det(PPz_x1)/det(P1))+log(det(Px2)*det(PPz_x2)/det(P2));
                        C{kk,8}=log(det(Px1)/det(InfoT_x1))+log(det(Px2)/det(InfoT_x2));
                        C{kk,9}=log(det(PPz_x1)/det(InfoTzver_x1))+log(det(PPz_x2)/det(InfoTzver_x2));
                        C{kk,10}=trace(Px1)*trace(PPz_x1)/trace(P1)+trace(Px2)*trace(PPz_x2)/trace(P2);
                        C{kk,11}=trace(Px1)/trace(InfoT_x1)+trace(Px2)/trace(InfoT_x2);
                        C{kk,12}=trace(PPz_x1)/trace(InfoTzver_x1)+trace(PPz_x2)/trace(InfoTzver_x2);
                        C{kk,13}=det(Px1)*det(PPz_x1)-det(P1)+det(Px2)*det(PPz_x2)-det(P2);
                        C{kk,14}=det(Px1)-det(InfoT_x1)+det(Px2)-det(InfoT_x2);
                        C{kk,15}=det(PPz_x1)-det(InfoTzver_x1)+det(PPz_x2)-det(InfoTzver_x2);
                        C{kk,16}=trace(Px1)*trace(PPz_x1)-trace(P1)+trace(Px2)*trace(PPz_x2)-trace(P2);
                        C{kk,17}=trace(Px1)-trace(InfoT_x1)+trace(Px2)-trace(InfoT_x2);
                        C{kk,18}=trace(PPz_x1)-trace(InfoTzver_x1)+trace(PPz_x2)-trace(InfoTzver_x2);
                        C{kk,19}=norm(Px1,'fro')^2*norm(PPz_x1,'fro')^2/norm(P1,'fro')^2+norm(Px2,'fro')^2*norm(PPz_x2,'fro')^2/norm(P2,'fro')^2;
                        C{kk,20}=norm(Px1,'fro')^2/norm(InfoT_x1,'fro')^2+norm(Px2,'fro')^2/norm(InfoT_x2,'fro')^2;
                        C{kk,21}=norm(PPz_x1,'fro')^2/norm(InfoTzver_x1,'fro')^2+norm(PPz_x2,'fro')^2/norm(InfoTzver_x2,'fro')^2;
                        C{kk,22}=norm(Px1,'fro')^2*norm(PPz_x1,'fro')^2-norm(P1,'fro')^2+norm(Px2,'fro')^2*norm(PPz_x2,'fro')^2-norm(P2,'fro')^2;
                        
                        C{kk,23}=(trace(Px1)+trace(Px2))*(trace(PPz_x1)+trace(PPz_x2))/(trace(P1)+trace(P2));
                        C{kk,24}=(trace(Px1)*trace(PPz_x1)+trace(Px2)*trace(PPz_x2))/(trace(P1)+trace(P2));
                        
                        C{kk,25}=(trace(PPz_x1)+trace(PPz_x2))/(trace(InfoTzver_x1)+trace(InfoTzver_x2));
                        
                        C{kk,26}=(norm(Px1,'fro')^2+norm(Px2,'fro')^2)*(norm(PPz_x1,'fro')^2+norm(PPz_x2,'fro')^2)/(norm(P1,'fro')^2+norm(P2,'fro')^2) ;
                        C{kk,27}=(norm(Px1,'fro')^2*norm(PPz_x1,'fro')^2+norm(Px2,'fro')^2*norm(PPz_x2,'fro')^2)/(norm(P1,'fro')^2+norm(P2,'fro')^2);
                        C{kk,28}=(norm(Px1,'fro')^2+norm(Px2,'fro')^2)/(norm(InfoT_x1,'fro')^2+norm(InfoT_x2,'fro')^2) ;
                        C{kk,29}=(norm(PPz_x1,'fro')^2+norm(PPz_x2,'fro')^2)/(norm(InfoTzver_x1,'fro')^2+norm(InfoTzver_x2,'fro')^2);
                        
                        
                        C{kk,30}=det(FIM);
                        C{kk,31}=trace(FIM);
                        
                        
                        % plotting the covariances
                        
                        %             for i=1:1:size(xx,1)
                        %                 for j=1:1:size(xx,2)
                        %                     pp1_prr(i,j)=mvnpdf([xx(i,j),yy(i,j)],mx1,Px1);
                        %                     pp1_pst(i,j)=mvnpdf([xx(i,j),yy(i,j)],mx1,Px1);
                        %
                        %                     pp2_prr(i,j)=mvnpdf([xx(i,j),yy(i,j)],mx1,Px1);
                        %                     pp2_pst(i,j)=mvnpdf([xx(i,j),yy(i,j)],mx1,Px1);
                        %                 end
                        %             end
                        figure(kk)
                        
                        subplot(1,2,1)
                        
                        plot_1sigellip([0,0],Px1,'r',2.5);
                        hold on
                        plot_1sigellip([0,0],InfoT_x1,'b--',2.5);
                        title(strcat('Target 1 : ',num2str([vz1_x1,vz2_x1,vz3_x1])))
                        hold off
                        subplot(1,2,2)
                        
                        plot_1sigellip([0,0],Px2,'r',2.5);
                        hold on
                        plot_1sigellip([0,0],InfoT_x2,'b--',2.5);
                        title(strcat('Target 2 : ',num2str([vz1_x2,vz2_x2,vz3_x2])))
                        hold off
                        pause(0.2)
                        ss=num2str([vz1_x1,vz2_x1,vz3_x1,vz1_x2,vz2_x2,vz3_x2]);
                        ss(find(ss==' '))='_';
                        saveas(gcf, strcat('COVellip_Quad_',ss), 'png')
                        pause(0.3)
                        
                        
                        kk=kk+1;
                    end
                end
            end
        end
    end
end

% fid = fopen('ConfEg22.csv', 'w') ;
% fprintf(fid, '%s,', C{1,1:end-1}) ;
% fprintf(fid, '%s\n', C{1,end}) ;
% fclose(fid) ;
% 
% dlmwrite('ConfEg22.csv', C(2:end,:), '-append') ;
% C


%%



% Then we define the optimization problem
mpol('v',1,4)
vindfun=@(ns,nt)(nt-1)*nz+ns;
[pd,qd]=MI_num_dem_poly_det_zver(v,{Px1,Px2},{Pxz_x1,Pxz_x2},{Pz_x1,Pz_x2},hn,nz,vindfun);
K=[v(1)-v(1)^2==0;v(2)-v(2)^2==0;v(3)-v(3)^2==0;v(4)-v(4)^2==0];
K=[K;v(1)+v(3)==1;v(2)+v(4)==1];
Prob = msdp(max(pd),mom(qd)==1,K);
[status,obj] = msol(Prob) ;
vsolz = round(double(v))


% Then we define the optimization problem
mpol('v',1,4)
[pt,qt]=MI_num_dem_poly_trace_zver(v,{Px1,Px2},{Pxz_x1,Pxz_x2},{Pz_x1,Pz_x2},hn,nz,vindfun);
K=[v(1)-v(1)^2==0;v(2)-v(2)^2==0;v(3)-v(3)^2==0;v(4)-v(4)^2==0];
K=[K;v(1)+v(3)==1;v(2)+v(4)==1];
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
pp1=zeros(size(xx));
pp2=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        pp1(i,j)=mvnpdf([xx(i,j),yy(i,j)],mx1,Px1);
        pp2(i,j)=mvnpdf([xx(i,j),yy(i,j)],mx2,Px2);
        
    end
end
figure
contour(xx,yy,pp1,10);
hold on
contour(xx,yy,pp2,10);

plot(mx1(1),mx1(2),'ro',mx2(1),mx2(2),'rx',s1(1),s1(2),'ks',s2(1),s2(2),'ks')



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

