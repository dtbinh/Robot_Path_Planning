%% Using contextual information
%  by psedo measurements

%% first generate the arena
figure
hold on
axis([0 10 0 10])
[xu,yu] = ginput(7);
plot(xu,yu)
axis([0 10 0 10])
[xd,yd] = ginput(7);
plot(xd,yd)
axis([0 10 0 10])
[xcu,ycu] = ginput(7);
plot(xcu,ycu)
axis([0 10 0 10])
[xcd,ycd] = ginput(7);
plot(xcd,ycd)
axis([0 10 0 10])

save('EgNaval','xu','yu','xd','yd','xcu','ycu','xcd','ycd')
%%
load('EgNaval')
plot(xu,yu,'k',xd,yd,'k',xcu,ycu,'k',xcd,ycd,'k')
hold on
dx=0.5;
Xc=0:dx:10;
Yc=0:dx:10;


for i=1:1:length(X)
    plot([Xc(i),Xc(i)],[0,10],'r--')
    plot([0,10],[Yc(i),Yc(i)],'r--')
end

Xg=dx/2:dx:10-dx/2;
Yg=dx/2:dx:10-dx/2;
k=1;
XY=zeros(length(Xg)^2,5);
for i=1:1:length(Xg)
    for j=1:1:length(Xg)
        XY(k,:)=[Xg(i),Yg(j),i,j,5];
        k=k+1;
    end
end
plot(XY(:,1),XY(:,2),'b+')
%%  Traficability values
% XY(:,3)=5*ones(size(XY,1),1); % default value in the regions where u can go
Pu=[xu,yu;10,10;0,10;xu(1),yu(1)];
Pd=[xd,yd;10,0;xd(1),yd(1)];
Pc=[xcu,ycu;xcd(end:-1:1),ycd(end:-1:1);0,0;xcu(1),ycu(1)];
% figure(2)
% plot(Pu(:,1),Pu(:,2),Pd(:,1),Pd(:,2),Pc(:,1),Pc(:,2))

N=size(XY,1);

[in,on] = inpolygon(XY(:,1),XY(:,2),Pu(:,1),Pu(:,2));
XY(find((in+on)==1),5)=10;

[in,on] = inpolygon(XY(:,1),XY(:,2),Pd(:,1),Pd(:,2));
XY(find((in+on)==1),5)=10;

[in,on] = inpolygon(XY(:,1),XY(:,2),Pc(:,1),Pc(:,2));
XY(find((in+on)==1),5)=1;

M=zeros(length(Xg),length(Xg));
for i=1:1:length(Xg)
    for j=1:1:length(Xg)
        indi=find(XY(:,3)==i);
        indj=find(XY(:,4)==j);
        b = intersect(indi,indj);
        M(i,j)=XY(b,5);
    end
end
% M=M(end-1,:);
% M=M';
%% Running the filtering algorithm
clc
close all
xk=[0.5,2,0.1,0]';
Pk=diag([0.3^2,0.3^2,(0.01)^2,(0.01)^2]);
hx=@(x)x;
gx=@(x,xc)10*norm(x(1:2)-xc(1:2));

Q=diag([0.001^2,0.001^2,(0.001)^2,(0.001)^2]);
Rmanv=diag([1^2,1^2]);

dt=0.2;
for k=2:1:300
    
    %     propagation
    [Xk,wk]=UT_sigmapoints(xk,Pk,2);
    Xk1=zeros(size(Xk));
    %     Xk1(1,:)= KIRB_UM_eg_dyn_disc_manv(Xk(1,:)',dt,XY,dx);
    for i=1:1:length(wk)
        Xk1(i,:)=KIRB_UM_eg_dyn_disc(Xk(i,:)',dt);
        %         Xk1(i,:)=KIRB_CT_eg_dyn_disc(Xk(i,:)',dt);
        %         Xk1(i,:)=KIRB_UM_eg_dyn_disc_manv(Xk(i,:)',dt,XY,dx);
    end
    
    [xk1,Pk1]=MeanCov(Xk1,wk);
    Pk1=Pk1+Q;
    [Xk1,wk1]=UT_sigmapoints(xk1,Pk1,2);
    
    
    %update ... use maneuverability data as measurements
    
    ym=[];
    Vc=get_side_ind(xk1,XY,dx);
    for j=1:1:9 % 9 measurements
        ym=vertcat(ym,Vc(j,4:5)');
    end
    Y=zeros(length(wk1),9*2);
    RR=0;
    for i=1:1:length(wk1)
        V=get_side_ind(Xk1(i,:),XY,dx);
        RS=[];
        for j=1:1:9
            Y(i,2*(j-1)+1:2*j)=Xk1(i,1:2);
            if V(j,6)>=0
                RS=blkdiag(RS,100*(V(j,3)*gx(Xk1(i,1:2),V(j,4:5))^2)*Rmanv);
            else
                RS=blkdiag(RS,1000*(V(j,3)*gx(Xk1(i,1:2),V(j,4:5))^2)*Rmanv);
            end
        end
        RR=RR+wk1(i)*RS;
    end
    
    [mz,Pz]=MeanCov(Y,wk1);
    Pz=Pz+RR;
    Pcc=CrossCov(Xk1,xk1,Y,mz,wk1);
    [xk,Pk]=KalmanUpdate(xk1,Pk1,mz,Pz,Pcc,ym);
    %     keyboard
    
    plot(xu,yu,'k',xd,yd,'k',xcu,ycu,'k',xcd,ycd,'k')
    hold on
    plot_1sig_contour(xk,Pk)
    for i=1:1:length(Xc)
        plot([Xc(i),Xc(i)],[0,10],'r--')
        plot([0,10],[Yc(i),Yc(i)],'r--')
    end
    plot(XY(:,1),XY(:,2),'b+')
    axis([0 10 0 10])
    hold off
%         keyboard
    pause(0.1)
end