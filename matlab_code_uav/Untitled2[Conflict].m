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

% save('EgNaval','xu','yu','xd','yd','xcu','ycu','xcd','ycd')
%%
load('EgNaval')
plot(xu,yu,'k',xd,yd,'k',xcu,ycu,'k',xcd,ycd,'k')
hold on
dx=0.5;
X=0:dx:10;
Y=0:dx:10;


for i=1:1:length(X)
    plot([X(i),X(i)],[0,10],'r--')
    plot([0,10],[Y(i),Y(i)],'r--')
end

Xg=dx/2:dx:10-dx/2;
Yg=dx/2:dx:10-dx/2;
k=1;
XY=zeros(length(Xg)^2,5);
for i=1:1:length(Xg)
    for j=1:1:length(Xg)
        XY(k,:)=[Xg(i),Yg(j),i,j,0.5];
        k=k+1;
    end
end
plot(XY(:,1),XY(:,2),'b+')
%%  Traficability values
% XY(:,3)=5*ones(size(XY,1),1); % default value in the regions where u can go
Pu=[xu,yu;10,10;0,10;xu(1),yu(1)];
Pd=[xd,yd;10,0;xd(1),yd(1)];
Pc=[xcu,ycu;xcd(end:-1:1),ycd(end:-1:1);0,0;xcu(1),ycu(1)];
figure(2)
plot(Pu(:,1),Pu(:,2),Pd(:,1),Pd(:,2),Pc(:,1),Pc(:,2))

N=size(XY,1);

[in,on] = inpolygon(XY(:,1),XY(:,2),Pu(:,1),Pu(:,2));
XY(find((in+on)==1),5)=0;

[in,on] = inpolygon(XY(:,1),XY(:,2),Pd(:,1),Pd(:,2));
XY(find((in+on)==1),5)=0;

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

% %Traficability surface
% sf = fit(XY(:,1:2),XY(:,5),'linearinterp')
% [xx,yy]=meshgrid(0:0.1:10);
% 
% Mp=zeros(size(xx));
% for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%         Mp(i,j)=sf([xx(i,j),yy(i,j)]);
%     end
% end
% mesh(xx,yy,Mp)
% hold on
% plot3(XY(:,1),XY(:,2),XY(:,5),'bo','MarkerSize',6)
%% Running the filtering algorithm
clc
close all
xk=[1,1,0.4,0]';
Pk=diag([0.1^2,0.1^2,(0.05)^2,(0.05)^2]);

Q=diag([0.1^2,0.1^2,(0.01)^2,(0.01)^2]);
Rmanv=(1*pi/180)^2;
figure(1)
clf
dt=0.5;
for k=2:1:100
    
    %     propagation
    [Xk,wk]=UT_sigmapoints(xk,Pk,2);
    Xk1=zeros(size(Xk));
    for i=1:1:length(wk)
        Xk1(i,:)=KIRB_UM_eg_dyn_disc_manv(Xk(i,:)',dt,XY,dx);
    end
    [xk1,Pk1]=MeanCov(Xk1,wk);
    Pk1=Pk1+Q;
    [Xk1,wk1]=UT_sigmapoints(xk1,Pk1,2);
    
    disp('OK prop')
%update ... use maneuverability data as measurements
%     Y=zeros(length(wk),1);
%     RR=0;
%     WK=wk;
%     for i=1:1:length(wk)
%         [a,b]=min(sqrt(sum((repmat(Xk1(i,1:2),N,1)-XY(:,1:2)).^2,2)));
%         thdt=atan2(XY(b,2)-Xk1(i,2),XY(b,1)-Xk1(i,1))/dt;
%         Y(i,1)=thdt;
%         manv=XY(b,3);
%         RR=RR+wk1(i)*Rmanv;
%         WK(i)=wk1(i)*(10-manv);
%     end
%     WK=WK/sum(WK);
%     [mz,Pz]=MeanCov(Y,WK);
%     Pz=Pz+RR;
%     Pcc=CrossCov(Xk1,xk1,Y,mz,wk1);
%     [xk,Pk]=KalmanUpdate(xk1,Pk1,mz,Pz,Pcc,0);
%     keyboard

    xk=xk1;
    Pk=Pk1;
    
    
    plot(xu,yu,'k',xd,yd,'k',xcu,ycu,'k',xcd,ycd,'k')
    hold on
    plot(Xk1(:,1),Xk1(:,2),'go')
    plot_1sig_contour(xk,Pk)
    for i=1:1:length(Xg)
        plot([Xg(i),Xg(i)],[0,10],'r--')
        plot([0,10],[Yg(i),Yg(i)],'r--')
    end
    plot(XY(:,1),XY(:,2),'b+')
    axis([0 10 0 10])
%     keyboard
    pause(0.1)
    clf
end