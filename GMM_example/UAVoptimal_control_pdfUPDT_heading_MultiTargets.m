%Optimal way point trajecrtory of UAV using dynamic programming
% radial penalty and limited view and pdf update 
% heading angle of sensor UAV
% Now there are 2 targets !!!!!!!!!
% clear
clc

%% setup grid and index ir!!

Nn=10; %the number of nodes on each axis

xc=100;% opp. x corner to origin
yc=100;% opp. y corner to origin

dx=xc/(Nn-1);
dy=yc/(Nn-1);
X=zeros(Nn^2,2); % first 2 cols for position, and last 4 cols for admiss. control number from 1,2,3,4,,,control u=1-right, 2-up, 3-left, 4-down 0-not there
U=repmat([1,2,3,4],Nn^2,1);
x=0;
y=0;
k=1;
for j=1:1:Nn %y incr
    for i=1:1:Nn %x incr
  x=(i-1)*dx;
  y=(j-1)*dy;
        if i==1
            x=0;
           % X(k,:)=[x,y];
            U(k,3)=0; %no left control
        end
        if i==Nn
            %X(k,:)=[x,y];
            U(k,1)=0;
        end
        if j==1
           % X(k,:)=[x,y];
            U(k,4)=0;
        end
        if j==Nn
           % X(k,:)=[x,y];
            U(k,2)=0;
        end
           X(k,:)=[x,y];
        k=k+1;
    end
end
%% sensor model and target dynamics

model.fn=4;
model.fx=@(x,paradt)x;
model.para_dt(1)=0;

GMM.method='gh';
GMM.ghpts=4;
GMM.wtup=0;
model.frequpdt=1;


model.R=0.008*eye(2);
model.hn=2;
Rmax=27;%field of view
alpha=20*pi/180;% max angle of view
model.hx=@(x,para,xc1)[atan2(x(2)-xc1(2),x(1)-xc1(1));atan2(x(4)-xc1(2),x(3)-xc1(1))];
model.gx=@(x,para,xc1)[sensor_limview_penalty(x(1:2),xc1,Rmax,para),0;0,sensor_limview_penalty(x(3:4),xc1,Rmax,para)];
%% MI parameters computation
model.marg_nx=0;
model.InfoQuadN=0;
model.InfoMethod='cut4';
model.Infoghpts=3;

GMM.mu=[20,45,65,80;
        20,40,70,80];
P1=[5,0,0,0;0,50,0,0;0,0,30,0;0,0,0,5];
P2=[50,0,0,0;0,5,0,0;0,0,5,0;0,0,0,30];

GMM.P=[reshape(P1,1,model.fn^2);reshape(P2,1,model.fn^2)];
GMM.w=[0.5;0.5];
GMM=GMM_moms(GMM);
sigma = cat(3,P1,P2);
 obj = gmdistribution(GMM.mu,sigma,GMM.w');
 obj1=gmdistribution(GMM.mu(1:2,1:2),cat(3,P1(1:2,1:2),P2(1:2,1:2)),GMM.w');
 obj2=gmdistribution(GMM.mu(1:2,3:4),cat(3,P1(3:4,3:4),P2(3:4,3:4)),GMM.w');
 % true target position
%  x_truth=[33,80];
 x_truth=[ 19.6734   41.2381   74.2133   78.0418];%random(obj);%mvnrnd(GMM.mean',GMM.cov, 1)
%% get the information value and optimal config at every node
[xx1,yy1]=meshgrid(GMM.mean(1,1)-30:1:GMM.mean(1,1)+30,GMM.mean(2,1)-30:1:GMM.mean(2,1)+30);
[xx2,yy2]=meshgrid(GMM.mean(3,1)-30:1:GMM.mean(3,1)+30,GMM.mean(4,1)-30:1:GMM.mean(4,1)+30);

CC1=zeros(size(xx1));
CC2=zeros(size(xx2));

for i=1:1:size(xx1,1)
    for j=1:1:size(xx1,2)
        CC1(i,j)=pdf(obj1,[xx1(i,j) yy1(i,j)]);
         CC2(i,j)=pdf(obj2,[xx2(i,j) yy2(i,j)]);
%  CC(i,j)=pdf(obj,[xx(i,j) yy(i,j)]);
    end
end
contour(xx1,yy1,CC1,10)
hold on
contour(xx2,yy2,CC2,10)
 axis([-10,xc+10,-10,yc+10])
 hold off
 
I=zeros(size(X,1),3); %first col is information and next to cols are [phi,alpha]
parfor i=1:1:size(X,1)
    i
     tic
    I(i,:)=max_Info(X(i,:),GMM,model);
     toc
%     contour(xx1,yy1,CC1,10)
%     hold on
%     contour(xx2,yy2,CC2,10)
%     plot(X(1:i,1),X(1:i,2),'bo')
%     plot_sens_view(X(i,:),I(i,3),I(i,2),Rmax,'b')
%     axis([-10,xc+10,-10,yc+10])
% pause(0.1)
end
 save('InfoFilterMULTItarget','I')
load('InfoFilterMULTItarget')
%% cost function parameters
B=0;
si=@(x)(x-[0,100])*B*(x-[0,100])';
model.si=si;
A=0;
model.A=A;
C=2;
model.C=C;
Nt=27;% 10 steps 
J=zeros(Nt,Nn^2,4); % optima cost to go function for all nodes
mu=zeros(Nt,Nn^2,4); % optimal policy for all nodes

%% Dynamic programming.
 for i=1:1:Nn^2
     for j=1:1:4
    J(Nt,i,j)=si(X(i,:));
     end
 end
 
for k=Nt-1:-1:1
   
    for ci=1:1:Nn^2
        for ch=1:1:4
        JJ=0;
        u=0;
        p=1;
        for j=-1:1:1
            ns=gridmove_heading(X,xc,yc,Nn,[ci,ch],j);
            if isnan(ns(1))~=1
                u(p)=j;
                JJ(p)=-1*I(ci,1)+C*u(p)^2+J(k+1,ns(1),ns(2));%+(X(ci,:)-GMM.mean')*A*(X(ci,:)-GMM.mean')'
                p=p+1;
            else
                JJ(p)=1e10;
                u(p)=j;
                p=p+1;
            end
        end
                [J(k,ci,ch),ind]=min(JJ);
                mu(k,ci,ch)=u(min(ind));
        end
    end
end


%% Simulating the optimal trajectory and filtered pdf of the target

ci=6;% initial condition
ch=2; % initial heading
cs=[ci,ch];
 close all

mov(1:Nt-1) = struct('cdata', [],...
                        'colormap', []);
Xtraj=X(ci,:);

Hc=zeros(1,Nt-1);
% Hp=0.5*log((2*pi*exp(1))^model.fn*det(reshape(GMM.P,model.fn,model.fn)));
mu1=GMM.mean;
P1=GMM.cov;
p=0;
newpath=[1,1,0];
newcol=cell(1,4);
newcol{1}='bo-';
newcol{2}='rs-';
newcol{3}='g^-';
newcol{4}='md-';
for k=1:1:Nt-1
    close
% figure
    GMM.P=real(GMM.P);
    GMM.cov=real(GMM.cov);
    % calculating the entropy
     Hc(k)=0.5*log((2*pi*exp(1))^model.fn*det(GMM.cov));
     mu2=GMM.mean;
%      ii=find(GMM.cov<0.5 & GMM.cov>0);
     P2=GMM.cov;
%      P2(ii)=0.5;
     DKL=real(log(sqrt(det(P2)/det(P1)))-0.5*model.fn+0.5*(mu1-mu2)'*inv(P2)*(mu1-mu2)+0.5*trace(inv(P2)*P1))
     if DKL>=20
         disp('Computing New Optimal Trajectory')
        mu1=mu2;
        P1=P2;
         mu=InfoOptimalTraj(X,Nn,xc,yc,k,Nt,GMM,model);  
         p=k-1;
         newpath=vertcat(newpath,[newpath(end,1)+1,k-1,0]);
         
     end
     newpath(end,end)=k;
    % plot the current target pdf
%     if k>8 && rem(k,2)==0
%      figure
%    for i=1:1:size(xx1,1)
%     for j=1:1:size(xx1,2)
%         PP=reshape(GMM.P,model.fn,model.fn);
%         if min(eig(PP(1:2,1:2)))<0.1
%             PP(1:2,1:2)=0.1*eye(2);
%         end
%         if min(eig(PP(3:4,3:4)))<0.1
%             PP(3:4,3:4)=0.1*eye(2);
%         end 
% %         CC1(i,j)=mvnpdf([xx1(i,j),yy1(i,j)],GMM.mu(1:2),PP(1:2,1:2));
% %         CC2(i,j)=mvnpdf([xx2(i,j),yy2(i,j)],GMM.mu(3:4),PP(3:4,3:4));
%   CC(i,j)=pdf(obj,[xx(i,j) yy(i,j)]);
%     end
%    end

PP1=reshape(GMM.P(1,:),model.fn,model.fn);
PP2=reshape(GMM.P(2,:),model.fn,model.fn);

if min(eig(PP1(1:2,1:2)))<0.5 || min(eig(PP2(1:2,1:2)))<0.5
    PP1(1:2,1:2)=0.5*eye(2);
    PP2(1:2,1:2)=0.5*eye(2);
end
if min(eig(PP1(3:4,3:4)))<0.5 || min(eig(PP2(3:4,3:4)))<0.5
    PP1(3:4,3:4)=0.5*eye(2);
    PP2(3:4,3:4)=0.5*eye(2);
end
obj1=gmdistribution(GMM.mu(1:2,1:2),cat(3,PP1(1:2,1:2),PP2(1:2,1:2)),GMM.w');
obj2=gmdistribution(GMM.mu(1:2,3:4),cat(3,PP1(3:4,3:4),PP2(3:4,3:4)),GMM.w');
for i=1:1:size(xx1,1)
    for j=1:1:size(xx1,2)
        CC1(i,j)=pdf(obj1,[xx1(i,j) yy1(i,j)]);
         CC2(i,j)=pdf(obj2,[xx2(i,j) yy2(i,j)]);
%  CC(i,j)=pdf(obj,[xx(i,j) yy(i,j)]);
    end
end
contour(xx1,yy1,CC1,10)
hold on
contour(xx2,yy2,CC2,10)
  %  plot the sensor
%     if k>=8 && k<=14 && rem(k,2)==0
%     hold on

    plot(x_truth(1),x_truth(2),'o',x_truth(3),x_truth(4),'o','MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',10)
    for hg=1:1:size(newpath,1)
            plot(Xtraj(newpath(hg,2):newpath(hg,3),1),Xtraj(newpath(hg,2):newpath(hg,3),2),newcol{newpath(hg,1)},'linewidth',2,'MarkerSize',6)
    end
    plot_sens_view(X(ci,:),I(ci,3),I(ci,2),Rmax,'b')
    axis([-10,xc+10,-10,yc+10])
    xlabel('x-axis')
    ylabel('y-axis')
  %   legend('Target state pdf','True Target positions')
     plot_prop_paper
%      end
    %     figure(2)
    
    %move the sensor optimally
    cs=gridmove_heading(X,xc,yc,Nn,cs,mu(k-p,ci,ch))
    ci=cs(1);
    ch=cs(2);
%     ci=gridmove(Nn,ci,mu(k,ci));
    Xtraj=vertcat(Xtraj,X(ci,:));
    
    % take measurement and update the target pdf
    I(ci,:)=max_Info(X(ci,:),GMM,model);
    
    sens.phi=I(ci,2);
    sens.alpha=I(ci,3);
    sens.xc1=X(ci,:);
    [GMM,model,sens,ym]=GMM_Filt_statTarg(GMM,model,sens,1,x_truth);
%     
%     ii=find(GMM.P>0 & GMM.P<0.1);
%     GMM.P(ii)=0.1;
    
    pause(0.5)
    saveas(gcf,strcat('twoStat_Targets_',num2str(k),'.jpg'))
    mov(k) = getframe;
    pause(0.3)
    
end
figure
plot(Hc)
xlabel('Time step')
ylabel('Joint Entropy')
   plot_prop_paper
movie(mov)
 movie2avi(mov,'MOVsensSTATtarg_LimView123','compression','none')
%% just for test section
% for k=1:1:size(X,1)
%     plot(X(1:k,1),X(1:k,2),'bo-')
%     axis([-10,xc+10,-10,yc+10])
%     pause(0.1)
% end
