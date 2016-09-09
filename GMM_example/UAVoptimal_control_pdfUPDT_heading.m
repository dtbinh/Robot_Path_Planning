%Optimal way point trajecrtory of UAV using dynamic programming
% radial penalty and limited view and pdf update 
% heading angle of sensor UAV
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

model.fn=2;
model.fx=@(x,paradt)x;
model.para_dt(1)=0;

GMM.method='gh';
GMM.ghpts=5;
GMM.wtup=0;
model.frequpdt=1;


model.R=1;
model.hn=1;
Rmax=30;%field of view
alpha=20*pi/180;% max angle of view
model.hx=@(x,para,xc1)atan2(x(2)-xc1(2),x(1)-xc1(1));
model.gx=@(x,para,xc1)sensor_limview_penalty(x,xc1,alpha,Rmax,para);
model.marg_nx=0;

GMM.mu=[35,75];
P1=[5,0;0,50];

GMM.P=reshape(P1,1,model.fn^2);
GMM.w=[1];
GMM=GMM_moms(GMM);
sigma = cat(3,P1);
 obj = gmdistribution(GMM.mu,sigma,GMM.w');
 
 % true target position
%  x_truth=[33,80];
 x_truth=[ 36.9470   72.3006]%=mvnrnd(GMM.mean',GMM.cov, 1)
%% get the information value and optimal config at every node
[xx,yy]=meshgrid(GMM.mu(1,1)-30:0.5:GMM.mu(1,1)+30,GMM.mu(1,2)-30:0.5:GMM.mu(1,2)+30);

CC=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        CC(i,j)=mvnpdf([xx(i,j),yy(i,j)],GMM.mean',GMM.cov);
%  CC(i,j)=pdf(obj,[xx(i,j) yy(i,j)]);
    end
end
contour(xx,yy,CC,15)
 axis([-10,xc+10,-10,yc+10])
 
 
% I=zeros(size(X,1),3); %first col is information and next to cols are [phi,alpha]
% parfor i=1:1:size(X,1)
%     i
%      tic
%     I(i,:)=max_Info(X(i,:),GMM,model);
%      toc
% %     contour(xx,yy,CC,15)
% %     hold on
% %     plot(X(1:i,1),X(1:i,2),'bo')
% %     plot_sens_view(X(i,:),I(i,3),I(i,2),Rmax,'b')
% %     axis([-10,xc+10,-10,yc+10])
% % pause(0.1)
% end
% save('InfoFilter','I')
load('InfoFilter')
%% cost function parameters
B=0;
si=@(x)(x-[0,100])*B*(x-[0,100])';
A=0;
C=150;
Nt=25;% 10 steps 
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
                JJ(p)=-100*I(ci,1)+(X(ci,:)-GMM.mean')*A*(X(ci,:)-GMM.mean')'+C*u(p)^2+J(k+1,ns(1),ns(2));
                p=p+1;
            end
        end
                [J(k,ci,ch),ind]=min(JJ);
                mu(k,ci,ch)=u(min(ind));
        end
    end
end


%% Simulating the optimal trajectory and filtered pdf of the target

ci=Nn;% initial condition
ch=2; % initial heading
cs=[ci,ch];
% close all

mov(1:Nt-1) = struct('cdata', [],...
                        'colormap', []);
Xtraj=X(ci,:);
H=zeros(1,Nt-1);
for k=1:1:Nt-1

    % calculating the entropy
    H(k)=0.5*log((2*pi*exp(1))^model.fn*det(GMM.cov));
    
    % plot the current target pdf
    if k>8 && rem(k,2)==0
    figure(k)
    CC=zeros(size(xx));
    for i=1:1:size(xx,1)
        for j=1:1:size(xx,2)
        CC(i,j)=mvnpdf([xx(i,j),yy(i,j)],GMM.mean',GMM.cov);
        end
    end
        contour(xx,yy,CC,15)
  %  plot the sensor
%     if k>=8 && k<=14 && rem(k,2)==0
    hold on
    plot(x_truth(1),x_truth(2),'o','MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',10)
    plot(Xtraj(:,1),Xtraj(:,2),'bo-','linewidth',2)
    plot_sens_view(X(ci,:),I(ci,3),I(ci,2),Rmax,'b')
    axis([-10,xc+10,-10,yc+10])
%     legend('Target state pdf','True Target position','UV Trajectory')
    plot_prop_paper
     end
    %     figure(2)
    
    %move the sensor optimally
    cs=gridmove_heading(X,xc,yc,Nn,cs,mu(k,ci,ch));
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
    
%     pause(0.5)
%     mov(k) = getframe;
    pause(0.3)
%     close
end
figure
plot(H)
   plot_prop_paper
% movie(mov)
%  movie2avi(mov,'MOVsensSTATtarg_LimView','compression','none')
%% just for test section
% for k=1:1:size(X,1)
%     plot(X(1:k,1),X(1:k,2),'bo-')
%     axis([-10,xc+10,-10,yc+10])
%     pause(0.1)
% end
