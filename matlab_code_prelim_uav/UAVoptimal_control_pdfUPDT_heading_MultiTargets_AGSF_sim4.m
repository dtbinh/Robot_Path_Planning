%Optimal way point trajecrtory of UAV using dynamic programming
% radial penalty and limited view and pdf update 
% heading angle of sensor UAV
% Multi gaussian multi moving target AGSF
clear
clc
close all
%% setup grid and index ir!!

Nn=8; %the number of nodes on each axis

xc=100;% opp. corner to origin
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

model.fn=5;
model.fx=@KIRB_CT_eg_dyn_disc;
T=1;
model.para_dt(1)=T;
L1=0.02;
L2=0.000001;
model.Q=L1*[T^3/3,0,T^2/2,0,0;
        0,T^3/3,0,T^2/2,0;
        T^2/2,0,T,0,0;
        0,T^2/2,0,T,0;
        0,0,0,0,T*L2/L1];
model.sQ=sqrtm(model.Q);



GMM.method='cut6';
GMM.ghpts=3;
GMM.wtup=0;



model.R=[0.5^2,0;0,(1*pi/180)^2];
model.hn=2;
Rmax=50;%field of view
alpha=30*pi/180;% max angle of view
model.hx=@(x,para,xc1)[sqrt((x(2)-xc1(2))^2+(x(1)-xc1(1))^2);atan2(x(2)-xc1(2),x(1)-xc1(1))];
model.gx=@(x,para,xc1)sensor_fov_R(x,para,xc1,Rmax);



%% MI parameters computation
model.marg_nx=[1,2];
model.InfoQuadN=0;
model.InfoMethod='cut6';
model.Infoghpts=4;
model.rmax=Rmax;

GMM.mu=[10,10,-1,5,-5*pi/180;10,10,2,3,-3*pi/180];
P1=diag([1,1,0.5,0.5,0.001]);
P2=diag([1,1,0.5,0.5,0.001]);

GMM.P=[reshape(P1,1,model.fn^2);reshape(P2,1,model.fn^2)];
GMM.w=[0.5,0.5]';
GMM=GMM_moms(GMM);
sigma = cat(3,P1,P2);
 obj = gmdistribution(GMM.mu,sigma,GMM.w');
 
    sens.phi=pi/2;
    sens.alpha=alpha;
    sens.xc1=[100,0];

%% cost function parameters
B=0;
si=@(x)(x-[0,100])*B*(x-[0,100])';
model.si=si;
A=0;
model.A=A;
C=1;
model.C=C;
Nt=25;% 10 steps 
J=zeros(Nt,Nn^2,4); % optima cost to go function for all nodes
mu=zeros(Nt,Nn^2,4); % optimal policy for all nodes

%% Dynamic programming.

[mu,Hnm,PDFk,I]=InfoOptimalTraj_AGSF(X,Nn,xc,yc,0,Nt,GMM,model,sens);

orig_PDFk=PDFk;
%% Generating true trajectory

%x0 = random(obj)';
x0=[11.9444    9.4282   -1.1768    3.8903   -0.1024]';
[t,x_mc]=ode45_discc(@(x)KIRB_CT_eg_dyn_disc(x,model.para_dt),0,T,(Nt)*T,x0,1e-200);

sens.x_mc=x_mc;

%% Simulating the optimal trajectory and filtered pdf of the target

ci=Nn;% initial condition
ch=2; % initial heading
cs=[ci,ch];

 close all

mov(1:Nt-1) = struct('cdata', [],...
                        'colormap', []);
Xtraj=X(ci,:);
Hc=zeros(1,Nt-1);

p=0;

newpath=[1,1,0];
newcol=cell(1,10);
newcol{1}='bo-';
newcol{2}='rs-';
newcol{3}='g^-';
newcol{4}='md-';
newcol{5}='bo-';
newcol{6}='rs-';
newcol{7}='g^-';
newcol{8}='md-';
newcol{9}='bo-';
newcol{10}='rs-';

model.frequpdt=ones(1,Nt);
model.timeupdt=1;

sensor_config=zeros(Nt,4);
%% plot intiial conditions
% figure
% plot_pdf(GMM,xc,yc,model)
% hold 
% plot(x_mc(1,1),x_mc(1,2),'k*-')

%     for hg=1:1:size(newpath,1)
%             plot(Xtraj(1,1),Xtraj(1,2),newcol{newpath(hg,1)},'linewidth',2,'MarkerSize',6)
%     end

%     plot_sens_view(X(ci,:),I(1,ci,3),I(1,ci,2),Rmax,'b')
%     axis([-10,xc+10,-10,yc+10])
%     xlabel('x-axis')
%     ylabel('y-axis')
%      plot_prop_paper  
new_PDFk=cell(1,Nt);
sens.phi=I(1,ci,2);
sens.alpha=I(1,ci,3);
for k=1:1:Nt-1
% calculating the entropy and recompute traj if necessary
new_PDFk{k}=GMM;
sensor_config(k,:)=[ci,ch,sens.phi,sens.alpha];
     Hc(k)=entropy_gmm(GMM);
     DKL=KL_gmm(PDFk{k-p},GMM)
     if DKL>=10
         disp('Computing New Optimal Trajectory')
        
         [mu,HH,PDFk,I]=InfoOptimalTraj_AGSF(X,Nn,xc,yc,k,Nt,GMM,model,sens);

         p=k-1;
         newpath=vertcat(newpath,[newpath(end,1)+1,k,0]);
         
     end
   
    %move the sensor optimally and the pdf
    cs=gridmove_heading(X,xc,yc,Nn,cs,mu(k-p,ci,ch));
    ci=cs(1);
    ch=cs(2);
    model.frequpdt=zeros(1,Nt);
    model.timeupdt=1;
    [GMM,model,sens]=GMM_Filt_MOVETarg(GMM,model,sens,k);
     % at time k+1 max info
    I(k+1-p,ci,:)=max_Info(X(ci,:),GMM,model);
    sens.phi=I(k+1-p,ci,2);
    sens.alpha=I(k+1-p,ci,3);
    sens.xc1=X(ci,:);
    % now update the pdf with measurement
    model.frequpdt=ones(1,Nt);
    model.timeupdt=0;    
    [GMM,model,sens]=GMM_Filt_MOVETarg(GMM,model,sens,k+1);
    Xtraj=vertcat(Xtraj,X(ci,:));
    newpath(end,end)=k+1;
    
    figure
    plot_pdf(GMM,xc,yc,model)
    hold on
    plot(x_mc(1:k+1,1),x_mc(1:k+1,2),'k*-')
k+1


    for hg=1:1:size(newpath,1)
            plot(Xtraj(newpath(hg,2):newpath(hg,3),1),Xtraj(newpath(hg,2):newpath(hg,3),2),newcol{newpath(hg,1)},'linewidth',2,'MarkerSize',6)
     end

    plot_sens_view(X(ci,:),I(k+1-p,ci,3),I(k+1-p,ci,2),Rmax,'b')
    axis([-10,xc+10,-10,yc+10])
    xlabel('x-axis')
    ylabel('y-axis')
     plot_prop_paper    

    pause(0.3)
end
movie(mov)
 movie2avi(mov,'MOVsensSTATtarg_LimView','compression','none')
% figure
% plot(1:1:length(Hnm),-Hnm,1:1:length(Hc),-Hc)
% xlabel('Time step')
% ylabel('Joint Entropy')
% legend('No meas','With meas')
%    plot_prop_paper
save('UAV_movTarg3','orig_PDFk','new_PDFk','Hnm','Hc','I','X','Xtraj','newpath','newcol','model','sens','x_mc','xc','yc','Nn','Nt','sensor_config')