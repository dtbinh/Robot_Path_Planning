%% solving the optimal constrol with information potential (DP formation)
%% load the required templates
clear
load('0_0_1p5708s40s-1s1s8s18s5s0s0')
clc
close all
Ntar=3;
prior_GMM.mu=cell(Ntar,1);
prior_GMM.P=cell(Ntar,1);
prior_GMM.w=cell(Ntar,1);

prior_GMM.mu{1}=[300,300];
prior_GMM.P{1}=[100,50,50,100];
prior_GMM.w{1}=[1];


prior_GMM.mu{2}=[70,200,0,4,-0.03;75,210,0,4,-0.03];
prior_GMM.P{2}=[reshape(diag([1,1,0.1,0.1,(0.0001)^2]),1,25);reshape(diag([1,1,0.1,0.1,(0.0001)^2]),1,25)];
prior_GMM.w{2}=[0.5;0.5];

prior_GMM.mu{3}=[200,80];
prior_GMM.P{3}=[100,0,0,10];
prior_GMM.w{3}=[1];

prior_GMM.mu{4}=[200,200];
prior_GMM.P{4}=[100,0,0,10];
prior_GMM.w{4}=[1];

prior_GMM.mu{5}=[270,30];
prior_GMM.P{5}=[10,5,5,200];
prior_GMM.w{5}=[1];

% load('shit_GMM')

% measurment model
h=cell(Ntar,1);
g=cell(Ntar,1);
% h{1}=@(x,s)[(x(1)-s(1));(x(2)-s(2))];
% h{2}=@(x,s)[(x(1)-s(1));(x(2)-s(2))];

% h{1}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];
% h{2}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];
% h{3}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];
% h{4}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];
% h{5}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];

h{1}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));
h{2}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));
h{3}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));
h{4}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));
h{5}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));

% h{1}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));
% h{2}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));

r=60; % is the field of view 

% variance penalty model
R=cell(Ntar,1);
R{1}=diag([0.1745].^2);%+(r/4)^2;
R{2}=diag([0.1745].^2);%+(r/4)^2;
R{3}=diag([0.1745].^2);
R{4}=diag([0.1745].^2);
R{5}=diag([0.1745].^2);

g{1}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);
g{2}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);
g{3}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);
g{4}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);
g{5}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);

G{1}=@(x,s)1;
G{2}=@(x,s)1;
G{3}=@(x,s)1;
G{4}=@(x,s)1;
G{5}=@(x,s)1;

nz=cell(2,1);
nz{1}=1;
nz{2}=1;
nz{3}=1;
nz{4}=1;
nz{5}=1;

Q{1}=5*[0,0;0,0];
Q{2}=diag([0.1^2,0.1^2,0.02^2,0.02^2,(0.0001)^2]);
Q{3}=5*[0,0;0,0];
Q{4}=5*[0,0;0,0];
Q{5}=5*[0,0;0,0];

para{1}=0;
para{2}=0;
para{3}=0;
para{4}=0;
para{5}=0;



Nu=9;
Ns=2;
d=40;
xlim=400;
ylim=400;
t0=0;
dt=5;
tf=100;

targ_dyn{1}=@(x,para)eye(2)*x;
targ_dyn{2}=@(x,para)CT_turn_model(x,para,dt);
targ_dyn{3}=@(x,para)eye(2)*x;
targ_dyn{4}=@(x,para)eye(2)*x;
targ_dyn{5}=@(x,para)eye(2)*x;


s0{1,1}=[1,1,2];
s0{2,1}=[10,2,2];

sts=[1 2];
targ_pos=[0,0];
for i=1:1:Ntar
    PP=reshape(prior_GMM.P{i}(1,:),length(prior_GMM.mu{i}(1,:)),length(prior_GMM.mu{i}(1,:)));
targ_pos(i,:)=mvnrnd(prior_GMM.mu{i}(1:2)',PP(1:2,1:2));
end
[Ygrid,Xgrid]=meshgrid(0:d:xlim);

prrGmm.mu=prior_GMM.mu{2};
prrGmm.P=prior_GMM.P{2};
prrGmm.w=prior_GMM.w{2};
obj=GMM2obj(prrGmm);
 YY= random(obj);
 targ_pos(2,:)=YY(1:2);
 mov2=YY;
% [M,J]=dubin_DP(t0,dt,tf,[1,1,2],Nu,Ns,Ntar,d,xlim,ylim,@UAV_grid_motion,targ_dyn,prior_GMM,nz,h,g,R,Q,para,templates1,sts);

%% moving the UAV around th grid
for i=1:1:Ns
c0{i,1}=s0{i,1};
XY{i,1}=[Xgrid(s0{i,1}(1,1),s0{i,1}(1,2)),Ygrid(s0{i,1}(1,1),s0{i,1}(1,2))];
end
Hist.u{1}=0;
Hist.priorGMM{1}=prior_GMM;
cc=0;
tt=2;
thresh=2; %threshold to ignore the info from this
coll{1,1}='r';
coll{2,1}='b';
while(cc<=6)
    cc=cc+1;
disp('computing new trajectory')    
% tic
[M,J]=dubin_DP(t0,dt,tf,c0,Nu,Ns,Ntar,d,xlim,ylim,@UAV_grid_motion,targ_dyn,prior_GMM,nz,h,g,R,Q,para,templates1,sts,thresh,'moving');
% toc
U=zeros(1,Ns);
for k=1:1:5%(length(t0:dt:tf)-1)
    k
    figure(1)
            plot_UAV_GMM(c0,d,xlim,prior_GMM,targ_pos,sts,r)
            hold on;
            for sdsd=1:1:Ns
            plot(XY{sdsd}(:,1),XY{sdsd}(:,2),coll{sdsd},'linewidth',2)
            end
            plot(mov2(:,1),mov2(:,2),'ko-')
            box on
            hold off;
%             plot_prop_paper
            xlabel('x')
            ylabel('y')
%             keyboard
            pause(0.5)
            axis square
            saveas(gcf,strcat('UAV2_move_targs_bear',num2str(tt-1)), 'pdf')
            saveas(gcf,strcat('UAV2_move_targs_bear',num2str(tt-1)), 'fig')
            % move the sensors
    for ss=1:1:Ns
        u=M{ss,k}(c0{ss}(1,1),c0{ss}(1,2),c0{ss}(1,3));
        J{ss,k}(c0{ss}(1,1),c0{ss}(1,2),c0{ss}(1,3));
        U(ss)=u;
        [Xf,V]=UAV_grid_motion(c0{ss}(1,1),c0{ss}(1,2),c0{ss}(1,3),u,xlim,d,templates1);
        XX=templates2{u,1}(:,[1,2]);
        Y=XX;
        if c0{ss}(1,3)==1
            Y(:,1)=XX(:,2);
            Y(:,2)=-XX(:,1);
        end
        if c0{ss}(1,3)==2
            Y=XX;
        end
        if c0{ss}(1,3)==3
              Y(:,1)=-XX(:,2);
              Y(:,2)=XX(:,1);
        end
        if c0{ss}(1,3)==4
          Y(:,1)=-XX(:,1);
          Y(:,2)=-XX(:,2);
        end
        XY{ss}=vertcat(XY{ss},Y+repmat([Xgrid(c0{ss}(1,1),c0{ss}(1,2)),Ygrid(c0{ss}(1,1),c0{ss}(1,2))],size(templates2{u,1},1),1));
%           keyboard
        c0{ss}(1,1)=Xf(1);c0{ss}(1,2)=Xf(2);c0{ss}(1,3)=Xf(3); 
    end
   
% make a measurement and update
    % update prior
    for nt=1:1:Ntar
    prGMM1.mu=prior_GMM.mu{nt};prGMM1.P=prior_GMM.P{nt};prGMM1.w=prior_GMM.w{nt};   
    prGMM1=prior_prop_GMM(prGMM1,targ_dyn{nt},para{nt},Q{nt});
    prior_GMM.mu{nt}=prGMM1.mu;prior_GMM.P{nt}=prGMM1.P;prior_GMM.w{nt}=prGMM1.w;
          % propagate targets true positions
          if nt==2
              mov2=vertcat(mov2,targ_dyn{nt}(mov2(end,:)',para)');
              targ_pos(2,:)=mov2(end,1:2);
          else
          targ_pos(nt,:)=targ_dyn{nt}(targ_pos(nt,:)',para);
          end
     end
    % update with meas
    for nt=1:1:Ntar
        z=[];
        for ss=1:1:Ns
           pp=[Xgrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2))]; 
        if norm(targ_pos(nt,:)-pp)<r
            z=h{nt}(targ_pos(nt,:),pp)+0.4*sqrt(R{nt})*randn(nz{nt},1); 
            z
            prGMM1.mu=prior_GMM.mu{nt};prGMM1.P=prior_GMM.P{nt};prGMM1.w=prior_GMM.w{nt};   
            measGMM=post_GMM(prGMM1,z,h{nt},G{nt},0.4*R{nt},pp);
            prior_GMM.mu{nt}=measGMM.mu;prior_GMM.P{nt}=measGMM.P;prior_GMM.w{nt}=measGMM.w;
        end  
        end
        prGMM1.mu=prior_GMM.mu{nt};prGMM1.P=prior_GMM.P{nt};prGMM1.w=prior_GMM.w{nt};
        
    end
    Hist.u{tt}=U;
    Hist.priorGMM{tt}=prior_GMM;
    tt=tt+1;
end % implement first 10 steps
% keyboard
end %cont loop while
save('control_UAV','Hist','targ_dyn','g','h','G','R','Q','d','t0','dt','tf','para','nz','xlim','Ntar','Ns','Nu','para','r','s0')

%% plotiing the trajectory along with
% for 