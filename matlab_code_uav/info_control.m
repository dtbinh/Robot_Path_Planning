%% solving the optimal constrol with information potential (DP formation)
%% load the required templates
clear
load('0_0_1p5708s40s-1s1s8s18s5s0s0')
clc
close all
Ntar=6;
prior_GMM.mu=cell(Ntar,1);
prior_GMM.P=cell(Ntar,1);
prior_GMM.w=cell(Ntar,1);
prior_GMM.truth=cell(Ntar,1);
prior_GMM.MovStat=cell(Ntar,1);
prior_GMM.fn=cell(Ntar,1);

prior_GMM.mu{1}=[160,360];
prior_GMM.P{1}=[100,50,50,100];
prior_GMM.w{1}=[1];
prior_GMM.truth{1}=0;
prior_GMM.MovStat{1}='stationary';
prior_GMM.fn{1}=2;


prior_GMM.mu{2}=[50,120];
prior_GMM.P{2}=[100,0,0,100];
prior_GMM.w{2}=[1];
prior_GMM.truth{2}=0;
prior_GMM.MovStat{2}='stationary';
prior_GMM.fn{2}=2;

prior_GMM.mu{3}=[70,250];
prior_GMM.P{3}=[100,0,0,10];
prior_GMM.w{3}=[1];
prior_GMM.truth{3}=0;
prior_GMM.MovStat{3}='stationary';
prior_GMM.fn{3}=2;

prior_GMM.mu{4}=[320,200];
prior_GMM.P{4}=[100,0,0,10];
prior_GMM.w{4}=[1];
prior_GMM.truth{4}=0;
prior_GMM.MovStat{4}='stationary';
prior_GMM.fn{4}=2;

prior_GMM.mu{5}=[240,60];
prior_GMM.P{5}=[10,5,5,200];
prior_GMM.w{5}=[1];
prior_GMM.truth{5}=0;
prior_GMM.MovStat{5}='stationary';
prior_GMM.fn{5}=2;

prior_GMM.mu{6}=[90,100,1,1;110,100,1,1];
prior_GMM.P{6}=[reshape([10,5,0,0;5,100,0,0;0,0,0.1,0.05;0,0,0.05,0.1],1,16);reshape([100,5,0,0;5,20,0,0;0,0,0.1,0.05;0,0,0.05,0.1],1,16)];
prior_GMM.w{6}=[0.5;0.5];
prior_GMM.truth{6}=0;
prior_GMM.MovStat{6}='moving';
prior_GMM.fn{6}=4;

% load('shit_GMM')

Nu=9;
Ns=2;
d=40;
xlim=400;
ylim=400;
t0=0;
dt=5;
tf=50;

%% measurment model
h=cell(Ns,1);
g=cell(Ns,1);
nz=cell(Ns,1);


% h{1}=@(x,s)[(x(1)-s(1));(x(2)-s(2))];
% h{2}=@(x,s)[(x(1)-s(1));(x(2)-s(2))];

h{1}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];
h{2}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];
h{3}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];
h{4}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];
h{5}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];
h{6}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];


% h{1}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));
% h{2}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));
% h{3}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));
% h{4}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));
% h{5}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));
% h{6}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));

% h{1}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));
% h{2}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));

r=60; % is the field of view 

% variance penalty model
R=cell(Ns,1);
% R{1}=diag([5*pi/180].^2);%+(r/4)^2;
% R{2}=diag([5*pi/180].^2);%+(r/4)^2;
% R{3}=diag([5*pi/180].^2);
% R{4}=diag([5*pi/180].^2);
% R{5}=diag([5*pi/180].^2);
% R{6}=diag([5*pi/180].^2);

R{1}=diag([4].^2);%+(r/4)^2;
R{2}=diag([4].^2);%+(r/4)^2;
R{3}=diag([4].^2);
R{4}=diag([4].^2);
R{5}=diag([4].^2);
R{6}=diag([4].^2);

g{1}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);
g{2}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);
g{3}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);
g{4}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);
g{5}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);
g{6}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);


% g for making true measurments
G{1}=@(x,s)1;
G{2}=@(x,s)1;
G{3}=@(x,s)1;
G{4}=@(x,s)1;
G{5}=@(x,s)1;
G{6}=@(x,s)1;

nz{1}=1;
nz{2}=1;
nz{3}=1;
nz{4}=1;
nz{5}=1;
nz{6}=1;

Q{1}=0.001*[1,0;0,1];
Q{2}=0.001*[1,0;0,1];
Q{3}=0.001*[1,0;0,1];
Q{4}=0.001*[1,0;0,1];
Q{5}=0.001*[1,0;0,1];
Q{6}=blkdiag(2,2,0.1,0.1);


para{1}=0;
para{2}=0;
para{3}=0;
para{4}=0;
para{5}=0;
para{6}=dt;

targ_dyn{1}=@(x,para)eye(2)*x;
targ_dyn{2}=@(x,para)eye(2)*x;
targ_dyn{3}=@(x,para)eye(2)*x;
targ_dyn{4}=@(x,para)eye(2)*x;
targ_dyn{5}=@(x,para)eye(2)*x;
targ_dyn{6}=@KIRB_UM_eg_dyn_disc;


s0{1,1}=[1,1,2];
s0{2,1}=[10,2,2];

sts=[1 2];

for i=1:1:Ntar
%     x0=mvnrnd(prior_GMM.mu{i}',reshape(prior_GMM.P{i},prior_GMM.fn{i},prior_GMM.fn{i}));
%     [t,x1]=ode45_discc(targ_dyn{i},t0,dt,tf,x0,zeros(prior_GMM.fn{i}));
GMM.mu=prior_GMM.mu{i};
GMM.P=prior_GMM.P{i};
GMM.w=prior_GMM.w{i};
obj=GMM2obj(GMM);
x0 = random(obj);
prior_GMM.truth{i}=x0;
end


[Ygrid,Xgrid]=meshgrid(0:d:xlim);

% plot_UAV_GMM2(d,xlim,prior_GMM,sts,60,Ntar)
keyboard
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
while(cc<=15)
    cc=cc+1;
    disp('computing new trajectory')
    % tic
    [M,J]=dubin_DP(t0,dt,tf,c0,Nu,Ns,Ntar,d,xlim,ylim,@UAV_grid_motion,targ_dyn,prior_GMM,nz,h,g,R,Q,para,templates1,sts,thresh);
    % toc
    U=zeros(1,Ns);
    for k=1:1:2%(length(t0:dt:tf)-1)
        k
        figure(1)
        plot_UAV_GMM(c0,d,xlim,prior_GMM,sts,r,Ntar)
        hold on;
        for sdsd=1:1:Ns
            plot(XY{sdsd}(:,1),XY{sdsd}(:,2),coll{sdsd},'linewidth',2)
        end
        box on
        hold off;
        %             plot_prop_paper
        xlabel('x')
        ylabel('y')
        %             keyboard
        pause(0.5)
        axis square
                     saveas(gcf,strcat('UAV2_move_targs_range2',num2str(tt-1)), 'pdf')
                     saveas(gcf,strcat('UAV2_move_targs_range2',num2str(tt-1)), 'fig')
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
            
            %     [txx,xx1]=ode45_discc(targ_dyn{nt},t0,dt,dt,prior_GMM.truth{nt},zeros(prior_GMM.fn{i}));
%             keyboard
            prior_GMM.truth{nt}=vertcat(prior_GMM.truth{nt},targ_dyn{nt}(prior_GMM.truth{nt}(end,:)',para{nt})');
        end
        % update with meas
        for nt=1:1:Ntar
            z=[];
            true_pos=prior_GMM.truth{nt}(end,:);
            SS=[];
            for ss=1:1:Ns
                pp=[Xgrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2))];
                if norm(true_pos(sts)-pp)<r
                    sqrt(R{ss})*randn(nz{ss},1)
                    z=vertcat(z,h{ss}(true_pos,pp)+sqrt(R{ss})*randn(nz{ss},1));
                    SS=vertcat(SS,[ss,pp]);% store the sensors that were used
                end
            end
            if isempty(z)==0
%                 if length(z)>1
%                     keyboard
%                 end
            prGMM1.mu=prior_GMM.mu{nt};prGMM1.P=prior_GMM.P{nt};prGMM1.w=prior_GMM.w{nt};
            measGMM=post_GMM2_modf_meas(prGMM1,z,h,G,R,SS);
            prior_GMM.mu{nt}=measGMM.mu;prior_GMM.P{nt}=measGMM.P;prior_GMM.w{nt}=measGMM.w;
%             prGMM1.mu=prior_GMM.mu{nt};prGMM1.P=prior_GMM.P{nt};prGMM1.w=prior_GMM.w{nt};
            end
        end
        Hist.u{tt}=U;
        Hist.priorGMM{tt}=prior_GMM;
        tt=tt+1;
    end % implement first 10 steps
%      keyboard
end %cont loop while
% save('control_UAV','Hist','targ_dyn','g','h','G','R','Q','d','t0','dt','tf','para','nz','xlim','Ntar','Ns','Nu','para','r','s0')

%% plotiing the trajectory along with
%  for i=1:1:25
%     openfig(strcat('UAV2_stat_targs_bear',num2str(i),'.fig')); 
%      saveas(gcf,strcat('UAV2_stat_targs_bear',num2str(i),'.jpg'))
%  end