%% Using contextual information
%  by psedo measurements

load('0_0_1p5708s40s-1s1s8s18s5s0s0')


%% first generate the arena
% figure
% hold on
% axis([0 10 0 10])
% [xu,yu] = ginput(7);
% plot(xu,yu)
% axis([0 10 0 10])
% [xd,yd] = ginput(7);
% plot(xd,yd)
% axis([0 10 0 10])
% [xcu,ycu] = ginput(7);
% plot(xcu,ycu)
% axis([0 10 0 10])
% [xcd,ycd] = ginput(7);
% plot(xcd,ycd)
% axis([0 10 0 10])

% save('EgNaval','xu','yu','xd','yd','xcu','ycu','xcd','ycd')
%%

load('EgNaval')


plot(xu,yu,'k',xd,yd,'k',xcu,ycu,'k',xcd,ycd,'k')
hold on
dx=0.2;
Xc=0:dx:10;
Yc=0:dx:10;


for i=1:1:length(Xc)
    plot([Xc(i),Xc(i)],[0,10],'r')
    plot([0,10],[Yc(i),Yc(i)],'r')
end

Xg=dx/2:dx:10-dx/2;
Yg=dx/2:dx:10-dx/2;
k=1;
XY=zeros(length(Xg)^2,5);
for i=1:1:length(Xg)
    for j=1:1:length(Xg)
        XY(k,:)=[Xg(i),Yg(j),i,j,0.7];
        k=k+1;
    end
end



% plot(XY(:,1),XY(:,2),'b+')
%% generate the truth
%
% [xt5,yt5] = ginput(15);
% plot(xt5,yt5)
% save('Xtruth','xt1','yt1')
% Xtruth{5}=[xt5,yt5];
%
% save('Xtruth','Xtruth')
load('Xtruth')
for i=1:1:length(Xtruth)
    Xtruth{i}=Xtruth{i}*40;
end

% DP=load('Xtruth2_good','Xtruth');
% Xtruth=DP.Xtruth;
% clear DP;
%%  Traficability values
% XY(:,3)=5*ones(size(XY,1),1); % default value in the regions where u can go
Pu=[xu,yu;10,10;0,10;xu(1),yu(1)];
Pd=[xd,yd;10,0;xd(1),yd(1)];
Pc=[xcu,ycu;xcd(end:-1:1),ycd(end:-1:1);0,0;xcu(1),ycu(1)];

dx=dx*40;
XY(:,1:2)=40*XY(:,1:2);
Xg=Xg*40; Yg=Yg*40;
xu=xu*40; yu=yu*40;
xd=xd*40; yd=yd*40;
xcu=xcu*40; ycu=ycu*40;
xcd=xcd*40; ycd=ycd*40;
Pu=Pu*40;
Pd=Pd*40;
Pc=Pc*40;

figure(2)
plot(Pu(:,1),Pu(:,2),Pd(:,1),Pd(:,2),Pc(:,1),Pc(:,2))

N=size(XY,1);

[in,on] = inpolygon(XY(:,1),XY(:,2),Pu(:,1),Pu(:,2));
XY(find((in+on)==1),5)=0.01;

[in,on] = inpolygon(XY(:,1),XY(:,2),Pd(:,1),Pd(:,2));
XY(find((in+on)==1),5)=0.01;

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
%% Time props for the targets
dt=5;
Time.t0=0;
Time.tf=500;
Time.Tvec=Time.t0:dt:Time.tf;
Time.dt=dt;
Time.nt=length(Time.Tvec);

%% target positions
TargetPos.truth=cell(length(Xtruth),1);
for i=1:1:floor(length(Xtruth)/2)
    TargetPos.truth{i} =[ interp1(linspace(Time.t0,Time.tf,length(Xtruth{i}(:,1))),Xtruth{i}(:,1),Time.Tvec)',...
        interp1(linspace(Time.t0,Time.tf,length(Xtruth{i}(:,1))),Xtruth{i}(:,2),Time.Tvec)'];
end
for i=floor(length(Xtruth)/2)+1:1:length(Xtruth)
    TargetPos.truth{i} =[ interp1(linspace(Time.t0,Time.tf,length(Xtruth{i}(:,1))),Xtruth{i}(end:-1:1,1),Time.Tvec)',...
        interp1(linspace(Time.t0,Time.tf,length(Xtruth{i}(:,1))),Xtruth{i}(end:-1:1,2),Time.Tvec)'];
end

TargetPos.Ntar=length(TargetPos.truth);

% traffic values
TargetPos.Traffic=XY;
TargetPos.Traffic_xlim=400;
TargetPos.Traffic_ylim=400;
TargetPos.Trafficdx=dx;
TargetPos.Pu=Pu;
TargetPos.Pd=Pd;
TargetPos.Pc=Pc;

%Target initial filter positions


for i=1:1:TargetPos.Ntar
    TargetPos.fn(i)=4;
    
    TargetPos.xf{i}=zeros(Time.nt,TargetPos.fn(i));
    TargetPos.Pf{i}=zeros(Time.nt,TargetPos.fn(i)^2);
    TargetPos.xf{i}(1,:)=[TargetPos.truth{i}(1,1),TargetPos.truth{i}(1,2),0.1*(TargetPos.truth{i}(2,1)-TargetPos.truth{i}(1,1))/dt,0.1*(TargetPos.truth{i}(2,2)-TargetPos.truth{i}(1,2))/dt];
    TargetPos.Pf{i}(1,:)=reshape(diag([5^2,5^2,(0.5)^2,(0.5)^2]),1,16);
    TargetPos.Q{i}=diag([1^2,1^2,(0.1)^2,(0.1)^2]);
    
end

TargetPos.thresh=0.01;


%% Sensors and generating measuremetns
SensorModel.Move.Ns=2;

Trajectories=cell(SensorModel.Move.Ns,1);


%UAV motion grid

[Ygrid,Xgrid]=meshgrid(0:40:400);

SensorModel.Move.Griddx=40;
SensorModel.Move.Gridx0=0;
SensorModel.Move.Gridxlim=400;
SensorModel.Move.Gridy0=0;
SensorModel.Move.Gridylim=400;

SensorModel.Move.GridTemplate=templates1;
SensorModel.Move.GridTemplate2=templates2;


SensorModel.Move.Xgrid=Xgrid;
SensorModel.Move.Ygrid=Ygrid;

SensorModel.Move.Grid2Abs=@(i,j,SensorModel)[SensorModel.Move.Xgrid(i,j),SensorModel.Move.Ygrid(i,j)];
% SensorModel.Move.Abs2Grid=@(x,y,SensorModel)[SensorModel.Move.Xgrid(i,j),SensorModel.Move.Ygrid(i,j)];


SensorModel.Move.GridU=9; % the number of grid-grid contrls available
SensorModel.Move.GridNdir=4; % direction state space 1-north/up, 2-right/east

SensorModel.Move.Gridpos{1}=zeros(Time.nt,3);
SensorModel.Move.Gridpos{2}=zeros(Time.nt,3);
SensorModel.Move.Gridpos{1}(1,:)=[2,3,1];
SensorModel.Move.Gridpos{2}(1,:)=[6,8,1];
Trajectories{1}=[SensorModel.Move.Grid2Abs(SensorModel.Move.Gridpos{1}(1,1),SensorModel.Move.Gridpos{1}(1,2),SensorModel),pi/2,1];
Trajectories{2}=[SensorModel.Move.Grid2Abs(SensorModel.Move.Gridpos{2}(1,1),SensorModel.Move.Gridpos{2}(1,2),SensorModel),pi/2,1];

SensorModel.Move.GridMotion=@(ix,jy,dir,u)UAV_grid_motion(ix,jy,dir,u,SensorModel);

SensorModel.Move.GridOptiControl{1}=zeros(Time.nt,1);
SensorModel.Move.GridOptiControl{2}=zeros(Time.nt,1);



SensorModel.Move.FOV{1}=[pi,80,pi];
SensorModel.Move.FOV{2}=[pi,70,pi];

SensorModel.Move.trajcol{1}='r';
SensorModel.Move.trajcol{2}='b';

SensorModel.Move.hx{1}=@(xtarg,xsenspos,FOV)hxrb(xtarg,xsenspos,FOV); % xtarg,SensorModel.Move.pos{1},SensorModel.Move.FOV{1}
SensorModel.Move.hx{2}=@(xtarg,xsenspos,FOV)hxrb(xtarg,xsenspos,FOV);

% SensorModel.Move.R{1}=(0.2*pi/180)^2;
% SensorModel.Move.R{2}=(0.2*pi/180)^2;
SensorModel.Move.R{1}=diag([5^2,(2*pi/180)^2]);
SensorModel.Move.R{2}=diag([5^2,(2*pi/180)^2]);

SensorModel.Move.hn(1)=2;
SensorModel.Move.hn(2)=2;

SensorModel.Move.Meas{1}=zeros(Time.nt,SensorModel.Move.hn(1));
SensorModel.Move.Meas{2}=zeros(Time.nt,SensorModel.Move.hn(2));


% Static sensors
SensorModel.Static.hn(1)=2;
SensorModel.Static.hn(2)=2;

SensorModel.Static.pos{1}=[3.25,2.75]*40;
SensorModel.Static.pos{2}=[6.25,8.75]*40;

SensorModel.Static.FOV{1}=[pi/8,200,pi/2]; %[half angle of FOV, max radius,dir of view meas from +ve x axis := (-pi,pi)]
SensorModel.Static.FOV{2}=[pi/8,300,-pi/2];

SensorModel.Static.hx{1}=@(xtarg,xsenspos,FOV)hxrb(xtarg,xsenspos,FOV);  %(xtarg,SensorModel.Static.pos{1},FOV)
SensorModel.Static.hx{2}=@(xtarg,xsenspos,FOV)hxrb(xtarg,xsenspos,FOV);

SensorModel.Static.R{1}=diag([2^2,(1*pi/180)^2]);
SensorModel.Static.R{2}=diag([2^2,(1*pi/180)^2]);

SensorModel.Static.Ns=length(SensorModel.Static.pos);

SensorModel.Static.Meas{1}=zeros(Time.nt,SensorModel.Static.hn(1));
SensorModel.Static.Meas{2}=zeros(Time.nt,SensorModel.Static.hn(2));

SensorModel.Static.Task=cell(Time.nt,1);  % task [sensor,targ] then make the measurement
SensorModel.Move.Task=cell(Time.nt,1);    % task [sensor,targ] then make the measurement


SensorModel.InfoMap=cell(Time.nt,1);
SensorModel.Move.DP_MU=cell(SensorModel.Move.Ns,Time.nt);
for k=1:1:Time.nt
    SensorModel.InfoMap{k}=zeros(size(SensorModel.Move.Xgrid));
end





%% Running the filtering algorithm

TargetPos.filtermethod='ut';
SensorModel.filtermethod='ut';
TargetPos.Dyn=1; % 1 is for normal, 2 is for estimating alpha also
SensorModel.Move.Mode='active';   % passive means no sensor constraints, active means sensors constraints, off = no sensor
SensorModel.Static.Mode='active';


clc
close all
Nk=7;
figure(1)
plot_scene_UAVnaval(TargetPos,Trajectories,SensorModel,1)
hold off
prevk=1;

for k=1:1:Time.nt-1
    

    if k==1
        disp('computing new trajectory')
        [MU,J,SensorModel]=dubin_DP(k,k+Nk,Time,SensorModel,TargetPos);
        disp('Done with computing new trajectory')
        prevk=1;
    elseif k-prevk>=3
        disp('computing new trajectory')
        [MU,J,SensorModel]=dubin_DP(k,k+Nk,Time,SensorModel,TargetPos);
        disp('Done with computing new trajectory')
        prevk=k;
    end
    
    % keyboard
    for ss=1:1:SensorModel.Move.Ns

            [SensorModel,Trajectories]=UAV_traj_pplied_control(ss,SensorModel,Trajectories,k);        
    end
    disp('done with uav motion')
    for ntar=1:1:TargetPos.Ntar
        % Propagrate
        [TargetPos,xf,Pf]=TargetPriorProp(ntar,k,TargetPos,Time);
        
        %SensTasking
        Tsk=SimpleSensortask_static(xf,Pf,SensorModel);
        if isempty(Tsk)==0
            SensorModel.Static.Task{(k+1)}=[Tsk,ntar*ones(length(Tsk),1)];
        end
        Tsk=SimpleSensortask_move([],[],k+1,xf,Pf,SensorModel);
        if isempty(Tsk)==0
            SensorModel.Move.Task{(k+1)}=[Tsk,ntar*ones(length(Tsk),1)];
        end
        %Update
        [TargetPos,xk,Pk]=TargetPostUpdate(ntar,k+1,'true',SensorModel,TargetPos);
        
    end
    disp('done with prop+update')
    figure(1)
    clf
    plot_scene_UAVnaval(TargetPos,Trajectories,SensorModel,k+1)
    hold off;
%     keyboard
        pause(0.1)
%         saveas(gca,strcat('Sim1_UAV_naval_',num2str(k)),'jpg')
        pause(0.1)
    %     if k<Time.nt
    %     clf
    %     end
end
save('SimData_ppt','TargetPos','Trajectories','SensorModel','Time')
close all
nk=10;
figure(3)
for kk=1:1:Time.nt-1
    for ck=5:1:nk
        figure(3)
        clf
    plot_scene_UAVnaval_cont(TargetPos,Trajectories,SensorModel,nk,ck,kk)
%     title(num2str([kk,ck]))
        pause(0.1)
        saveas(gca,strcat('Sim1_UAV_naval_',num2str(kk),'_',num2str(ck)),'jpg')
        pause(0.5)
    end
end