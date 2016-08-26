function Kirb_sims_12thJuly2012_allfilters_call()
clear all
clc 
close all

     matlabpool open 10
%% ------------------------------------------------------------------------
% time
%--------------------------------------------------------------------------
time.t0 = 0;
time.dt = 5;
time.tf = 125+90+125+30+125;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);
%__________________________________________________________________________

%% ------------------------------------------------------------------------
% model
T=time.dt;
L1=0.2;
L2=0.01;
Q_CT=L1*[T^3/3,0,T^2/2,0,0;
        0,T^3/3,0,T^2/2,0;
        T^2/2,0,T,0,0;
        0,T^2/2,0,T,0;
        0,0,0,0,T*L2/L1];
%     Q_UM=Q_CT(1:4,1:4);
 R=diag([(7e2)^2,(7e2)^2]);    
     
model.fn = 5;               % state space dimensionality
model.fx = @KIRB_CT_eg_dyn_disc;
model.fx_jac=@KIRB_CT_eg_dyn_jac_disc;
model.hn = 2;               % measurement dimensionality
model.hx =@KIRB_eg_meas_disc;
model.hx_jac=@KIRB_eg_meas_jac_disc;
model.Q = Q_CT;
model.sQ = sqrtm(model.Q);
model.R = R;
model.sR=sqrtm(R);
model.Qtruth=1e-200;
model.x0tr=[25000,10000,-120,0]';
model.para_dt(1)=T;

%% ------------------------------------------------------------------------
% particle filter settings
%--------------------------------------------------------------------------
pf.no_particles = 5000;
pf.no_bins = 1000;               % when computing the histogram

pf.resample = true;             % resampling after measurement update
pf.neff = 0.8;                  % treshold for the number of effective samples

%% --------------------------------------------------------------
%filter props
x0tr_CT=[25000,10000,-120,0,0.5]';
P0tr_CT=diag([1000^2,1000^2,100,100,(2*pi/180)^2]);
filter.paras_ukf_kappa=1;
filter.paras_gh_pts=4;
filter.freq=1; %% This is actually the number of 'dt' steps after which a meas updt is done
filter.P0_filt_start=P0tr_CT;
filter.save_pf_data='false';
%switch on filters
filter.EKF='true';
filter.CKF='true';
filter.UKF='true';
filter.KF='false';
filter.GHKF='true';
filter.PF='true';
filter.CUT4KF='true';
filter.CUT6KF='true';
filter.CUT8KF='true';
filter.truth=[];
filter.ymeas=[];
%% Allocate space for NNN no. of iterations
NNN=50;

xNNN_mc=zeros(time.nSteps,model.fn,NNN);
YNNN_mc=zeros(time.nSteps,model.hn,NNN);
xNNN_ukf=zeros(time.nSteps,model.fn,NNN);
xNNN_ckf=zeros(time.nSteps,model.fn,NNN);
xNNN_cut4=zeros(time.nSteps,model.fn,NNN);
xNNN_cut6=zeros(time.nSteps,model.fn,NNN);
xNNN_cut8=zeros(time.nSteps,model.fn,NNN);
xNNN_gh=zeros(time.nSteps,model.fn,NNN);
xNNN_ekf=zeros(time.nSteps,model.fn,NNN);
xNNN_mupf=zeros(time.nSteps,model.fn,NNN);
xNNN_mopf=zeros(time.nSteps,model.fn,NNN);
PNNN_ukf=zeros(time.nSteps,model.fn^2,NNN);
PNNN_ckf=zeros(time.nSteps,model.fn^2,NNN);
PNNN_cut4=zeros(time.nSteps,model.fn^2,NNN);
PNNN_cut6=zeros(time.nSteps,model.fn^2,NNN);
PNNN_cut8=zeros(time.nSteps,model.fn^2,NNN);
PNNN_gh=zeros(time.nSteps,model.fn^2,NNN);
PNNN_ekf=zeros(time.nSteps,model.fn^2,NNN);
PNNN_pf=zeros(time.nSteps,model.fn^2,NNN);
PF_complete_data=cell([1,NNN]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the NNN data sets in a cell to sent to parfor lloop

    % firstly generating the truth
     [t,x_mc1]=ode45_discc(@(x)KIRB_UM_eg_dyn_disc(x,model.para_dt),time.t0,time.dt,125,model.x0tr,1e-200);
    [t,x_mc2]=ode45_discc(@(x)KIRB_CT_eg_dyn_disc(x,model.para_dt),t(end),time.dt,t(end)+90,[x_mc1(end,:)';1*pi/180],1e-200);
    [t,x_mc3]=ode45_discc(@(x)KIRB_UM_eg_dyn_disc(x,model.para_dt),t(end),time.dt,t(end)+125,x_mc2(end,1:end-1)',1e-200);
    [t,x_mc4]=ode45_discc(@(x)KIRB_CT_eg_dyn_disc(x,model.para_dt),t(end),time.dt,t(end)+30,[x_mc3(end,:)';-3*pi/180],1e-200);
    [t,x_mc5]=ode45_discc(@(x)KIRB_UM_eg_dyn_disc(x,model.para_dt),t(end),time.dt,t(end)+125,x_mc4(end,1:end-1)',1e-200);
    x_mc=[[x_mc1(1:end-1,1:4),zeros(size(x_mc1,1)-1,1)];x_mc2(1:end-1,1:5);[x_mc3(1:end-1,1:4),zeros(size(x_mc3,1)-1,1)];x_mc4(1:end-1,1:5);[x_mc5(1:end,1:4),zeros(size(x_mc5,1),1)]];

    data_set=cell([4,NNN]);
    
    %Start from a random intiial point and random measuremnts generated
    % from the truth.
for i=1:1:NNN
    filter.x0_filt_start=mvnrnd(x0tr_CT,P0tr_CT)';
    
    %generating the measurement
    ym=zeros(size(x_mc,1),model.hn);
    for ii=1:1:time.nSteps
        ym(ii,:)=(model.hx(x_mc(ii,:)',model.para_dt)+model.sR*randn(model.hn,1))';
    end
%         figure(1)
% plot(x_mc(:,1),x_mc(:,2),ym(:,1),ym(:,2),'k*')
    filter.ymeas=ym;
    filter.truth=x_mc;
    at = now;
aa=datestr(at);
aa(find(aa==' '))='_';
aa(find(aa==':'))='-';
    model.name=strcat('KIRB1_',aa,'_NNN_',num2str(i));
    
    data_set{1,i}=model;
    data_set{2,i}=filter;
    data_set{3,i}=pf;
    data_set{4,i}=time;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run in parallel
parfor j=1:NNN

 j
results_filter=Allfilters(data_set{1,j},data_set{2,j},data_set{3,j},data_set{4,j});

  xNNN_mc(:,:,j)=results_filter.mc;
  YNNN_mc(:,:,j)=results_filter.ymeas;
  xNNN_mopf(:,:,j)=results_filter.mopf;
  xNNN_mupf(:,:,j)=results_filter.mupf;
  PNNN_pf(:,:,j)=results_filter.Ppf;
  xNNN_gh(:,:,j)=results_filter.mu_gh;
  xNNN_ekf(:,:,j)=results_filter.mu_ekf;
  xNNN_cut8(:,:,j)=results_filter.mu_cut8;
  xNNN_cut6(:,:,j)=results_filter.mu_cut6;
  xNNN_cut4(:,:,j)=results_filter.mu_cut4;
  xNNN_ckf(:,:,j)=results_filter.mu_ckf;
  xNNN_ukf(:,:,j)=results_filter.mu_ukf;
  PNNN_gh(:,:,j)=results_filter.P_gh;
  PNNN_ekf(:,:,j)=results_filter.P_ekf;
  PNNN_cut8(:,:,j)=results_filter.P_cut8;
  PNNN_cut6(:,:,j)=results_filter.P_cut6;
  PNNN_cut4(:,:,j)=results_filter.P_cut4;
  PNNN_ckf(:,:,j)=results_filter.P_ckf;
  PNNN_ukf(:,:,j)=results_filter.P_ukf; 
  PF_complete_data{j}=results_filter.PFdata;
end

% Averaging over all the runs
est_fin_ukf=zeros(time.nSteps,3);
est_fin_ckf=zeros(time.nSteps,3);
est_fin_cut4=zeros(time.nSteps,3);
est_fin_cut6=zeros(time.nSteps,3);
est_fin_cut8=zeros(time.nSteps,3);
est_fin_gh=zeros(time.nSteps,3);
est_fin_ekf=zeros(time.nSteps,3);
est_fin_mupf=zeros(time.nSteps,3);
est_fin_mopf=zeros(time.nSteps,3);
avg_traj_ukf=zeros(time.nSteps,model.fn);
avg_traj_ckf=zeros(time.nSteps,model.fn);
avg_traj_cut4=zeros(time.nSteps,model.fn);
avg_traj_cut6=zeros(time.nSteps,model.fn);
avg_traj_cut8=zeros(time.nSteps,model.fn);
avg_traj_gh=zeros(time.nSteps,model.fn);
avg_traj_ekf=zeros(time.nSteps,model.fn);
avg_traj_mupf=zeros(time.nSteps,model.fn);
avg_traj_mopf=zeros(time.nSteps,model.fn);
avg_traj_mc=zeros(time.nSteps,model.fn);
avg_traj_meas=zeros(time.nSteps,model.hn);

for r=1:1:time.nSteps
        for c=1:1:3

           if c==1
           est_fin_ukf(r,c)= sqrt(mean((xNNN_ukf(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_ukf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
           est_fin_ckf(r,c) =sqrt(mean((xNNN_ckf(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_ckf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
           est_fin_cut4(r,c)=sqrt(mean((xNNN_cut4(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_cut4(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
           est_fin_cut6(r,c)=sqrt(mean((xNNN_cut6(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_cut6(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
           est_fin_cut8(r,c)=sqrt(mean((xNNN_cut8(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_cut8(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
           est_fin_gh(r,c)=  sqrt(mean((xNNN_gh(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_gh(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
           est_fin_ekf(r,c)=  sqrt(mean((xNNN_ekf(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_ekf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
           est_fin_mupf(r,c)=  sqrt(mean((xNNN_mupf(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_mupf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
           est_fin_mopf(r,c)=  sqrt(mean((xNNN_mopf(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_mopf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
           end
           if c==2
           est_fin_ukf(r,c)= sqrt(mean((xNNN_ukf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_ukf(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
           est_fin_ckf(r,c) =sqrt(mean((xNNN_ckf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_ckf(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
           est_fin_cut4(r,c)=sqrt(mean((xNNN_cut4(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_cut4(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
           est_fin_cut6(r,c)=sqrt(mean((xNNN_cut6(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_cut6(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
           est_fin_cut8(r,c)=sqrt(mean((xNNN_cut8(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_cut8(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
           est_fin_gh(r,c)=  sqrt(mean((xNNN_gh(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_gh(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
            est_fin_ekf(r,c)=  sqrt(mean((xNNN_ekf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_ekf(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
           est_fin_mupf(r,c)=  sqrt(mean((xNNN_mupf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_mupf(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
           est_fin_mopf(r,c)=  sqrt(mean((xNNN_mopf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_mopf(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
           end 
           if c==3
           est_fin_ukf(r,c)= sqrt(mean((xNNN_ukf(r,5,:)-xNNN_mc(r,5,:)).^2));
           est_fin_ckf(r,c) =sqrt(mean((xNNN_ckf(r,5,:)-xNNN_mc(r,5,:)).^2));
           est_fin_cut4(r,c)=sqrt(mean((xNNN_cut4(r,5,:)-xNNN_mc(r,5,:)).^2));
           est_fin_cut6(r,c)=sqrt(mean((xNNN_cut6(r,5,:)-xNNN_mc(r,5,:)).^2));
           est_fin_cut8(r,c)=sqrt(mean((xNNN_cut8(r,5,:)-xNNN_mc(r,5,:)).^2));
           est_fin_gh(r,c)=  sqrt(mean((xNNN_gh(r,5,:)-xNNN_mc(r,5,:)).^2));
           est_fin_ekf(r,c)=  sqrt(mean((xNNN_ekf(r,5,:)-xNNN_mc(r,5,:)).^2));
           est_fin_mupf(r,c)=  sqrt(mean((xNNN_mupf(r,5,:)-xNNN_mc(r,5,:)).^2));
           est_fin_mopf(r,c)=  sqrt(mean((xNNN_mopf(r,5,:)-xNNN_mc(r,5,:)).^2));
           end
        end
end
for r=1:1:time.nSteps
        for c=1:1:model.fn
avg_traj_ukf(r,c)=mean(xNNN_ukf(r,c,:));
avg_traj_ckf(r,c)=mean(xNNN_ckf(r,c,:));
avg_traj_cut4(r,c)=mean(xNNN_cut4(r,c,:));
avg_traj_cut6(r,c)=mean(xNNN_cut6(r,c,:));
avg_traj_cut8(r,c)=mean(xNNN_cut8(r,c,:));
avg_traj_gh(r,c)=mean(xNNN_gh(r,c,:));
avg_traj_ekf(r,c)=mean(xNNN_ekf(r,c,:));
avg_traj_mupf(r,c)=mean(xNNN_mupf(r,c,:));
avg_traj_mopf(r,c)=mean(xNNN_mopf(r,c,:));
avg_traj_mc(r,c)=mean(xNNN_mc(r,c,:));
if c<=model.hn
avg_traj_meas(r,c)=mean(YNNN_mc(r,c,:));
end
        end
end

save(model.name,'model','time','pf','filter','PF_complete_data','avg_traj_ukf','avg_traj_ckf','avg_traj_cut4','avg_traj_cut6','avg_traj_cut8','avg_traj_gh','avg_traj_ekf','avg_traj_mopf','avg_traj_mupf','avg_traj_mc','est_fin_ukf','est_fin_ckf','est_fin_cut4','est_fin_cut6','est_fin_cut8','est_fin_gh','est_fin_ekf','est_fin_mupf','est_fin_mopf','xNNN_mc','YNNN_mc','xNNN_ukf','xNNN_ckf','xNNN_cut4','xNNN_cut6','xNNN_cut8','xNNN_gh','xNNN_ekf','xNNN_mupf','xNNN_mopf','PNNN_ukf','PNNN_ckf','PNNN_cut4','PNNN_cut6','PNNN_cut8','PNNN_gh','PNNN_ekf','PNNN_pf')
  matlabpool close
%% Plotting
% t=time.tspan;
% figure(1)
% plot(t,est_fin_mupf(:,1),'-.',t,est_fin_mopf(:,1),':',t,est_fin_ckf(:,1),t,est_fin_ukf(:,1),t,est_fin_cut4(:,1),t,est_fin_cut6(:,1),t,est_fin_cut8(:,1),t,est_fin_gh(:,1),t,est_fin_ekf(:,1),'LineWidth',2)
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh','ekf')
% xlabel('t (s)')
% ylabel('RMSE_{pos} (m)')
% plot_prop_paper
% 
% figure(2)
% plot(t,est_fin_mupf(:,2),'-.',t,est_fin_mopf(:,2),':',t,est_fin_ckf(:,2),t,est_fin_ukf(:,2),t,est_fin_cut4(:,2),t,est_fin_cut6(:,2),t,est_fin_cut8(:,2),t,est_fin_gh(:,2),t,est_fin_ekf(:,2),'LineWidth',2)
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh','ekf')
% xlabel('t (s)')
% ylabel('RMSE_{vel} (m/s)')
% plot_prop_paper
% 
% figure(3)
% plot(t,(180/pi)*est_fin_mupf(:,3),'-.',t,(180/pi)*est_fin_mopf(:,3),':',t,(180/pi)*est_fin_ckf(:,3),t,(180/pi)*est_fin_ukf(:,3),t,(180/pi)*est_fin_cut4(:,3),t,(180/pi)*est_fin_cut6(:,3),t,(180/pi)*est_fin_cut8(:,3),t,(180/pi)*est_fin_gh(:,3),t,(180/pi)*est_fin_ekf(:,3),'LineWidth',2)
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh','ekf')
% xlabel('t (s)')
% ylabel('RMSE_{\Omega} (deg/s)')
% plot_prop_paper
% 
% figure(4)
% plot(avg_traj_mupf(:,1),avg_traj_mupf(:,2),'-.',avg_traj_mopf(:,1),avg_traj_mopf(:,2),':',avg_traj_ckf(:,1),avg_traj_ckf(:,2),avg_traj_ukf(:,1),avg_traj_ukf(:,2),avg_traj_cut4(:,1),avg_traj_cut4(:,2),avg_traj_cut6(:,1),avg_traj_cut6(:,2),avg_traj_cut8(:,1),avg_traj_cut8(:,2,end),avg_traj_gh(:,1),avg_traj_gh(:,2),avg_traj_ekf(:,1),avg_traj_ekf(:,2),avg_traj_mc(:,1),avg_traj_mc(:,2),'--',YNNN_mc(:,1),YNNN_mc(:,2),'k*','LineWidth',2)
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh','ekf','MC','sensor')
% xlabel('x (m)')
% ylabel('y (m)')
% plot_prop_paper
%% Compute 2-norm over time Numbers

fid = fopen(strcat(model.name,'.txt'), 'wt');
fprintf(fid, '%s % f % f % f % f % f\n' ,'Q=', model.Q);
fprintf(fid, '\n')
fprintf(fid, '%s % f %f\n' ,'R=', model.R);
fprintf(fid, '\n')
fprintf(fid, '%s % f' ,'GH=', filter.paras_gh_pts);
fprintf(fid, '\n')
fprintf(fid, '%s % f' ,'PF=', pf.no_particles);
fprintf(fid, '\n')
fprintf(fid, '%s % f % f % f % f % f' ,'x0tr=', x0tr_CT');
fprintf(fid, '\n')
fprintf(fid, '%s % f % f % f % f % f\n' ,'P0tr=', P0tr_CT);
fprintf(fid, '\n')

fprintf(fid, '\n')
fprintf(fid, ' %s %s %s %s %s %s %s %s %s\n' , ['PF-MEAN  ','    PF-MODE ','     CKF ','      UKF ','       CUT4 ','          CUT6 ','       CUT8 ','       GH ','       EKF ']);
fprintf(fid, '\n')
% in position
fprintf(fid, '%s','Position')
fprintf(fid, '\n')
aa=1/1000*(sqrt([sum((est_fin_mopf(:,1)-est_fin_mupf(:,1)).^2),sum((est_fin_ckf(:,1)-est_fin_mupf(:,1)).^2),sum((est_fin_ukf(:,1)-est_fin_mupf(:,1)).^2),sum((est_fin_cut4(:,1)-est_fin_mupf(:,1)).^2),sum((est_fin_cut6(:,1)-est_fin_mupf(:,1)).^2),sum((est_fin_cut8(:,1)-est_fin_mupf(:,1)).^2),sum((est_fin_gh(:,1)-est_fin_mupf(:,1)).^2),sum((est_fin_ekf(:,1)-est_fin_mupf(:,1)).^2)]));
bb=sqrt([sum(est_fin_mupf(:,1).^2),sum(est_fin_mopf(:,1).^2),sum(est_fin_ckf(:,1).^2),sum(est_fin_ukf(:,1).^2),sum(est_fin_cut4(:,1).^2),sum(est_fin_cut6(:,1).^2),sum(est_fin_cut8(:,1).^2),sum(est_fin_gh(:,1).^2),sum(est_fin_ekf(:,1).^2)]);
fprintf(fid, '% 12.8f' , bb);
fprintf(fid, '\n')
% in vel
fprintf(fid, '%s','Velocity')
fprintf(fid, '\n')
aa=1/1000*(sqrt([sum((est_fin_mopf(:,2)-est_fin_mupf(:,2)).^2),sum((est_fin_ckf(:,2)-est_fin_mupf(:,2)).^2),sum((est_fin_ukf(:,2)-est_fin_mupf(:,2)).^2),sum((est_fin_cut4(:,2)-est_fin_mupf(:,2)).^2),sum((est_fin_cut6(:,2)-est_fin_mupf(:,2)).^2),sum((est_fin_cut8(:,2)-est_fin_mupf(:,2)).^2),sum((est_fin_gh(:,2)-est_fin_mupf(:,2)).^2),sum((est_fin_ekf(:,2)-est_fin_mupf(:,2)).^2)]));
bb=sqrt([sum(est_fin_mupf(:,2).^2),sum(est_fin_mopf(:,2).^2),sum(est_fin_ckf(:,2).^2),sum(est_fin_ukf(:,2).^2),sum(est_fin_cut4(:,2).^2),sum(est_fin_cut6(:,2).^2),sum(est_fin_cut8(:,2).^2),sum(est_fin_gh(:,2).^2),sum(est_fin_ekf(:,2).^2)]);
fprintf(fid, '% 12.8f' , bb);
fprintf(fid, '\n')
% in omega
fprintf(fid, '%s','Omega')
fprintf(fid, '\n')
aa=(sqrt([sum((est_fin_mopf(:,3)-est_fin_mupf(:,3)).^2),sum((est_fin_ckf(:,3)-est_fin_mupf(:,3)).^2),sum((est_fin_ukf(:,3)-est_fin_mupf(:,3)).^2),sum((est_fin_cut4(:,3)-est_fin_mupf(:,3)).^2),sum((est_fin_cut6(:,3)-est_fin_mupf(:,3)).^2),sum((est_fin_cut8(:,3)-est_fin_mupf(:,3)).^2),sum((est_fin_gh(:,3)-est_fin_mupf(:,3)).^2),sum((est_fin_ekf(:,3)-est_fin_mupf(:,3)).^2)]));
bb=sqrt([sum(est_fin_mupf(:,3).^2),sum(est_fin_mopf(:,3).^2),sum(est_fin_ckf(:,3).^2),sum(est_fin_ukf(:,3).^2),sum(est_fin_cut4(:,3).^2),sum(est_fin_cut6(:,3).^2),sum(est_fin_cut8(:,3).^2),sum(est_fin_gh(:,3).^2),sum(est_fin_ekf(:,3).^2)])
fprintf(fid, '% 12.8f' , bb);
fprintf(fid, '\n')
fclose(fid);
end