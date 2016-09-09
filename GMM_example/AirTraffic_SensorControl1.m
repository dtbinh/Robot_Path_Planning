%  function AirTraffic_SensorControl1()
clear all
clc 
close all

%  matlabpool open 10
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
%  R=diag([(200)^2,(1*pi/180)^2]);    
    R=diag([(1)^2,(1)^2]);     
sens.alpha1=10*pi/180;
sens.alpha2=10*pi/180;
 
  sens.xc1=[-1.6,-2]*1e4;
 sens.xc2=[1,-2]*1e4;
%  sens.alpha=10*pi/180;
 
model.fn = 5;               % state space dimensionality
model.fx = @KIRB_CT_eg_dyn_disc;
model.fx_jac=@KIRB_CT_eg_dyn_jac_disc;
model.hn = 2;               % measurement dimensionality
model.hx =@(x,sens)[sqrt((x(1)-sens.xc1(1))^2+(x(2)-sens.xc1(2))^2);sqrt((x(1)-sens.xc2(1))^2+(x(2)-sens.xc2(2))^2)];
model.gx =@(x,sens)[sensor_radial_penalty(x,sens.xc1,[sens.para(1),sens.para(2)]),0;0,sensor_radial_penalty(x,sens.xc2,[sens.para(3),sens.para(4)])];

% Mutual Info calcualtion props
model.marg_nx=[1,2]; %yes marginalise the first 2 states
model.InfoQuadN=1;
model.InfoMethod='gh';
model.Infoghpts=3;

model.hx_jac=@KIRB_eg_meas_jac_disc2;
model.Q = Q_CT;
model.sQ = sqrtm(model.Q);
model.R = R;
model.sR=sqrtm(R);
model.Qtruth=1e-200;
model.x0tr=[25000,10000,-120,0]';
% model parameters that change with time and sent into the model equations 
% for every evaluation.
model.para_dt(1)=T;
% model.para_dt(2)=pi/4;% second parameter value is sensor orientation 'phi'
model.frequpdt=ones(1,time.nSteps);
%% ------------------------------------------------------------------------
% particle filter settings
%--------------------------------------------------------------------------
pf.no_particles = 5000;
pf.no_bins = 1000;               % when computing the histogram

pf.resample = true;             % resampling after measurement update
pf.neff = 0.8;                  % treshold for the number of effective samples


%% --------------------------------------------------------------
%filter props
x0tr_CT=[25000,10000,-120,0,0.0005]';
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
%% --------------------------------------------------------------------
% Gaussian mixture components
ng=3;
% X0=mvnrnd(x0tr_CT',P0tr_CT,floor( pf.no_particles/4));
% options = statset('MaxIter',5000);
% obj = gmdistribution.fit(X0,GMM.ng,'SharedCov',true,'Options',options);
% GMM.mu=obj.mu;
% GMM.mu=GMM.mu+mvnrnd(x0tr_CT',0.00001*P0tr_CT,GMM.ng);
GMM.mu=[x0tr_CT';mvnrnd(x0tr_CT',P0tr_CT,ng-1)];
GMM.P=repmat(reshape(P0tr_CT,1,model.fn^2),ng,1);
GMM.w=[1,0,0]';
GMM.method='GH';
GMM.kappa=1;
GMM.ghpts=3;
GMM=GMM_moms(GMM);
%use adaptive weights or not
GMM.wtup=1;
%optimize the sensor config 1
GMM.opt=1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the NNN data sets in a cell to sent to parfor lloop

    % firstly generating the truth
    [t,x_mc1]=ode45_discc(@(x)KIRB_UM_eg_dyn_disc(x,model.para_dt),time.t0,time.dt,125,model.x0tr,1e-200);
    [t,x_mc2]=ode45_discc(@(x)KIRB_CT_eg_dyn_disc(x,model.para_dt),t(end),time.dt,t(end)+90,[x_mc1(end,:)';1*pi/180],1e-200);
    [t,x_mc3]=ode45_discc(@(x)KIRB_UM_eg_dyn_disc(x,model.para_dt),t(end),time.dt,t(end)+125,x_mc2(end,1:end-1)',1e-200);
    [t,x_mc4]=ode45_discc(@(x)KIRB_CT_eg_dyn_disc(x,model.para_dt),t(end),time.dt,t(end)+30,[x_mc3(end,:)';-3*pi/180],1e-200);
    [t,x_mc5]=ode45_discc(@(x)KIRB_UM_eg_dyn_disc(x,model.para_dt),t(end),time.dt,t(end)+125,x_mc4(end,1:end-1)',1e-200);
    x_mc=[[x_mc1(1:end-1,1:4),zeros(size(x_mc1,1)-1,1)];x_mc2(1:end-1,1:5);[x_mc3(1:end-1,1:4),zeros(size(x_mc3,1)-1,1)];x_mc4(1:end-1,1:5);[x_mc5(1:end,1:4),zeros(size(x_mc5,1),1)]];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_gmm=zeros(time.nSteps,model.fn);
zk=zeros(time.nSteps,model.hn);
x_gmm(1,:)=GMM.mean';
mu_gh=mvnrnd(x0tr_CT,P0tr_CT,1)';
P_gh=P0tr_CT;
x_gh=zeros(time.nSteps,model.fn);
x_gh(1,:)=mu_gh';
% 
sens.para(1)=atan2(x0tr_CT(2)-sens.xc1(2),x0tr_CT(1)-sens.xc1(1));
sens.para(3)=atan2(x0tr_CT(2)-sens.xc2(2),x0tr_CT(1)-sens.xc2(1));
sens.para(2)=sens.alpha1;
sens.para(4)=sens.alpha2;




plot(x_mc(:,1),x_mc(:,2),'r--',x0tr_CT(1),x0tr_CT(2),'b*-',x0tr_CT(1),x0tr_CT(2),'k+','Linewidth',2)
 plot_sens_view(sens.xc1,sens.para(2),sens.para(1),5*10^4,'b')
plot_sens_view(sens.xc2,sens.para(4),sens.para(3),5*10^4,'r')
 
axis([-1.7,2.8,-2.3,1.5]*10^4)
pause(1)

for j=2:1:time.nSteps
j
   GMM.w    
%        zm=ym(j,:)';
[GMM,model,sens,zm]=GMMF(GMM,model,sens,j,x_mc);
%         [GMM,model,para,zm]=GMMF1(GMM,model,para,j,x_mc);
%          [mu_gh,P_gh]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_gh,P_gh,zm,'gh',filter.paras_gh_pts); 
        zk(j,:)=fsolve(@(z)[(z(1)-sens.xc1(1))^2+(z(2)-sens.xc1(2))^2-zm(1)^2;(z(1)-sens.xc2(1))^2+(z(2)-sens.xc2(2))^2-zm(2)^2],GMM.mean(1:2));

    
    

x_gmm(j,:)=GMM.mean';
x_gh(j,:)=mu_gh';

NN=@(x,mu,P,D)1/sqrt((2*pi)^D*det(P))*exp(-0.5*(x-mu)'*(P\(x-mu))); %Gaussian pdf
 cc= NN(sqrtm(GMM.cov(1:2,1:2))*[1;0]+GMM.mean(1:2),GMM.mean(1:2),GMM.cov(1:2,1:2),2);
 sq=sqrtm(GMM.cov(1:2,1:2));
 [xx,yy]=meshgrid(GMM.mean(1,1)-2*sq(1,1):5:GMM.mean(1,1)+2*sq(1,1),GMM.mean(2,1)-2*sq(2,2):5:GMM.mean(2,1)+2*sq(2,2));

% CC=zeros(size(xx));
% for i=1:1:size(xx,1)
%     for k=1:1:size(xx,2)
%         CC(i,k)=mvnpdf([xx(i,k),yy(i,k)],GMM.mean(1:2)',GMM.cov(1:2,1:2));
% %  CC(i,k)=pdf(obj,[xx(i,k) yy(i,k)]);
%     end
% end
% contour(xx,yy,CC,[cc,cc])
% hold on
plot(x_mc(:,1),x_mc(:,2),'r--',x_gmm(1:j,1),x_gmm(1:j,2),'b*-',zk(1:j,1),zk(1:j,2),'k+','Linewidth',2)
 plot_sens_view(sens.xc1,sens.para(2),sens.para(1),5*10^4,'b')
plot_sens_view(sens.xc2,sens.para(4),sens.para(3),5*10^4,'r')
 
axis([-1.7,2.8,-2.3,1.5]*10^4)
pause(0.1)
end




%  end
 
 %%%%%%%%%%%%%%%
%  function [GMM,model,para,ym]=GMMF1(GMM,model,para,timeStep,x_mc)
% %% GMM contains all the gaussian mixutres
% wx=GMM.w;
% ng=length(wx);
% mux=GMM.mu;
% Px=GMM.P;
% method=GMM.method;
% nx=model.fn;
% nz=model.hn;
% % kappa=filter.paras_ukf_kappa;
%   switch lower(method)
%       case 'ut'
%           r1=sqrt(nx+GMM.kappa);
%           X=[zeros(1,nx);r1*eye(nx);-r1*eye(nx)];
%           W=[GMM.kappa/(nx+GMM.kappa);1/(2*(nx+GMM.kappa))*ones(2*nx,1)];
%       case 'ckf'    
%           X=sqrt(nx)*[eye(nx);-eye(nx)];
%           W=1/(2*nx)*ones(2*nx,1);
%       case 'gh'
%           [X,W]=GH_pts(zeros(nx,1),eye(nx),GMM.ghpts);
%       case 'cut4'
%           [X,W]=conjugate_dir_gausspts(zeros(nx,1),eye(nx));
%       case 'cut6'
%           [X,W]=conjugate_dir_gausspts_till_6moment_scheme2(zeros(nx,1),eye(nx));
%       case 'cut8'
%           [X,W]=conjugate_dir_gausspts_till_8moment(zeros(nx,1),eye(nx));
%   end
% %% time evolution step
% mu1x=zeros(size(mux));
% P1x=zeros(size(Px));
% for i=1:1:ng
%     mu=mux(i,:)';
%     sqP=sqrtm(reshape(Px(i,:),nx,nx));
%     
%     for j=1:1:length(W)
%        Xj=sqP*X(j,:)'+mu;
%        fx=model.fx(Xj,model.para_dt(1));
%        mu1x(i,:)=mu1x(i,:)+W(j)*fx';
%        P1x(i,:)=P1x(i,:)+W(j)*reshape(fx*fx',1,nx^2);
%     end
%     P1x(i,:)=reshape(reshape(P1x(i,:),nx,nx)-mu1x(i,:)'*mu1x(i,:)+model.Q,1,nx^2);
% end
% 
%    
% %% Weight update
% if GMM.wtup==1
%     w1x=wx;
%     GMM.mu=mu1x;
%     GMM.P=P1x;
%     GMM.w=w1x;
% else
%     w1x=wx;
%     GMM.mu=mu1x;
%     GMM.P=P1x;
%     GMM.w=w1x;
% end
% 
% %% Measurement Step- Kalman Update of each component
% if model.frequpdt(timeStep)==0
%     GMM.mu=mu1x;
%     GMM.P=P1x;
%     GMM.w=w1x;
%     
%     GMM.mean=0;
%     GMM.cov=0;
%     for i=1:1:ng
%         GMM.mean=GMM.mean+w1x(i)*mu1x(i,:)';
%         GMM.cov=GMM.cov+w1x(i)*(reshape(P1x(i,:),nx,nx)+mu1x(i,:)'*mu1x(i,:));
%     end
%     GMM.cov=GMM.cov-GMM.mean*GMM.mean';
% return
% end
% %% optimize sensor config
% GMM.opt
% if GMM.opt==1
% 
% ang1=(para(1)-10*pi/180):2*pi/180:(para(1)+10*pi/180);
% ang2=(para(2)-10*pi/180):2*pi/180:(para(2)+10*pi/180);
% DKLp=0;
% i=1;
%     for phi1=ang1
%             for phi2=ang2
%         DKLn=ExptKL_modf_MeasOPT_gsm(GMM,model,[phi1,phi2]);
%         if DKLn>DKLp
%         para(1)=phi1;
%         para(2)=phi2;    
%         DKLp=DKLn;
%         end
%         i=i+1
%             end
%      end
% end
% 
% 
% 
% 
% 
% %% generate the measurement with this configuration
%  ym=model.hx(x_mc(timeStep,:)',para)+model.gx(x_mc(timeStep,:)',para)*mvnrnd(zeros(1,model.hn),model.R, 1)';
% %  ym=[sqrt((x_mc(timeStep,1)-model.xc1(1))^2+(x_mc(timeStep,2)-model.xc1(2))^2),atan2((x_mc(timeStep,2)-model.xc1(2)),(x_mc(timeStep,1)-model.xc1(1)))]';
% % keyboard
% %% Kalmana like update
% mu2x=zeros(size(mux));
% P2x=zeros(size(Px));
% muz=zeros(ng,nz);
% Pz=zeros(ng,nz^2);
% Pxz=zeros(ng,nx*nz);
% for i=1:1:ng
%     mu=mu1x(i,:)';
%     sqP=sqrtm(reshape(P1x(i,:),nx,nx));
%     
%     for j=1:1:length(W)
%        Xj=sqP*X(j,:)'+mu;
%        hx=model.hx(Xj,para);
%        GG=model.gx(Xj,para)*model.R*model.gx(Xj,para)';
% 
%        muz(i,:)=muz(i,:)+W(j)*hx';
%        Pz(i,:)=Pz(i,:)+W(j)*reshape(GG+hx*hx',1,nz^2);
%        
%        Pxz(i,:)=Pxz(i,:)+W(j)*reshape(Xj*hx',1,nx*nz);
%     end
%     Pzi=reshape(Pz(i,:),nz,nz)-muz(i,:)'*muz(i,:);
%     Pxzi=reshape(Pxz(i,:),nx,nz)-mu*muz(i,:);
%     Pz(i,:)=reshape(Pzi,1,nz^2);
%     Pxz(i,:)=reshape(Pxzi,1,nx*nz);
%     
%     Ki=Pxzi/Pzi;
%     mu2x(i,:)=mu1x(i,:)+(Ki*(ym-muz(i,:)'))';
%     P2x(i,:)=reshape(reshape(P1x(i,:),nx,nx)-Ki*Pzi*Ki',1,nx^2);
% end
%     GMM.mu=mu2x;
%     GMM.P=P2x;
%     GMM.w=w1x;
%     
% 
%     GMM.mean=0;
%     GMM.cov=0;
%  
%     for i=1:1:ng
%         GMM.mean=GMM.mean+w1x(i)*mu2x(i,:)';
%         GMM.cov=GMM.cov+w1x(i)*(reshape(P2x(i,:),nx,nx)+mu2x(i,:)'*mu2x(i,:));
%     end
%     GMM.cov=GMM.cov-GMM.mean*GMM.mean';
% 
%  end