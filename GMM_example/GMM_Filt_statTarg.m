function [GMM,model,sens,ym]=GMM_Filt_statTarg(GMM,model,sens,timeStep,x_mc)
%% GMM contains all the gaussian mixutres
wx=GMM.w;
ng=length(wx);
mux=GMM.mu;
Px=GMM.P;
method=GMM.method;
nx=model.fn;
nz=model.hn;
% kappa=filter.paras_ukf_kappa;
  switch lower(method)
      case 'ut'
          r1=sqrt(nx+GMM.kappa);
          X=[zeros(1,nx);r1*eye(nx);-r1*eye(nx)];
          W=[GMM.kappa/(nx+GMM.kappa);1/(2*(nx+GMM.kappa))*ones(2*nx,1)];
      case 'ckf'    
          X=sqrt(nx)*[eye(nx);-eye(nx)];
          W=1/(2*nx)*ones(2*nx,1);
      case 'gh'
          [X,W]=GH_pts(zeros(nx,1),eye(nx),GMM.ghpts);
      case 'cut4'
          [X,W]=conjugate_dir_gausspts(zeros(nx,1),eye(nx));
      case 'cut6'
          [X,W]=conjugate_dir_gausspts_till_6moment_scheme2(zeros(nx,1),eye(nx));
      case 'cut8'
          [X,W]=conjugate_dir_gausspts_till_8moment(zeros(nx,1),eye(nx));
  end
%% time evolution step
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

%% Measurement Step- Kalman Update of each component
if model.frequpdt(timeStep)==0
    GMM.mu=mux;
    GMM.P=Px;
    GMM.w=wx;
    
    GMM.mean=0;
    GMM.cov=0;
    for i=1:1:ng
        GMM.mean=GMM.mean+w1x(i)*mu1x(i,:)';
        GMM.cov=GMM.cov+w1x(i)*(reshape(P1x(i,:),nx,nx)+mu1x(i,:)'*mu1x(i,:));
    end
    GMM.cov=GMM.cov-GMM.mean*GMM.mean';
return
end

% if GMM.opt==1
% %     Options= optimset('MaxFunEvals',2000);
% %      [phifmin,DKLfmin]=fminunc(@(phi)-ExptKL_modf_MeasOPT_gsm(GMM,model,phi),[model.para_dt(2),model.para_dt(3)],Options);
% %       model.para_dt(2)=phifmin(1);
% %       model.para_dt(3)=phifmin(2);
% elseif GMM.opt==3
%     GMM.mean=0;
%     for i=1:1:ng
%         GMM.mean=GMM.mean+GMM.w(i)*GMM.mu(i,:)';
%     end
%     A=atan2(GMM.mean(2)-sens.xc1(2),GMM.mean(1)-sens.xc1(1));
% 
%     ang=A-20*pi/180:2*pi/180:A+20*pi/180;
%     alp=[2,5,20]*pi/180;
% %     mm=length(ang);
% %     DKL=zeros(1,mm);
% DKLp=0;
%     for phi=ang
%         for alpha=alp
%         DKLn=ExptKL_modf_MeasOPT_gsm(GMM,model,sens,[phi,alpha]);
%         if DKLn>DKLp
%         sens.phi=phi;
%         sens.alpha=alpha;    
%         DKLp=DKLn;
%         end
%         end
%     end
% %     
% %     [aa,bb]=max(DKL);
% %      model.para_dt(2)=ang(bb);
%      
% elseif GMM.opt==2
%     GMM.mean=0;
%     for i=1:1:ng
%         GMM.mean=GMM.mean+GMM.w(i)*GMM.mu(i,:)';
%     end
% 
% %     sens.phi=atan2(GMM.mean(2)-sens.xc1(2),GMM.mean(1)-sens.xc1(1));
% sens.phi=atan2(x_mc(timeStep,2)-sens.xc1(2),x_mc(timeStep,1)-sens.xc1(1));
%     sens.alpha=5*pi/180;
% 
% end
%% generate the measurement with this configuration
 ym=model.hx(x_mc(timeStep,:)',[sens.phi,sens.alpha],sens.xc1)+model.gx(x_mc(timeStep,:)',[sens.phi,sens.alpha],sens.xc1)*mvnrnd(zeros(1,model.hn),model.R, 1)';
ym;
%% Kalmana like update
mu2x=zeros(size(mux));
P2x=zeros(size(Px));
muz=zeros(ng,nz);
Pz=zeros(ng,nz^2);
Pxz=zeros(ng,nx*nz);
for i=1:1:ng
    mu=mux(i,:)';
    sqP=sqrtm(reshape(Px(i,:),nx,nx));
    
    for j=1:1:length(W)
       Xj=sqP*X(j,:)'+mu;
       hx=model.hx(Xj,[sens.phi,sens.alpha],sens.xc1);
       GG=model.gx(x_mc(timeStep,:)',[sens.phi,sens.alpha],sens.xc1)*model.R*model.gx(x_mc(timeStep,:)',[sens.phi,sens.alpha],sens.xc1)';

       muz(i,:)=muz(i,:)+W(j)*hx';
       Pz(i,:)=Pz(i,:)+W(j)*reshape(GG+hx*hx',1,nz^2);
       
       Pxz(i,:)=Pxz(i,:)+W(j)*reshape(Xj*hx',1,nx*nz);
    end
    Pzi=reshape(Pz(i,:),nz,nz)-muz(i,:)'*muz(i,:);
    Pxzi=reshape(Pxz(i,:),nx,nz)-mu*muz(i,:);
    Pz(i,:)=reshape(Pzi,1,nz^2);
    Pxz(i,:)=reshape(Pxzi,1,nx*nz);
    
    Ki=Pxzi/Pzi;
    mu2x(i,:)=mux(i,:)+(Ki*(ym-muz(i,:)'))';
    P2x(i,:)=reshape(reshape(Px(i,:),nx,nx)-Ki*Pzi*Ki',1,nx^2);
end
    GMM.mu=mu2x;
    GMM.P=P2x;
    GMM.w=wx;
    

    GMM.mean=0;
    GMM.cov=0;
 
    for i=1:1:ng
        GMM.mean=GMM.mean+wx(i)*mu2x(i,:)';
        GMM.cov=GMM.cov+wx(i)*(reshape(P2x(i,:),nx,nx)+mu2x(i,:)'*mu2x(i,:));
    end
    GMM.cov=GMM.cov-GMM.mean*GMM.mean';

