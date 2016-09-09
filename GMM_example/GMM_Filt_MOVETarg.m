function [GMM,model,sens]=GMM_Filt_MOVETarg(GMM,model,sens,timeStep)
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
if model.timeupdt==1
mu1x=zeros(size(mux));
P1x=zeros(size(Px));
for i=1:1:ng
    mu=mux(i,:)';
    sqP=sqrtm(reshape(Px(i,:),nx,nx));
    
    for j=1:1:length(W)
       Xj=sqP*X(j,:)'+mu;
       fx=model.fx(Xj,model.para_dt(1));
       mu1x(i,:)=mu1x(i,:)+W(j)*fx';
       P1x(i,:)=P1x(i,:)+W(j)*reshape(fx*fx',1,nx^2);
    end
    P1x(i,:)=reshape(reshape(P1x(i,:),nx,nx)-mu1x(i,:)'*mu1x(i,:)+model.Q,1,nx^2);
end
else
mu1x=mux;
P1x=Px;    
end

   
%% Weight update
if GMM.wtup==1
  NN=@(x,mu,P,D)1/sqrt((2*pi)^D*det(P))*exp(-0.5*(x-mu)'*(P\(x-mu))); %Gaussian pdf
   %evaluating M
   M=zeros(ng,ng);
   N=zeros(ng,ng);
   for i=1:1:ng
       for j=1:1:ng
           
           P=reshape(P1x(i,:)+P1x(j,:),nx,nx);
           M(i,j)=NN(mu1x(i,:)',mu1x(j,:)',P,nx);
           
           mu=mux(j,:)';
           sqP=sqrtm(reshape(Px(j,:),nx,nx));   
           for k=1:1:length(W)
               Xk=sqP*X(k,:)'+mu;
               N(i,j)=N(i,j)+W(k)*NN(model.fx(Xk,model.para_dt(1)),mu1x(i,:)',reshape(P1x(i,:),nx,nx)+model.Q,nx);
           end
       end
   end
   f=N*wx;
   Aeq=ones(1,ng);
   beq=1;
   lb=zeros(1,ng);
   ub=ones(1,ng);
    w1x= quadprog(M,-f,[],[],Aeq,beq,lb,ub);
    GMM.mu=mu1x;
    GMM.P=P1x;
    GMM.w=w1x;
else

    GMM.mu=mu1x;
    GMM.P=P1x;
    GMM.w=wx;
end



%% Measurement Step- Kalman Update of each component
if model.frequpdt(timeStep)==0
    
    GMM.mean=0;
    GMM.cov=0;
    for i=1:1:ng
        GMM.mean=GMM.mean+GMM.w(i)*GMM.mu(i,:)';
        GMM.cov=GMM.cov+GMM.w(i)*(reshape(GMM.P(i,:),nx,nx)+GMM.mu(i,:)'*GMM.mu(i,:));
    end
    GMM.cov=GMM.cov-GMM.mean*GMM.mean';
return
end

%% generate the measurement with this configuration
x_mc=sens.x_mc;
ym=model.hx(x_mc(timeStep,:)',[sens.phi,sens.alpha],sens.xc1)+model.gx(x_mc(timeStep,:)',[sens.phi,sens.alpha],sens.xc1)*mvnrnd(zeros(1,model.hn),model.R, 1)';
sens.ym=ym;
sens.gg=model.gx(x_mc(timeStep,:)',[sens.phi,sens.alpha],sens.xc1);
sens.hh=model.hx(x_mc(timeStep,:)',[sens.phi,sens.alpha],sens.xc1);

if (sens.gg(1,1)==100 && sens.gg(2,2)==100)

return
end
%% Kalmana like update
mu2x=zeros(size(mux));
P2x=zeros(size(Px));
muz=zeros(ng,nz);
Pz=zeros(ng,nz^2);
Pxz=zeros(ng,nx*nz);
for i=1:1:ng
    mu=GMM.mu(i,:)';
    sqP=sqrtm(reshape(GMM.P(i,:),nx,nx));
    
    for j=1:1:length(W)
       Xj=sqP*X(j,:)'+mu;
       hx=model.hx(Xj,[sens.phi,sens.alpha],sens.xc1);
       GG=model.gx(Xj,[sens.phi,sens.alpha],sens.xc1)*model.R*model.gx(Xj,[sens.phi,sens.alpha],sens.xc1)';

       muz(i,:)=muz(i,:)+W(j)*hx';
       Pz(i,:)=Pz(i,:)+W(j)*reshape(GG+hx*hx',1,nz^2);
       
       Pxz(i,:)=Pxz(i,:)+W(j)*reshape(Xj*hx',1,nx*nz);
    end
    Pzi=reshape(Pz(i,:),nz,nz)-muz(i,:)'*muz(i,:);
    Pxzi=reshape(Pxz(i,:),nx,nz)-mu*muz(i,:);
    Pz(i,:)=reshape(Pzi,1,nz^2);
    Pxz(i,:)=reshape(Pxzi,1,nx*nz);
    
    Ki=Pxzi/Pzi;
    mu2x(i,:)=GMM.mu(i,:)+(Ki*(ym-muz(i,:)'))';
    P2x(i,:)=reshape(reshape(GMM.P(i,:),nx,nx)-Ki*Pzi*Ki',1,nx^2);
end
    GMM.mu=mu2x;
    GMM.P=P2x;

    GMM.mean=0;
    GMM.cov=0;
 
    for i=1:1:ng
        GMM.mean=GMM.mean+GMM.w(i)*mu2x(i,:)';
        GMM.cov=GMM.cov+GMM.w(i)*(reshape(P2x(i,:),nx,nx)+mu2x(i,:)'*mu2x(i,:));
    end
    GMM.cov=GMM.cov-GMM.mean*GMM.mean';

