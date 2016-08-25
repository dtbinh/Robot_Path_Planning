function GMM=post_GMM2_modf_meas(prior_GMM,z,h,g,R,SS)
GMM.w=prior_GMM.w;


ng=length(prior_GMM.w);
nz=length(z);
nx=size(prior_GMM.mu,2);

GMM.mu=zeros(size(prior_GMM.mu));
GMM.P=zeros(size(prior_GMM.P));
GMMz.mu=zeros(ng,nz);
GMMz.P=zeros(ng,nz^2);
Pxz=zeros(ng,nx*nz);

[X,W]=GH_pts(zeros(nx,1),eye(nx),2);
% keyboard
for i=1:1:ng
    mu=prior_GMM.mu(i,:)';
    sqP=sqrtm(reshape(prior_GMM.P(i,:),nx,nx));
    
    for j=1:1:length(W)
       Xj=sqP*X(j,:)'+mu;
%        keyboard
       [hx,GG]=comb_meas_props(Xj,h,g,R,SS);
%        hx=h(Xj,s);
%        GG=g(Xj,s)*R*g(Xj,s)';

       GMMz.mu(i,:)=GMMz.mu(i,:)+W(j)*hx';
%         keyboard
       GMMz.P(i,:)=GMMz.P(i,:)+W(j)*reshape(GG+hx*hx',1,nz^2);
       
       Pxz(i,:)=Pxz(i,:)+W(j)*reshape(Xj*hx',1,nx*nz);
    end
    GMMz.P(i,:)=reshape(reshape(GMMz.P(i,:),nz,nz)-GMMz.mu(i,:)'*GMMz.mu(i,:),1,nz^2);
    Pxzi=reshape(Pxz(i,:),nx,nz)-mu*GMMz.mu(i,:);
%     Pz(i,:)=reshape(Pzi,1,nz^2);
%     Pxz(i,:)=reshape(Pxzi,1,nx*nz);
    
    Ki=Pxzi/reshape(GMMz.P(i,:),nz,nz);
    GMM.mu(i,:)=prior_GMM.mu(i,:)+(Ki*(z-GMMz.mu(i,:)'))';
    GMM.P(i,:)=reshape(reshape(prior_GMM.P(i,:),nx,nx)-Ki*reshape(GMMz.P(i,:),nz,nz)*Ki',1,nx^2);
end
%     GMM.mu=mu2x;
%     GMM.P=P2x;

    

    GMM.mean=0;
    GMM.cov=0;
 
    for i=1:1:ng
        GMM.mean=GMM.mean+prior_GMM.w(i)*GMM.mu(i,:)';
        GMM.cov=GMM.cov+prior_GMM.w(i)*(reshape(GMM.P(i,:),nx,nx)+GMM.mu(i,:)'*GMM.mu(i,:));
    end
    GMM.cov=GMM.cov-GMM.mean*GMM.mean';