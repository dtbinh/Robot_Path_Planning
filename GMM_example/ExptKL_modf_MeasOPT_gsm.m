function DKL=ExptKL_modf_MeasOPT_gsm(GMM,model,para)
% mux is the vector of means
%Px is the cov as a row vector for each component
%wx is the weights of the gauss components
mux=GMM.mu;
Px=GMM.P;
wx=GMM.w;


%%
ng=size(mux,1); % no. of gaussian components
nx=model.fn; %size of target dynamics
nz=model.hn; %size of sensor model
nw=length(model.R);
marg_nx=model.marg_nx;
%% marginalization of the summ of gaussians
if model.marg_nx~=0
mux=GMM.mu(:,marg_nx);
Px=zeros(ng,length(marg_nx)^2);
for i=1:1:ng
    P=reshape(GMM.P(i,:),nx,nx);
    Px(i,:)=reshape(P(marg_nx,marg_nx),1,length(marg_nx)^2);
end
nx=length(model.marg_nx);
end
    
%%


if model.InfoQuadN==1
Nquad=3;
[X,W]=GH_pts(zeros(nx,1),eye(nx),Nquad);
[Z,Wz]=GH_pts(zeros(nz,1),eye(nz),Nquad);
else
    switch lower(model.InfoMethod)
      case 'ut'
          r1=sqrt(nx+model.Infokappa);
          X=[zeros(1,nx);r1*eye(nx);-r1*eye(nx)];
          W=[model.Infokappa/(nx+model.Infokappa);1/(2*(nx+model.Infokappa))*ones(2*nx,1)];
          
          r1=sqrt(nz+model.Infokappa);
          Z=[zeros(1,nz);r1*eye(nz);-r1*eye(nz)];
          Wz=[model.Infokappa/(nz+model.Infokappa);1/(2*(nz+model.Infokappa))*ones(2*nz,1)];
          
      case 'ckf'    
          X=sqrt(nx)*[eye(nx);-eye(nx)];
          W=1/(2*nx)*ones(2*nx,1);
          Z=sqrt(nz)*[eye(nz);-eye(nz)];
          Wz=1/(2*nz)*ones(2*nz,1);
      case 'gh'
          [X,W]=GH_pts(zeros(nx,1),eye(nx),model.Infoghpts);
          [Z,Wz]=GH_pts(zeros(nz,1),eye(nz),model.Infoghpts);
      case 'cut4'
          [X,W]=conjugate_dir_gausspts(zeros(nx,1),eye(nx));
          [Z,Wz]=conjugate_dir_gausspts(zeros(nz,1),eye(nz));
      case 'cut6'
          [X,W]=conjugate_dir_gausspts_till_6moment_scheme2(zeros(nx,1),eye(nx));
          [Z,Wz]=conjugate_dir_gausspts_till_6moment_scheme2(zeros(nz,1),eye(nz));
      case 'cut8'
          [X,W]=conjugate_dir_gausspts_till_8moment(zeros(nx,1),eye(nx));
          [Z,Wz]= conjugate_dir_gausspts_till_8moment(zeros(nz,1),eye(nz));
    end 
end

th=para;
h=@(x,th)model.hx(x,th);
g=@(x,th)model.gx(x,th);
R=model.R;

A1=zeros(ng,1);
muzi=zeros(ng,nz);
Mzi=zeros(ng,nz^2);
Pzi=zeros(ng,nz^2);
ZZ=zeros(ng,1);



for i=1:1:ng
     muxi=mux(i,:)';
     Pxi=reshape(Px(i,:),nx,nx);
     sqPxi=sqrtm(Pxi);
    for j=1:1:length(W)
        Xj=sqPxi*X(j,:)'+muxi;
        GG=g(Xj,th)*R*g(Xj,th)';
  
        hj=h(Xj,th);
        
        A1(i)=A1(i)+W(j)*log((2*pi*exp(1))^nz*det(GG));
        muzi(i,:)=muzi(i,:)+W(j)*hj';
        Mzi(i,:)=Mzi(i,:)+W(j)*reshape(GG+hj*hj',1,nz^2);

    
 
        
    end
    Pzi(i,:)=reshape(reshape(Mzi(i,:),nz,nz)-muzi(i,:)'*muzi(i,:),1,nz^2);
end
A1=wx'*A1;



%now calculating the entropy of z

Hi=zeros(ng,1);
for i=1:1:ng
     PPzi=reshape(Pzi(i,:),nz,nz);
     sqPzi=sqrtm(PPzi);
     
     for j=1:1:length(Wz)
        Zj=sqPzi*Z(j,:)'+muzi(i,:)';
        S=0;  
        for k=1:1:ng
            S=S+wx(k)*mvnpdf(Zj',muzi(k,:),reshape(Pzi(k,:),nz,nz));
        end
       Hi(i)=Hi(i)+Wz(j)*log(S);
        
     end
end
Hz=wx'*Hi;

DKL=-0.5*A1-Hz;



