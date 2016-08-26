function DKL=ExptKLMeasOPT_gsm(GMM,model,sens_conf)
% mux is the vector of means
%Px is the cov as a row vector for each component
%wx is the weights of the gauss components
mux=GMM.mu;
Px=GMM.P;
wx=GMM.w;
ng=size(mux,1); % no. of gaussian components
nx=model.fn; %size of target dynamics
nz=model.hn; %size of sensor model
nw=length(model.R);

Nquad=4;

th=sens_conf;
h=@(x,th)model.hx(x,th);
g=@(x,th)model.gx(x,th);
R=model.R;

A1=zeros(ng,1);
A2=zeros(ng,nz^2);
A3=zeros(ng,nz);
A4=zeros(ng,1);
muzi=zeros(ng,nz);
Mzi=zeros(ng,nz^2);
Pzi=zeros(ng,nz^2);

[X,W]=GH_pts(zeros(nx,1),eye(nx),Nquad);

for i=1:1:ng
     muxi=mux(i,:)';
     Pxi=reshape(Px(i,:),nx,nx);
     sqPxi=sqrtm(Pxi);
    for j=1:1:length(W)
        Xj=sqPxi*X(j,:)'+muxi;
        GG=g(Xj,th)*R*g(Xj,th)';
        inGG=inv(GG);
        hj=h(Xj,th);
        
        A1(i)=A1(i)+W(j)*log(det(GG));
        A2(i,:)=A2(i,:)+W(j)*reshape(inGG,1,nz^2);
        A3(i,:)=A3(i,:)+W(j)*(inGG*hj)';
        A4(i)=A4(i)+W(j)*hj'*inGG*hj;
        muzi(i,:)=muzi(i,:)+W(j)*hj';
        Mzi(i,:)=Mzi(i,:)+W(j)*reshape(GG+hj*hj',1,nz^2);
    end
    Pzi(i,:)=reshape(reshape(Mzi(i,:),nz,nz)-muzi(i,:)'*muzi(i,:),1,nz^2);
end
A1=wx'*A1;
A2=reshape(sum(repmat(wx,1,nz^2).*A2,1),nz,nz);
A3=sum(repmat(wx,1,nz).*A3,1)';
A4=wx'*A4;
muz=sum(repmat(wx,1,nz).*muzi,1)';
Mz=reshape(sum(repmat(wx,1,nz^2).*Mzi,1),nz,nz);
Pz=Mz-muz*muz';

%now calculating the entropy of z
[Z,W]=GH_pts(zeros(nz,1),eye(nz),Nquad);
Hi=zeros(ng,1);
for i=1:1:ng
     PPzi=reshape(Pzi(i,:),nz,nz);
     sqPzi=sqrtm(PPzi);
     
     for j=1:1:length(W)
        Zj=sqPzi*Z(j,:)'+muzi(i,:)';
        S=0;  
        for k=1:1:ng
            S=S+wx(k)*mvnpdf(Zj',muzi(k,:),reshape(Pzi(k,:),nz,nz));
        end
       Hi(i)=Hi(i)+W(j)*log(S);
        
     end
end
Hz=wx'*Hi;

Hzx=-nz/2*log(2*pi)-0.5*A1-0.5*trace(A2*Pz)-0.5*muz'*A2*muz+muz'*A3-0.5*A4;
DKL=Hz-Hzx;




