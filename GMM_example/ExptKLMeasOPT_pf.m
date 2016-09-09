function DKL=ExptKLMeasOPT_pf(Xx,wx,model,sens_conf)
ns=length(wx);
ng=3;
nx=model.fn; %size of target dynamics
nz=model.hn; %size of sensor model
nw=length(model.R);

options = statset('MaxIter',5000);

%% generating the samples for z
%and also evaluating the integrals
th=sens_conf;
h=@(x,th)model.hx(x,th);
g=@(x,th)model.gx(x,th);
R=model.R;
Xz=zeros(ns,nz);
A1=0;
A2=zeros(nz,nz);
A3=zeros(nz,1);
A4=0;
for i=1:1:ns
   GG=g(Xx(i,:)',th)*R*g(Xx(i,:)',th)';
   inGG=inv(GG);
   hj=h(Xx(i,:)',th);
    
     A1=A1+wx(i)*log(det(GG));
     A2=A2+wx(i)*inGG;
     A3=A3+wx(i)*(inGG*hj);
     A4=A4+wx(i)*hj'*inGG*hj;
     
    Xz(i,:)=(h(Xx(i,:)',th)+g(Xx(i,:)',th)*mvnrnd(zeros(1,nw),R, 1)')';
     
end

    W=repmat(wx,1,nz);
    muz=sum(W.*Xz,1)';
    
MU=repmat(muz',ns,1);
X=Xz-MU;
Pz=X'*(W.*X);

%% Fitting a gaussian distribution to Xz
obj = gmdistribution.fit(Xz,ng,'SharedCov',true,'Options',options);
% plot(Xz(:,1),Xz(:,2),'bo')


Hz=wx'*log(pdf(obj,Xz));
Hzx=-nz/2*log(2*pi)-0.5*A1-0.5*trace(A2*Pz)-0.5*muz'*A2*muz+muz'*A3-0.5*A4;
DKL=Hz-Hzx;
end


