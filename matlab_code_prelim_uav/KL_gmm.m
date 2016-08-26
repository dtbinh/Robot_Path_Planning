function H=KL_gmm(GMM1,GMM2)

ng1=length(GMM1.w);
ng2=length(GMM2.w);

nx=size(GMM1.mu,2);

[Z,W] = GH_points(zeros(nx,1),eye(nx),3);

Hi=zeros(ng1,1);
for i=1:1:ng1
     Pi=reshape(GMM1.P(i,:),nx,nx);
     sqPi=sqrtm(Pi);
     
     for j=1:1:length(W)
        Xj=sqPi*Z(j,:)'+GMM1.mu(i,:)';
        S1=0;
        S2=0;
        for k=1:1:ng1
            S1=S1+GMM1.w(k)*mvnpdf(Xj',GMM1.mu(k,:),reshape(GMM1.P(k,:),nx,nx));
        end
        for k=1:1:ng2
            S2=S2+GMM2.w(k)*mvnpdf(Xj',GMM2.mu(k,:),reshape(GMM2.P(k,:),nx,nx));
        end
       Hi(i)=Hi(i)+W(j)*log(S1/S2);
        
     end
end
H=GMM1.w'*Hi;