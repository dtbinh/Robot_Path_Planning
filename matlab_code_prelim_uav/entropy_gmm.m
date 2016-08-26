function H=entropy_gmm(GMM)

ng=length(GMM.w);
nx=size(GMM.mu,2);

[Z,W] = GH_points(zeros(nx,1),eye(nx),3);

Hi=zeros(ng,1);
for i=1:1:ng
     Pi=reshape(GMM.P(i,:),nx,nx);
     sqPi=sqrtm(Pi);
     
     for j=1:1:length(W)
        Xj=sqPi*Z(j,:)'+GMM.mu(i,:)';
        S=0;  
        for k=1:1:ng
            S=S+GMM.w(k)*mvnpdf(Xj',GMM.mu(k,:),reshape(GMM.P(k,:),nx,nx));
        end
       Hi(i)=Hi(i)+W(j)*log(S);
        
     end
end
H=-GMM.w'*Hi;