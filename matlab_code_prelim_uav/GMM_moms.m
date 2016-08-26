function GMM=GMM_moms(GMM)
ng=length(GMM.w);
nx=size(GMM.mu,2);
%% 
GMM.mean=zeros(nx,1);
GMM.cov=zeros(nx,nx);
for i=1:1:ng
    GMM.mean=GMM.mean+GMM.w(i)*GMM.mu(i,:)';
    GMM.cov=GMM.cov+GMM.w(i)*(reshape(GMM.P(i,:),nx,nx)+GMM.mu(i,:)'*GMM.mu(i,:));
end
GMM.cov=GMM.cov-GMM.mean*GMM.mean';
end