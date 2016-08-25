function plot_UAV_GMM2(d,xlim,prior_GMM,sts,r,Ntar)
[Ygrid,Xgrid]=meshgrid(0:d:xlim);
%plot only the target positions and pdf contours

for nt=1:1:Ntar
plot(prior_GMM.truth{nt}(:,1),prior_GMM.truth{nt}(:,2),'ks-','MarkerSize',6,'linewidth',2)
hold on
end

for nt=1:1:Ntar
        prGMM.mu=prior_GMM.mu{nt};
        prGMM.P=prior_GMM.P{nt};
        prGMM.w=prior_GMM.w{nt};
        if sts~=0
%             keyboard
        prGMM=marg_GMM(prGMM,sts);
        end
        [mu,P]=GMM_mean_cov(prGMM);
        lam=max(eig(P));
        
        if size(prGMM.mu,1)==1
        [Xg,Yg]=meshgrid(linspace(mu(1)-2*lam,mu(1)+2*lam,100),linspace(mu(2)-2*lam,mu(2)+2*lam,100));
        obj=GMM2obj(prGMM);
        p=zeros(size(Xg));
        for j=1:1:size(Xg,1)
            for k=1:1:size(Xg,2)
           p(j,k)=pdf(obj,[Xg(j,k) Yg(j,k)]);
            end
        end
        contour(Xg,Yg,p,12,'linewidth',1)
        alpha 0.4
        else
         [Xg,Yg]=meshgrid(linspace(mu(1)-2*lam,mu(1)+2*lam,200),linspace(mu(2)-2*lam,mu(2)+2*lam,200));
        obj=GMM2obj(prGMM);
        p=zeros(size(Xg));
        for j=1:1:size(Xg,1)
            for k=1:1:size(Xg,2)
           p(j,k)=pdf(obj,[Xg(j,k) Yg(j,k)]);
            end
        end
        contour(Xg,Yg,p,12,'linewidth',1)
        alpha 0.4
        end
end

hold off

end







