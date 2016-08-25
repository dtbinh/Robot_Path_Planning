function plot_UAV_GMM(c0,d,xlim,prior_GMM,sts,r,Ntar)
[Ygrid,Xgrid]=meshgrid(0:d:xlim);

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
end
for ss=1:1:size(c0,1)
            if c0{ss}(3)==1
            plot(Xgrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2)),'ko',[Xgrid(c0{ss}(1),c0{ss}(2)),Xgrid(c0{ss}(1),c0{ss}(2))+10],[Ygrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2))],'b','linewidth',2);
            plot_triangle(15,Xgrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2)),c0{ss}(3))
            elseif c0{ss}(3)==2
            plot(Xgrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2)),'ko',[Xgrid(c0{ss}(1),c0{ss}(2)),Xgrid(c0{ss}(1),c0{ss}(2))],[Ygrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2))+10],'b','linewidth',2);
            plot_triangle(15,Xgrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2)),c0{ss}(3))
            elseif c0{ss}(3)==3
            plot(Xgrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2)),'ko',[Xgrid(c0{ss}(1),c0{ss}(2)),Xgrid(c0{ss}(1),c0{ss}(2))-10],[Ygrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2))],'b','linewidth',2);
            plot_triangle(15,Xgrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2)),c0{ss}(3))
            elseif c0{ss}(3)==4
            plot(Xgrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2)),'ko',[Xgrid(c0{ss}(1),c0{ss}(2)),Xgrid(c0{ss}(1),c0{ss}(2))],[Ygrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2))-10],'b','linewidth',2);
            plot_triangle(15,Xgrid(c0{ss}(1),c0{ss}(2)),Ygrid(c0{ss}(1),c0{ss}(2)),c0{ss}(3))
            end
end
%  plot cirlc FOV
for ss=1:1:size(c0,1)
    xx=0;
    yy=0;
    k=1;
    for th=0:0.01:2*pi
    xx(k)=Xgrid(c0{ss}(1),c0{ss}(2))+r*cos(th);
    yy(k)=Ygrid(c0{ss}(1),c0{ss}(2))+r*sin(th);
    k=k+1;
    end
    plot(xx,yy,'k','linewidth',2)
    fill(xx,yy,'r')
    alpha 0.2
end
            axis([-10,xlim+10,-10,xlim+10])
            set(gca,'XTick',0:d:xlim)
            set(gca,'YTick',0:d:xlim)
            grid on
%% now plotiing the path history for each sensor
% for ss=1:1:size(c0,1)
%     X=[];
%     cc0=[s0{ss}(1,1),s0{ss}(1,2),s0{ss}(1,3)];
%     for t=1:1:size(Hist.u)
%         XX=temp{Hist.u{t}(1,ss),1};
%         
%         X=vertcat(X,XX);        
%     end
%     
% end

hold off

end










function plot_triangle(l,c1,c2,dir)

x1=0:0.01:l*cosd(60);
y1=tand(60).*x1;
x2=l*cosd(60):0.01:l;
y2=tand(120).*(x2-l);
x2=x2(end:-1:1);
y2=y2(end:-1:1);
x3=0:0.01:l;
y3=0.*x3;
x3=x3(end:-1:1);
y3=y3(end:-1:1);

X=[x1',y1';x2',y2';x3',y3'];
X(:,1)=X(:,1)-l*cosd(60);
X(:,2)=X(:,2)-l*cosd(60)*tand(30);
Y=zeros(size(X));
if dir==1
  Y(:,1)=X(:,2);
  Y(:,2)=-X(:,1);
end
if dir==2
  Y=X;
end
if dir==3
  Y(:,1)=-X(:,2);
  Y(:,2)=X(:,1);
end
if dir==4
  Y(:,1)=-X(:,1);
  Y(:,2)=-X(:,2);
end
plot(Y(:,1)+c1,Y(:,2)+c2,'k','linewidth',2)
end