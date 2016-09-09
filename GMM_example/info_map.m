%% calculating information potential
% load('prior_mixture')
close all

Ntar=5;
prior_GMM.mu=cell(Ntar,1);
prior_GMM.P=cell(Ntar,1);
prior_GMM.w=cell(Ntar,1);

prior_GMM.mu{1}=[300,300];
prior_GMM.P{1}=[200,0,0,200];
prior_GMM.w{1}=[1];

prior_GMM.mu{2}=[120,120];
prior_GMM.P{2}=[100,0,0,100];
prior_GMM.w{2}=[1];

prior_GMM.mu{3}=[100,250];
prior_GMM.P{3}=[100,0,0,10];
prior_GMM.w{3}=[1];

prior_GMM.mu{4}=[200,200];
prior_GMM.P{4}=[100,0,0,10];
prior_GMM.w{4}=[1];

prior_GMM.mu{5}=[270,100];
prior_GMM.P{5}=[10,5,5,200];
prior_GMM.w{5}=[1];


h=cell(Ntar,1);
g=cell(Ntar,1);
% h{1}=@(x,s)[(x(1)-s(1));(x(2)-s(2))];
% h{2}=@(x,s)[(x(1)-s(1));(x(2)-s(2))];

h{1}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];
h{2}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];
h{3}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];
h{4}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];
h{5}=@(x,s)[sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)];

% h{1}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));
% h{2}=@(x,s)atan2((x(2)-s(2)),(x(1)-s(1)));

r=60; % is the field of view 

% variance penalty model
R=cell(Ntar,1);
R{1}=diag([1].^2);%+(r/4)^2;
R{2}=diag([5].^2);%+(r/4)^2;
R{3}=diag([5].^2);
R{4}=diag([5].^2);
R{5}=diag([5].^2);

g{1}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);
g{2}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);
g{3}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);
g{4}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);
g{5}=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,1/(2*r^4)*sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)^8+1/2);

% g=@(x,s)conditional(sqrt((x(1)-s(1))^2+(x(2)-s(2))^2)<=r,1,NaN);
nz=cell(2,1);
nz{1}=1;
nz{2}=1;
nz{3}=1;
nz{4}=1;
nz{5}=1;

[s1,s2]=meshgrid(0:5:400);
ExpI=zeros(size(s2));
tic
% k=1;
for i=1:1:size(s1,1)
    parfor j=1:1:size(s1,2)
%         k
%         k=k+1;
        m=zeros(Ntar,1);
        
        for nt=1:1:Ntar
prGMM=struct('mu',0,'P',0,'w',0);
        prGMM.mu=prior_GMM.mu{nt};
        prGMM.P=prior_GMM.P{nt};
        prGMM.w=prior_GMM.w{nt};
        
        measGMM=meas_GMM(prGMM,nz{nt},h{nt},g{nt},R{nt},[s1(i,j),s2(i,j)]);
        m(nt)=Expt_GMM(measGMM,@(z)dist_eval_meas(prGMM,z,h{nt},g{nt},R{nt},[s1(i,j),s2(i,j)],1));
        end
        ExpI(i,j)=sum(m);
        
    end
end
toc
% figure(1)
% mesh(s1,s2,ExpI)
% xlabel('s1')
% ylabel('s2')
% alpha(0.5)
% hold on
% plot(XYt(:,1),XYt(:,2),'ro')
% contour(s1,s2,ExpI,0:0.005:0.3)
%% polynomial surface fitting
% phi=@(x)[1,x(1),x(2),x(1)^2,x(2)^2,x(1)*x(2),x(1)^3,x(2)^3,x(1)^2*x(2),x(1)*x(2)^2,x(1)^4,x(2)^4,x(1)^3*x(2),x(1)^2*x(2)^2,x(1)^1*x(2)^3,x(1)^5,x(2)^5,x(1)^4*x(2),x(1)^3*x(2)^2,x(1)^2*x(2)^3,x(1)^1*x(2)^4,x(1)^2,x(2)^6,x(1)^5*x(2)^1,x(1)^4*x(2)^2,x(1)^3*x(2)^3,x(1)^2*x(2)^4,x(1)^1*x(2)^5]';
% A=zeros(1,length(phi([1,2])));
% B=[0;0];
% k=1;
% for i=1:1:size(s1,1)
%     for j=1:1:size(s1,2)
%      B(k)=ExpI(i,j);
%      A(k,:)=phi([s1(i,j),s2(i,j)])';
%         k=k+1;
%     end
% end
%   C= A\B; 
  [xx,yy]=meshgrid(0:5:400);
  apprxSurf=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
%      apprxSurf(i,j)=C'*phi([s1(i,j),s2(i,j)]);
     apprxSurf(i,j)= interp2(s1,s2,ExpI,xx(i,j),yy(i,j),'cubic');
    end
end
figure
% mesh(s1,s2,ExpI)
xlabel('s1')
ylabel('s2')
hold on
   mesh(xx,yy,-(apprxSurf))
%   figure(2)
%    mesh(s1,s2,ExpI)

   
   %% max ent pdf fitting 
% k=1;
% f=[0;0];
% pts=[0,0];
% for i=1:1:size(s1,1)
%     for j=1:1:size(s1,2)
%      f(k)=ExpI(i,j);
%      pts(k,:)=[s1(i,j),s2(i,j)]/100;
%         k=k+1;
%     end
% end
% phi=@(x)[1,x(1),x(2),x(1)^2,x(2)^2,x(1)*x(2),x(1)^3,x(2)^3,x(1)^2*x(2),x(1)*x(2)^2]';
% % ,x(1)^4,x(2)^4,x(1)^3*x(2),x(1)^2*x(2)^2,x(1)^1*x(2)^3,x(1)^5,x(2)^5,x(1)^4*x(2),x(1)^3*x(2)^2,x(1)^2*x(2)^3,x(1)^1*x(2)^4,x(1)^6,x(2)^6,x(1)^5*x(2)^1,x(1)^4*x(2)^2,x(1)^3*x(2)^3,x(1)^2*x(2)^4,x(1)^1*x(2)^5
% options=optimset('Display','Iter','MaxFunEvals',10000,'MaxIter',2000);
% lam=lsqcurvefit(@(lam,s)max_ent_surf_fit(lam,s,phi),zeros(length(phi([1,2])),1),pts,f,[],[],options);
% MaxentapprxSurf=zeros(size(xx));
% for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
% %      apprxSurf(i,j)=C'*phi([s1(i,j),s2(i,j)]);
%      MaxentapprxSurf(i,j)= max_ent_surf_fit(lam,[xx(i,j)/100,yy(i,j)/100],phi);
%     end
% end
% mesh(xx,yy,MaxentapprxSurf)