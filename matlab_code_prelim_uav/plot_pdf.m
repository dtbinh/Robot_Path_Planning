function plot_pdf(GMM,xc,yc,model)


ng=length(GMM.w); % no. of gaussian components
nx=model.fn; %size of target dynamics
nz=model.hn; %size of sensor model
nw=length(model.R);
marg_nx=model.marg_nx;
moden=length(model.marg_nx);
%% marginalization of the summ of gaussians
if model.marg_nx~=0
GMM.mu=GMM.mu(:,marg_nx);
Px=zeros(ng,length(marg_nx)^2);
for i=1:1:ng
    P=reshape(GMM.P(i,:),nx,nx);
    Px(i,:)=reshape(P(marg_nx,marg_nx),1,length(marg_nx)^2);
end
nx=length(model.marg_nx);
end
GMM.P=Px;


%% plotting
sigma=[];
for i=1:1:length(ng)
sigma = cat(3,sigma,reshape(GMM.P(i,:),moden,moden));
end
 obj = gmdistribution(GMM.mu,sigma,GMM.w');
 
 
[xx1,yy1]=meshgrid(GMM.mean(1,1)-40:1:GMM.mean(1,1)+40,GMM.mean(2,1)-40:1:GMM.mean(2,1)+40);

CC1=zeros(size(xx1));
for i=1:1:size(xx1,1)
    for j=1:1:size(xx1,2)       
  CC1(i,j)=pdf(obj,[xx1(i,j) yy1(i,j)]);
    end
end
contour(xx1,yy1,CC1,7,'linewidth',2)
% hold on
% contour(xx2,yy2,CC2,5)
 axis([-10,xc+10,-10,yc+10])
%  hold off
%    for i=1:1:size(xx1,1)
%     for j=1:1:size(xx1,2)
%         PP=reshape(GMM.P,model.fn,model.fn);
%         if min(eig(PP(1:2,1:2)))<0.1
%             PP(1:2,1:2)=0.1*eye(2);
%         end
%         if min(eig(PP(3:4,3:4)))<0.1
%             PP(3:4,3:4)=0.1*eye(2);
%         end 
%         CC1(i,j)=mvnpdf([xx1(i,j),yy1(i,j)],GMM.mu(1:2),PP(1:2,1:2));
%         CC2(i,j)=mvnpdf([xx2(i,j),yy2(i,j)],GMM.mu(3:4),PP(3:4,3:4));
%  CC(i,j)=pdf(obj,[xx(i,j) yy(i,j)]);
%     end
%    end
