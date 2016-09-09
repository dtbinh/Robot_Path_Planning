%% an example showing that gaussian is not good enough
  y=[zeros(1,2);y1;y2];
  lam=lam_2xy(t,:);
[xx,yy]=meshgrid(-1:0.02:1);
 XX=zeros(size(xx));
 YY=zeros(size(xx));
 pdf=zeros(size(xx));
for i=1:1:size(xx,1)
         for j=1:1:size(xx,2)
             pdf(i,j)= pdf_MaxEnt([xx(i,j),yy(i,j)],lam,y);
             [XB,o]=transform_domain([xx(i,j),yy(i,j)],[-1,-1],[1,1],[10,10],[30,30]);
             XX(i,j)=XB(1);
             YY(i,j)=XB(2);
         end
end
 PDF=pdf./abs(prod(o));
     c=max(max(PDF));
     figure
     contour(XX,YY,PDF,linspace(0.1640/100,0.1640,20),'linewidth',2)
     figure
     mesh(XX,YY,PDF)
%      hold on
%      plot(DXY(:,1),DXY(:,2),'r.')
     
     