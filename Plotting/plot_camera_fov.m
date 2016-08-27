function plot_camera_fov(pos,fRg,cam_para,col)
% fRg: is the dir cosine matrix from ground to camera coordinate system
f=cam_para.focal;
dim=cam_para.dim;
% dim1 is along x and dim2 is along y
npix=cam_para.npix;

x=[dim(1),dim(2);-dim(1),dim(2);-dim(1),-dim(2);dim(1),-dim(2);dim(1),dim(2)];
X=zeros(1,5);
Y=zeros(1,5);

for i=1:1:5
    lf=[x(i,:)';f]/norm([x(i,:)';f]);
    lg=fRg\lf;
    c=-pos(3)/lg(3);
    rt=pos+c*lg;
    X(i)=rt(1);
    Y(i)=rt(2);
end
f=max(abs([X(1)-X(2),Y(2)-Y(3)]))/50;
xc=pos(1)-f:0.1:pos(1)+f;
ycp=pos(2)+sqrt(-(xc-pos(1)).^2+f^2);
ycm=pos(2)-sqrt(-(xc-pos(1)).^2+f^2);
hold on
plot(X,Y,'k','linewidth',2)

fill(X,Y,col,'FaceAlpha',0.1)
plot(xc,ycp,'k',xc,ycm,'k','linewidth',2)
hold off    
end
    