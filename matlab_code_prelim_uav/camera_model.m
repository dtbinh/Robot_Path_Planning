% camera model
function Ixy=camera_model(pos,fRg,Xt,cam_para)
% fRg: is the dir cosine matrix from ground to camera coordinate system
f=cam_para.focal;
dim=cam_para.dim;
% dim1 is along x and dim2 is along y
npix=cam_para.npix;

% number of measurments sets
nm=size(Xt,1)/3;
Ixy=0;
for i=1:1:nm
    x=Xt(3*i-2:1:3*i,1);
    lg=1/norm(pos-x)*(pos-x);
    lf=fRg*lg;
    xf=f/lf(3)*lf;
    if xf(1)>dim(1)/2 || xf(1)<-dim(1)/2
        xf(1)=NaN;
        xf(2)=NaN;
        xf(3)=NaN;
    end
    if xf(2)>dim(2)/2 || xf(2)<-dim(2)/2
        xf(1)=NaN;
        xf(2)=NaN;
        xf(3)=NaN;
    end
    
    Ixy=vertcat(Ixy,xf(1:2,1));
end
end