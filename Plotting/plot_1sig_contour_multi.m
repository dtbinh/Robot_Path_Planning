function plot_1sig_contour_multi(TargetPos,kk)
for ntar=1:1:TargetPos.Ntar
 Xf=TargetPos.xf{ntar}(kk,:)';
 Pf=reshape(TargetPos.Pf{ntar}(kk,:),TargetPos.fn(ntar),TargetPos.fn(ntar));
 
     switch lower(TargetPos.filtermethod)
        case 'ckf'
            qd_pts=@cubature_KF_points;
        case 'ut'
            qd_pts=@(m,P)UT_sigmapoints(m,P,2);
        case 'cut4'
            qd_pts=@conjugate_dir_gausspts;
        case 'cut6'
            qd_pts=@conjugate_dir_gausspts_till_6moment_scheme2;
        case 'cut8'
            qd_pts=@conjugate_dir_gausspts_till_8moment;
        case 'gh'
            qd_pts=@(m,P)GH_pts(m,P,para);
        otherwise
            error('smthg is wrong: DONT ask me what')
     end
    
 [Xpt,w]=qd_pts(Xf,Pf);
Y=zeros(length(0:0.1:360),2);
muv=Xf(3:4);
muv=muv/norm(muv);

mux=Xf(1:2);
mux=mux(:);
Px=Pf(1:2,1:2);

A=sqrtm(Px);
k=1;
for th=0:0.1:360
    X=[cosd(th);sind(th)];
    Y(k,:)=A*X+mux;
    k=k+1;
end
% muv=A*[1;0]+mux;
plot(TargetPos.xf{ntar}(1:kk,1),TargetPos.xf{ntar}(1:kk,2),'r--',TargetPos.truth{ntar}(1:kk,1),TargetPos.truth{ntar}(1:kk,2),'k*-',Xpt(:,1),Xpt(:,2),'go',mux(1),mux(2),'ks',Y(:,1),Y(:,2),'g','linewidth',2)
% Xpt [mux(1),muv(1)+mux(1)],[mux(2),muv(2)+mux(2)]
end



end