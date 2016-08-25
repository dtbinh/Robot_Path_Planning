function [TargetPos,xk1,Pk1]=TargetPriorProp(ntars,k,TargetPos,Time)
%from k to k+1
% can update all or some or 1
for targid=ntars
    Xfprev=TargetPos.xf{targid}(k,:)';
    Pfprev=reshape(TargetPos.Pf{targid}(k,:),TargetPos.fn(targid),TargetPos.fn(targid));
    
    
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
    
    [Xk,wk]=qd_pts(Xfprev,Pfprev);
    
    Xk1=zeros(size(Xk));
    for i=1:1:length(wk)
        if TargetPos.Dyn==1
            try
                Xk1(i,:)=KIRB_UM_eg_dyn_disc_manv(Xk(i,:)',Time.dt,TargetPos.Traffic,TargetPos.Trafficdx);
            catch
                keyboard
            end
        elseif TargetPos.Dyn==2
            Xk1(i,:)=KIRB_UM_eg_dyn_disc_manv_2(Xk(i,:)',Time.dt,TargetPos.Traffic,TargetPos.Trafficdx);
        end
    end
    
    [xk1,Pk1]=MeanCov(Xk1,wk);
    Pk1=Pk1+TargetPos.Q{targid};
    
    if norm(xk1(1:2)-Xfprev(1:2))>150
        keyboard
    end
    
    TargetPos.xf{targid}(k+1,:)=xk1;
    TargetPos.Pf{targid}(k+1,:)=reshape(Pk1,1,TargetPos.fn(targid)^2);
end


end
