function Tsk=SimpleSensortask_move(ii,jj,k,xf,Pf,SensorModel)
% just pair it with all the sensors that are within the FOV

    switch lower(SensorModel.filtermethod)
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
    
[Xk,wk]=qd_pts(xf(:),Pf);
% keyboard
Tsk=[];
if strcmp(SensorModel.Move.Mode,'off')
    return;
end
for n=1:1:SensorModel.Move.Ns
    G=0;
    for i=1:1:length(wk)
        if isempty(ii) || isempty(jj)
            pos=SensorModel.Move.Grid2Abs(SensorModel.Move.Gridpos{n}(k,1),SensorModel.Move.Gridpos{n}(k,2),SensorModel);
        else
            pos=SensorModel.Move.Grid2Abs(ii,jj,SensorModel);
        end
        try
            [y,g]=SensorModel.Move.hx{n}(Xk(i,:),pos,SensorModel.Move.FOV{n});
        catch
            keyboard
        end
        if isnan(g)==1
            G=G+1;
        end
    end
    if strcmp(SensorModel.Move.Mode,'passive')
        Tsk=vertcat(Tsk,n);
        continue;
    end
    if G/length(wk)<=0.75
        Tsk=vertcat(Tsk,n);
    end
end


end