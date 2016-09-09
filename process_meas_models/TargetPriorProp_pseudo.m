function [TargetPos,XK1,PK1]=TargetPriorProp_pseudo(ntars,curk,futk,TargetPos,Time,further)
%from curk to futk, assuming curk is the initial condition
% can update all or some or 1
% M=GetCellOptions(further,opt)
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
Tn=curk:1:futk;
Tvec=Time.Tvec(Tn);
Tvec=Tvec(:)';
XK1=cell(length(1:TargetPos.Ntar),1);
PK1=cell(length(1:TargetPos.Ntar),1);

for targid=ntars
    XK1{targid}=zeros(Time.nt,TargetPos.fn(targid));
    PK1{targid}=zeros(Time.nt,TargetPos.fn(targid)^2);
    
    
    for k=1:1:length(Tn)-1
        Tk=Tn(k);
        Tk1=Tn(k+1);
        Xfprev=TargetPos.xf{targid}(Tk,:)';
        try
        Pfprev=reshape(TargetPos.Pf{targid}(Tk,:),TargetPos.fn(targid),TargetPos.fn(targid));
        catch
            keyboard
        end
        
        if k==1
        XK1{targid}(Tk,:)=Xfprev;
        PK1{targid}(Tk,:)=reshape(Pfprev,1,TargetPos.fn(targid)^2);
        end
        
        [Xk,wk]=qd_pts(Xfprev,Pfprev);
        
        Xk1=zeros(size(Xk));
        for i=1:1:length(wk)
            Xk1(i,:)=TargetPos.dynModel{targid}(Xk(i,:)');
        end
        
        [xk1,Pk1]=MeanCov(Xk1,wk);
        
        [st,~]=GetCellOptions(further,'NoProcessNoise');
        if st==true
            
        else
            Pk1=Pk1+TargetPos.Q{targid};
        end
        
        if norm(xk1(1:2)-Xfprev(1:2))>150
            keyboard
        end
        
        [st,~]=GetCellOptions(further,'NoLoad2targModel');
        if st==true
            
        else
        TargetPos.xf{targid}(Tk1,:)=xk1;
        TargetPos.Pf{targid}(Tk1,:)=reshape(Pk1,1,TargetPos.fn(targid)^2);
        end
        
        XK1{targid}(Tk1,:)=xk1;
        PK1{targid}(Tk1,:)=reshape(Pk1,1,TargetPos.fn(targid)^2);
    end
    
    
    
end

end
