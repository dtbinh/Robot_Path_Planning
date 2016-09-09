function [TargetPos,XKF,PKF,XKU,PKU,flag]=TargetPostUpdate_pseudo(ntars,curk,type,SensorModel,TargetPos,Time,further)
% update the target at the same time curk
% update according to the tasks that have been assigned
% so if type=='true' have to generate a measurement and then update
%    if type=='pseudo' no measurement just update cov and prop mean

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


XKF=cell(length(1:TargetPos.Ntar),1);
PKF=cell(length(1:TargetPos.Ntar),1);

XKU=cell(length(1:TargetPos.Ntar),1);
PKU=cell(length(1:TargetPos.Ntar),1);

flag=ones(length(ntars),1);

for targid=ntars
    
    
    %     keyboard
    
    xfk=TargetPos.xf{targid}(curk,:)';
    Pfk=reshape(TargetPos.Pf{targid}(curk,:),TargetPos.fn(targid),TargetPos.fn(targid));
    
    XKF{targid}=xfk;
    PKF{targid}=Pfk;
    
    [Xk,wk]=qd_pts(xfk,Pfk);
    
    %Now assemble the R and ym according to the sensor tasking
    % Do this both for static and movable sensors
    
    %% No tasking considered, just do the measurement update with all the targets within the sensor FOV
    [st,~]=GetCellOptions(further,'DisregardTasking');
    if st==true
        RR=[];
        ym=[];
        yy=[];
        Y=[];
        ytp=[];
        
        [st,~]=GetCellOptions(further,'NoStaticSensors');
        if st==false
            for nsens=1:SensorModel.Static.Ns
                Gs=[];
                Gm=[];
                yy=zeros(length(wk),SensorModel.Static.hn(nsens));
                for i=1:1:length(wk)
                    [yy(i,:),g,tp]=SensorModel.Static.hx{nsens}(Xk(i,:),SensorModel.Static.pos{nsens},SensorModel.Static.FOV{nsens});
                    Gs=vertcat(Gs,g);
                end
                % true measurements
                [hh,g,tp]=SensorModel.Static.hx{nsens}(TargetPos.truth{targid}(curk,:),SensorModel.Static.pos{nsens},SensorModel.Static.FOV{nsens});
                Gm=vertcat(Gm,g);            
                yg=hh+sqrtm(SensorModel.Static.R{nsens})*randn(SensorModel.Static.hn(nsens),1);
                if length(yg)==2
                    x123=cos(yg(2));
                    y123=sin(yg(2));
                    yg(2)=atan2(y123,x123);
                end
                [ft,~]=GetCellOptions(further,'NoStaticSensorsFOV');
                if ft==true
                    ym=vertcat(ym,yg);
                    ytp=vertcat(ytp,tp);
                    Y=horzcat(Y,yy);
                    RR=blkdiag(RR,SensorModel.Static.R{nsens});
                else
                    if  sum(isnan(Gs))>0.75*length(Gs) || isnan(Gm)

                    else
                        ym=vertcat(ym,yg);
                        ytp=vertcat(ytp,tp);
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,SensorModel.Static.R{nsens});
                    end
                end
            end
        end
        
        [st,~]=GetCellOptions(further,'NoMoveSensros');
        if st==false
            for nsens=1:SensorModel.Move.Ns
                Gs=[];
                Gm=[];
                pos=SensorModel.Move.Grid2Abs(SensorModel.Move.Gridpos{nsens}(curk,1),SensorModel.Move.Gridpos{nsens}(curk,2),SensorModel);
                yy=zeros(length(wk),SensorModel.Move.hn(nsens));
                for i=1:1:length(wk)
                    [yy(i,:),gs,tp]=SensorModel.Move.hx{nsens}(Xk(i,:),pos,SensorModel.Move.FOV{nsens});
                    Gs=vertcat(Gs,gs);
                end
                % true measurements
                [hh,gm,tp]=SensorModel.Move.hx{nsens}(TargetPos.truth{targid}(curk,:),pos,SensorModel.Move.FOV{nsens});
                Gm=vertcat(Gm,gm);
                yg=hh+sqrtm(SensorModel.Move.R{nsens})*randn(SensorModel.Move.hn(nsens),1);
                if length(yg)==2
                    x123=cos(yg(2));
                    y123=sin(yg(2));
                    yg(2)=atan2(y123,x123);
                end
                if  sum(isnan(Gs))>0.75*length(Gs) || isnan(Gm)
                    
                else
                    ym=vertcat(ym,yg);
                    ytp=vertcat(ytp,tp);
                    Y=horzcat(Y,yy);
                    RR=blkdiag(RR,SensorModel.Move.R{nsens});
                end
            end
        end
        
        [st,~]=GetCellOptions(further,'pseudo');
        if st==true
            if isempty(Y)==0
                [mz,Pz]=MeanCov(Y,wk);
                Pz=Pz+RR;
                Pcc=CrossCov(Xk,xfk,Y,mz,wk);
                K=Pcc/Pz;
                
                xk=xfk;
                Pk=Pfk-K*Pz*K';
            else
                xk=xfk;
                Pk=Pfk;
                flag(targid)=NaN;
            end
        else
            if isempty(ym)==0
                [mz,Pz]=MeanCov(Y,wk);
                Pz=Pz+RR;
                Pcc=CrossCov(Xk,xfk,Y,mz,wk);
                K=Pcc/Pz;
                
                d=ym(:)-mz(:);
                ind=find(ytp==1);
                
                try
                    for sd=1:1:length(ind)
                        %                 qz=angle2quat(mz(ind(sd)),0,0);
                        %                 qm=angle2quat(ym(ind(sd)),0,0);
                        %                 d(sd)=quat2angle(quatmultiply(qm,quatinv(qz)));
                       
                        if ym(ind(sd))<0
                            ym(ind(sd))=ym(ind(sd))+2*pi;
                        end
                        if mz(ind(sd))<0
                            mz(ind(sd))=mz(ind(sd))+2*pi;
                        end
                        d(ind(sd))=ym(ind(sd))-mz(ind(sd));
                        x123=cos(d(ind(sd)));
                        y123=sin(d(ind(sd)));
                        d(ind(sd))=atan2(y123,x123);
                    end
                catch
                    keyboard
                end
                
                xk=xfk+K*(d);
                Pk=Pfk-K*Pz*K';
            else
                xk=xfk;
                Pk=Pfk;
                flag(targid)=NaN;
            end
        end
        
        if norm(xk(1:2)-xfk(1:2))>100 || min(sign(xk(1:2)))==-1 
            keyboard
        end
        
        TargetPos.xf{targid}(curk,:)=xk;
        TargetPos.Pf{targid}(curk,:)=reshape(Pk,1,TargetPos.fn(targid)^2);
        
        XKU{targid}=xk;
        PKU{targid}=Pk;
        
        continue
    end
    
    %% With tasking
    % SensorModel.Static.Task=cell(Time.nt,1);  % task [sensor,targ] then make the measurement
    % SensorModel.Move.Task=cell(Time.nt,1);
    Ns_static=SensorModel.Static.Task{curk}( find(SensorModel.Static.Task{curk}(:,2)==targid),1);% task [sensor,targ]
    Ns_move=SensorModel.Move.Task{curk}( find(SensorModel.Move.Task{curk}(:,2)==targid)  ,1);% task [sensor,targ]
    
    RR=[];
    ym=[];
    yy=[];
    Y=[];
    ytp=[];
%     keyboard
    [st,~]=GetCellOptions(further,'NoStaticSensors');
    if st==false
        for nsens=1:SensorModel.Static.Ns
            if length(find(Ns_static==nsens))>0
                
                Gs=[];
                Gm=[];
                yy=zeros(length(wk),SensorModel.Static.hn(nsens));
                for i=1:1:length(wk)
                    [yy(i,:),g,tp]=SensorModel.Static.hx{nsens}(Xk(i,:),SensorModel.Static.pos{nsens},SensorModel.Static.FOV{nsens});
                    Gs=vertcat(Gs,g);
                end
                % true measurements
                [hh,g,tp]=SensorModel.Static.hx{nsens}(TargetPos.truth{targid}(curk,:),SensorModel.Static.pos{nsens},SensorModel.Static.FOV{nsens});
                Gm=vertcat(Gm,g);
                
                yg=hh+sqrtm(SensorModel.Static.R{nsens})*randn(SensorModel.Static.hn(nsens),1);
                if length(yg)==2
                    x123=cos(yg(2));
                    y123=sin(yg(2));
                    yg(2)=atan2(y123,x123);
                end
                [ft,~]=GetCellOptions(further,'NoStaticSensorsFOV');
                if ft==true
                    ym=vertcat(ym,yg);
                    ytp=vertcat(ytp,tp);
                    Y=horzcat(Y,yy);
                    RR=blkdiag(RR,SensorModel.Static.R{nsens});
                else
                    if  sum(isnan(Gs))>0.75*length(Gs) || isnan(Gm)

                    else
                        ym=vertcat(ym,yg);
                        ytp=vertcat(ytp,tp);
                        Y=horzcat(Y,yy);
                        RR=blkdiag(RR,SensorModel.Static.R{nsens});
                    end
                end
            end
        end
    end
%     keyboard
    [st,~]=GetCellOptions(further,'NoMoveSensros');
    if st==false
        for nsens=1:SensorModel.Move.Ns
            if length(find(Ns_move==nsens))>0
                Gs=[];
                Gm=[];
                pos=SensorModel.Move.Grid2Abs(SensorModel.Move.Gridpos{nsens}(curk,1),SensorModel.Move.Gridpos{nsens}(curk,2),SensorModel);
                yy=zeros(length(wk),SensorModel.Move.hn(nsens));
                for i=1:1:length(wk)
                    [yy(i,:),gs,tp]=SensorModel.Move.hx{nsens}(Xk(i,:),pos,SensorModel.Move.FOV{nsens});
                    Gs=vertcat(Gs,gs);
                end
                % true measurements
                [hh,gm,tp]=SensorModel.Move.hx{nsens}(TargetPos.truth{targid}(curk,:),pos,SensorModel.Move.FOV{nsens});
                Gm=vertcat(Gm,gm);
                
                yg=hh+sqrtm(SensorModel.Move.R{nsens})*randn(SensorModel.Move.hn(nsens),1);
                if length(yg)==2
                    x123=cos(yg(2));
                    y123=sin(yg(2));
                    yg(2)=atan2(y123,x123);
                end
                if  sum(isnan(Gs))>0.75*length(Gs) || isnan(Gm)
                    
                else
                    ym=vertcat(ym,yg);
                    ytp=vertcat(ytp,tp);
                    Y=horzcat(Y,yy);
                    RR=blkdiag(RR,SensorModel.Move.R{nsens});
                end
            end
        end
    end
    
    
%     keyboard
    
    [st,~]=GetCellOptions(further,'pseudo');
    if st==true
        if isempty(Y)==0
            [mz,Pz]=MeanCov(Y,wk);
            try
            Pz=Pz+RR;
            catch
                keyboard
            end
            Pcc=CrossCov(Xk,xfk,Y,mz,wk);
            K=Pcc/Pz;
            
            xk=xfk;
            Pk=Pfk-K*Pz*K';
        else
            xk=xfk;
            Pk=Pfk;   
            flag(targid)=NaN;
        end
    else
        if isempty(ym)==0
            [mz,Pz]=MeanCov(Y,wk);
            Pz=Pz+RR;
            Pcc=CrossCov(Xk,xfk,Y,mz,wk);
            K=Pcc/Pz;
            
            d=ym(:)-mz(:);
            ind=find(ytp==1);
            for sd=1:1:length(ind)
                %                 qz=angle2quat(mz(ind(sd)),0,0);
                %                 qm=angle2quat(ym(ind(sd)),0,0);
                %                 d(sd)=quat2angle(quatmultiply(qm,quatinv(qz)));
                
%                 x123=cos(d(ind(sd)));
%                 y123=sin(d(ind(sd)));
%                 d(ind(sd))=atan2(y123,x123);

                if ym(ind(sd))<0
                    ym(ind(sd))=ym(ind(sd))+2*pi;
                end
                if mz(ind(sd))<0
                    mz(ind(sd))=mz(ind(sd))+2*pi;
                end
                d(ind(sd))=ym(ind(sd))-mz(ind(sd));
                x123=cos(d(ind(sd)));
                y123=sin(d(ind(sd)));
                d(ind(sd))=atan2(y123,x123);
            end

            xk=xfk+K*(d);
            Pk=Pfk-K*Pz*K';
        else
            xk=xfk;
            Pk=Pfk;    
            flag(targid)=NaN;
        end
    end
    
    
    
    if norm(xk(1:2)-xfk(1:2))>100 || min(sign(xk(1:2)))==-1 
        keyboard
    end
    
    TargetPos.xf{targid}(curk,:)=xk;
    TargetPos.Pf{targid}(curk,:)=reshape(Pk,1,TargetPos.fn(targid)^2);
    
    XKU{targid}=xk;
    PKU{targid}=Pk;
    
    
    
    %% Old Stuff.... Needs to be disregarde
    %     if isempty(SensorModel.Static.Task{k})
    %         indstatic=[];
    %     else
    %         indstatic=find(SensorModel.Static.Task{k}(:,2)==targid);
    %     end
    %     if isempty(SensorModel.Move.Task{k})
    %         indmove=[];
    %     else
    %         indmove=find(SensorModel.Move.Task{k}(:,2)==targid);
    %     end
    %     if isempty(indstatic) && isempty(indmove)
    %         ym=[];
    %         ytp=[];
    %         RR=[];
    %     else
    %         RR=[];
    %         ym=[];
    %         ytp=[];
    %         Y=[];
    %
    %         if sum(isempty(indstatic))==0
    %             for nsens=SensorModel.Static.Task{k}(indstatic,1)'
    %                 yy=zeros(length(wk),SensorModel.Static.hn(nsens));
    %                 for i=1:1:length(wk)
    %                     yy(i,:)=SensorModel.Static.hx{nsens}(Xk(i,:),SensorModel.Static.pos{nsens},SensorModel.Static.FOV{nsens});
    %
    %                 end
    %                 if strcmpi(type,'true')==1
    %                     [hh,~,tp]=SensorModel.Static.hx{nsens}(TargetPos.truth{targid}(k,:),SensorModel.Static.pos{nsens},SensorModel.Static.FOV{nsens});
    %                     ytp=vertcat(ytp,tp);
    %
    %                     yg=hh+sqrtm(SensorModel.Static.R{nsens})*randn(SensorModel.Static.hn(nsens),1);
    %                     if length(yg)==2
    %                         if yg(2)>pi
    %                             yg(2)=yg(2)-2*pi;
    %                         elseif yg(2)<-pi
    %                             yg(2)=yg(2)+2*pi;
    %                         end
    %                     end
    %                     ym=vertcat(ym,yg);
    %                 else
    %                     ym=vertcat(ym,NaN*ones(SensorModel.Static.hn(nsens),1));
    %                 end
    %                 Y=horzcat(Y,yy);
    %                 RR=blkdiag(RR,SensorModel.Static.R{nsens});
    %
    %             end
    %         end
    %         if sum(isempty(indmove))==0
    %             for nsens=SensorModel.Move.Task{k}(indmove,1)'
    %                 pos=SensorModel.Move.Grid2Abs(SensorModel.Move.Gridpos{nsens}(k,1),SensorModel.Move.Gridpos{nsens}(k,2),SensorModel);
    %                 yy=zeros(length(wk),SensorModel.Move.hn(nsens));
    %                 for i=1:1:length(wk)
    %                     try
    %                         yy(i,:)=SensorModel.Move.hx{nsens}(Xk(i,:),pos,SensorModel.Move.FOV{nsens});
    %                     catch
    %                         keyboard
    %                     end
    %                 end
    %
    %                 if strcmpi(type,'true')==1
    %
    %                     [hh,~,tp]=SensorModel.Move.hx{nsens}(TargetPos.truth{targid}(k,:),pos,SensorModel.Move.FOV{nsens});
    %                     ytp=vertcat(ytp,tp);
    %
    %                     yg=hh+sqrtm(SensorModel.Move.R{nsens})*randn(SensorModel.Move.hn(nsens),1);
    %                     if length(yg)==2
    %                         if yg(2)>pi
    %                             yg(2)=yg(2)-2*pi;
    %                         elseif yg(2)<-pi
    %                             yg(2)=yg(2)+2*pi;
    %                         end
    %                     end
    %                     ym=vertcat(ym,yg);
    %                 else
    %                     ym=vertcat(ym,NaN*ones(SensorModel.Move.hn(nsens),1));
    %                 end
    %                 Y=horzcat(Y,yy);
    %                 RR=blkdiag(RR,SensorModel.Move.R{nsens});
    %             end
    %         end
    %
    %     end
    %
    %
    %     if sum(isnan(ym))==0 && sum(isempty(ym))==0 && strcmpi(type,'true')==1
    %
    %         [mz,Pz]=MeanCov(Y,wk);
    %         Pz=Pz+RR;
    %         Pcc=CrossCov(Xk,xfk,Y,mz,wk);
    %         %         [xk,Pk]=KalmanUpdate(xfk,Pfk,mz,Pz,Pcc,ym);
    %         [xk,Pk]=KalmanUpdate_modf_angle(xfk,Pfk,mz,Pz,Pcc,ym,ytp);
    %
    %         if norm(xk(1:2)-xfk(1:2))>150
    %             keyboard
    %         end
    %         TargetPos.xf{targid}(k,:)=xk;
    %         TargetPos.Pf{targid}(k,:)=reshape(Pk,1,TargetPos.fn(targid)^2);
    %
    %     elseif  strcmpi(type,'pseudo')==1 && sum(isempty(ym))==0
    %
    %         [mz,Pz]=MeanCov(Y,wk);
    %         Pz=Pz+RR;
    %         Pcc=CrossCov(Xk,xfk,Y,mz,wk);
    %         %        try
    %         [xk,Pk]=KalmanUpdate(xfk,Pfk,mz,Pz,Pcc,mz); % ym is replaced by mz so there is bo innovation
    %         %         [xk,Pk]=KalmanUpdate_modf_angle(xfk,Pfk,mz,Pz,Pcc,mz);
    %
    %         %        catch
    %         %                    keyboard
    %         %        end
    %         if norm(xk(1:2)-xfk(1:2))>150
    %             keyboard
    %         end
    %         TargetPos.xf{targid}(k,:)=xk;
    %         TargetPos.Pf{targid}(k,:)=reshape(Pk,1,TargetPos.fn(targid)^2);
    %
    %
    %     elseif strcmp(lower(type),'true')==1
    %         %         disp('dude')
    %         %i.e. when true meas update has to happen buyt there is no
    %         %measurement assigned just do nothing
    %         %         TargetPos.xf{targid}(k,:)=xk;
    %         %         TargetPos.Pf{targid}(k,:)=reshape(Pk,1,TargetPos.fn(targid)^2);
    %
    %     end
end



end

