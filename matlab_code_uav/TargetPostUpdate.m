function [TargetPos,xk,Pk]=TargetPostUpdate(ntars,k,type,SensorModel,TargetPos)
% update the target at the same time k
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

xk=[];
Pk=[];
for targid=ntars
    
    xfk=TargetPos.xf{targid}(k,:)';
    Pfk=reshape(TargetPos.Pf{targid}(k,:),TargetPos.fn(targid),TargetPos.fn(targid));
    [Xk,wk]=qd_pts(xfk,Pfk);
    
    %Now assemble the R and ym according to the sensor tasking
    % Do this both for static and movable sensors
    if isempty(SensorModel.Static.Task{k})
        indstatic=[];
    else
        indstatic=find(SensorModel.Static.Task{k}(:,2)==targid);
    end
    if isempty(SensorModel.Move.Task{k})
        indmove=[];
    else
        indmove=find(SensorModel.Move.Task{k}(:,2)==targid);
    end
    if isempty(indstatic) && isempty(indmove)
        ym=[];
        ytp=[];
        RR=[];
    else
        RR=[];
        ym=[];
        ytp=[];
        Y=[];
        
        if sum(isempty(indstatic))==0
            for nsens=SensorModel.Static.Task{k}(indstatic,1)'
                yy=zeros(length(wk),SensorModel.Static.hn(nsens));
                for i=1:1:length(wk)
                    yy(i,:)=SensorModel.Static.hx{nsens}(Xk(i,:),SensorModel.Static.pos{nsens},SensorModel.Static.FOV{nsens});
                    
                end
                if strcmpi(type,'true')==1
                    [hh,~,tp]=SensorModel.Static.hx{nsens}(TargetPos.truth{targid}(k,:),SensorModel.Static.pos{nsens},SensorModel.Static.FOV{nsens});
                    ytp=vertcat(ytp,tp);
                    
                    yg=hh+sqrtm(SensorModel.Static.R{nsens})*randn(SensorModel.Static.hn(nsens),1);
                    if length(yg)==2
                        if yg(2)>pi
                            yg(2)=yg(2)-2*pi;
                        elseif yg(2)<-pi
                            yg(2)=yg(2)+2*pi;
                        end
                    end
                    ym=vertcat(ym,yg);
                else
                    ym=vertcat(ym,NaN*ones(SensorModel.Static.hn(nsens),1));
                end
                Y=horzcat(Y,yy);
                RR=blkdiag(RR,SensorModel.Static.R{nsens});
                
            end
        end
        if sum(isempty(indmove))==0
            for nsens=SensorModel.Move.Task{k}(indmove,1)'
                pos=SensorModel.Move.Grid2Abs(SensorModel.Move.Gridpos{nsens}(k,1),SensorModel.Move.Gridpos{nsens}(k,2),SensorModel);
                yy=zeros(length(wk),SensorModel.Move.hn(nsens));
                for i=1:1:length(wk)
                    try
                        yy(i,:)=SensorModel.Move.hx{nsens}(Xk(i,:),pos,SensorModel.Move.FOV{nsens});
                    catch
                        keyboard
                    end
                end
                
                if strcmpi(type,'true')==1
                    
                    [hh,~,tp]=SensorModel.Move.hx{nsens}(TargetPos.truth{targid}(k,:),pos,SensorModel.Move.FOV{nsens});
                    ytp=vertcat(ytp,tp);
                    
                    yg=hh+sqrtm(SensorModel.Move.R{nsens})*randn(SensorModel.Move.hn(nsens),1);
                    if length(yg)==2
                        if yg(2)>pi
                            yg(2)=yg(2)-2*pi;
                        elseif yg(2)<-pi
                            yg(2)=yg(2)+2*pi;
                        end
                    end
                    ym=vertcat(ym,yg);
                else
                    ym=vertcat(ym,NaN*ones(SensorModel.Move.hn(nsens),1));
                end
                Y=horzcat(Y,yy);
                RR=blkdiag(RR,SensorModel.Move.R{nsens});
            end
        end
        
    end
    

    if sum(isnan(ym))==0 && sum(isempty(ym))==0 && strcmpi(type,'true')==1
        
        [mz,Pz]=MeanCov(Y,wk);
        Pz=Pz+RR;
        Pcc=CrossCov(Xk,xfk,Y,mz,wk);
        %         [xk,Pk]=KalmanUpdate(xfk,Pfk,mz,Pz,Pcc,ym);
        [xk,Pk]=KalmanUpdate_modf_angle(xfk,Pfk,mz,Pz,Pcc,ym,ytp);
        
        if norm(xk(1:2)-xfk(1:2))>150
            keyboard
        end
        TargetPos.xf{targid}(k,:)=xk;
        TargetPos.Pf{targid}(k,:)=reshape(Pk,1,TargetPos.fn(targid)^2);
        
    elseif  strcmpi(type,'pseudo')==1 && sum(isempty(ym))==0
        
        [mz,Pz]=MeanCov(Y,wk);
        Pz=Pz+RR;
        Pcc=CrossCov(Xk,xfk,Y,mz,wk);
        %        try
        [xk,Pk]=KalmanUpdate(xfk,Pfk,mz,Pz,Pcc,mz); % ym is replaced by mz so there is bo innovation
        %         [xk,Pk]=KalmanUpdate_modf_angle(xfk,Pfk,mz,Pz,Pcc,mz);
        
        %        catch
        %                    keyboard
        %        end
        if norm(xk(1:2)-xfk(1:2))>150
            keyboard
        end
        TargetPos.xf{targid}(k,:)=xk;
        TargetPos.Pf{targid}(k,:)=reshape(Pk,1,TargetPos.fn(targid)^2);
        
        
    elseif strcmp(lower(type),'true')==1
%         disp('dude')
        %i.e. when true meas update has to happen buyt there is no
        %measurement assigned just do nothing
        %         TargetPos.xf{targid}(k,:)=xk;
        %         TargetPos.Pf{targid}(k,:)=reshape(Pk,1,TargetPos.fn(targid)^2);

    end
end



end

