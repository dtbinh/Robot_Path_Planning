function [M,J,SensorModel]=dubin_DP(curk,futk,Time,SensorModel,TargetPos)


Tk=curk:1:futk-1;



tic
for k=1:1:length(Tk)
    
    for i=1:1:size(SensorModel.Move.Xgrid,1)
        
        for j=1:size(SensorModel.Move.Xgrid,2)
            m=zeros(1,TargetPos.Ntar );
            for targid=1:1:TargetPos.Ntar
                xf=TargetPos.xf{targid}(Tk(k),:)';
                Pf=reshape(TargetPos.Pf{targid}(Tk(k),:),TargetPos.fn(targid),TargetPos.fn(targid));
                if max(eig(sqrtm(Pf)))<=TargetPos.thresh
                    disp(strcat('target id : ',num2str(targid),' information ignored'))
                    continue;
                end
                SensorModel.Static.Task{Tk(k)}=[];
                Tsk=SimpleSensortask_static(xf,Pf,SensorModel);
                if isempty(Tsk)==0
                    SensorModel.Static.Task{Tk(k)}=[Tsk,targid*ones(length(Tsk),1)];
                end
                %                 keyboard
                SensorModel.Move.Task{Tk(k)}=[];
                Tsk=SimpleSensortask_move(i,j,[],xf,Pf,SensorModel);
                if isempty(Tsk)==0
                    SensorModel.Move.Task{Tk(k)}=[Tsk,targid*ones(length(Tsk),1)];
                end
                %                 if isempty(Tsk)==0
                %                     SensorModel.Move.Task{Tk(k)}(2:end,:)=[];
                %                 end
                %                  keyboard
                SensorModel_copy=SensorModel;
                for ss=1:1:SensorModel_copy.Move.Ns
                    SensorModel_copy.Move.Gridpos{ss}(Tk(k),1:2)=[i,j];
                end
                %                     [i,j]
                %                     if i==1 && j==4 && targid == 2
                %                         keyboard
                %                     end
                [~,~,Pu]=TargetPostUpdate(targid,Tk(k),'pseudo',SensorModel_copy,TargetPos);
                %                 SensorModel.Move.Task{Tk(k)}=[];
%                 trace(sqrtm(Pf))
                if max(max(isempty(Pu)))==1 || trace(sqrtm(Pf))<5
                    m(targid)=0;
                else
                    m(targid)=0.5*log(det(Pf)/det(Pu));
                end
            end
            SensorModel.InfoMap{Tk(k)}(i,j)=sum(m);
            
        end
    end
    figure(5)
    mesh(SensorModel.Move.Xgrid,SensorModel.Move.Ygrid, SensorModel.InfoMap{Tk(k)})
    %     keyboard
    
    [TargetPos,~,~]=TargetPriorProp(1:1:TargetPos.Ntar,Tk(k),TargetPos,Time);
    
    
end
tt=toc;
disp(strcat('Info map generation took = ',num2str(tt)))
%% Run DP sequentially centralized
%  keyboard
tic
% setup the priority order
Ps=1:1:SensorModel.Move.Ns;
% initialize DP and run
J=cell(SensorModel.Move.Ns,length(Tk));
M=cell(SensorModel.Move.Ns,length(Tk));
for i=1:SensorModel.Move.Ns
    for j=1:length(Tk)
        J{i,j}=zeros(size(SensorModel.Move.Xgrid,1),size(SensorModel.Move.Xgrid,2),4);
        J{i,j}(:,:,1)=SensorModel.InfoMap{Tk(end)}';
        J{i,j}(:,:,2)=SensorModel.InfoMap{Tk(end)}';
        J{i,j}(:,:,3)=SensorModel.InfoMap{Tk(end)}';
        J{i,j}(:,:,4)=SensorModel.InfoMap{Tk(end)}';
        M{i,j}=zeros(size(SensorModel.Move.Xgrid,1),size(SensorModel.Move.Xgrid,2),4);
    end
end
Ic=zeros(length(Tk),1); % info collected by previous sensors
Xpath=cell(SensorModel.Move.Ns,1);
%

A=1e3;% information weighting
B=1e2;
for ss=1:1:SensorModel.Move.Ns
    for k=(length(Tk)-1):-1:1
        for i=1:1:size(SensorModel.Move.Xgrid,1)
            for j=1:1:size(SensorModel.Move.Xgrid,2)
                for dirn=1:1:4
                    JJ=1e25;
                    umin=1;
                    for U=1:1:SensorModel.Move.GridU
                        [Xf,V]=SensorModel.Move.GridMotion(i,j,dirn,U);
                        if isnan(Xf)==1
                            continue
                        end      % +1e2*collision_penalty(Xpath,i,j,size(Xpath,1),k)

                        Jc=-A*(Ic(k)+SensorModel.InfoMap{Tk(k)}(i,j))+1*V/10+10000*collision_penalty(Xpath,i,j,ss-1,k)+J{Ps(ss),k+1}(Xf(1),Xf(2),Xf(3));
                        if Jc<JJ
                            JJ=Jc;
                            umin=U;
                        end
                    end
                    if umin==0
                        keybvoard
                    end
                    J{Ps(ss),k}(i,j,dirn)=JJ;
                    M{Ps(ss),k}(i,j,dirn)=umin;
                end
            end
        end
        
    end % time step end
    %loading this path
    Xpath{ss}=SensorModel.Move.Gridpos{ss}(curk,:);
    for kk=1:1:(length(Tk))-1
        c= Xpath{ss}(kk,:);
        Ic(kk)=Ic(kk)+SensorModel.InfoMap{Tk(kk)}(c(1),c(2));
        u=M{ss,kk}(c(1),c(2),c(3));
        try
            [Xf,~]=SensorModel.Move.GridMotion(c(1),c(2),c(3),u);
        catch
            keyboard
        end
        Xpath{ss}(kk+1,:)=Xf;
        
        
    end
    
    % calculate the path and its information information along the optimal path
    
    
    
end % end for each sensor

tt=toc;
disp(strcat('Optimal Traj generation took = ',num2str(tt)))
%   keyboard
%% load the optimal path into the SensorModel.Move.GridOptiTraj

for kk=1:1:(length(Tk))
    for ss=1:1:SensorModel.Move.Ns
        SensorModel.Move.DP_MU{ss,Tk(kk)}=M{ss,kk};
    end
end


