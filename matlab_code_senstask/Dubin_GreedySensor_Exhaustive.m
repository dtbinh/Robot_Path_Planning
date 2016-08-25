function SensorModel=Dubin_GreedySensor_Exhaustive(curk,futk,Time,SensorModel,TargetPos,further)
% curk-1 is the initial pos of sensor and is known, planning from curk to
% futk,  SensorModel.Move.Gridpos{ss}(curk-1,:) is known.
% The system is up-to-date at curk-1
% Load the path into SensorModel.Move.Gridpos{ss}(curk+1 to futk,:)=[xg,yg,dir] 1
% is up, 2 is right

Tn=curk-1:1:futk;
Tvec=Time.Tvec(Tn);


%               Mutual info is calculated with the prior of the first step
               [TargetPos_0,~,~]=TargetPriorProp(1:TargetPos.Ntar,curk-1,curk,TargetPos,Time,{});
                X0=cell(TargetPos.Ntar,1);
                P0=cell(TargetPos.Ntar,1);
                for i=1:TargetPos.Ntar
                  X0{i}=TargetPos_0.xf{i}(curk-1,:)';
                  P0{i}=reshape(TargetPos_0.Pf{i}(curk-1,:),TargetPos_0.fn(i),TargetPos_0.fn(i));
                end
            
                for k=2:1:length(Tn)
                    for statens=1:1:SensorModel.Static.Ns
                        SensorModel.Static.Task{Tn(k)}=vertcat(SensorModel.Static.Task{Tn(k)},[statens*ones(TargetPos.Ntar,1),[1:TargetPos.Ntar]']);
                    end
                end
                
                
                
% get the 2-sigma grid nodes
TargSensNodes_index=[];
for targid=1:TargetPos.Ntar
    for k=Tn
        Xfk=XK1{targid}(k,1:2)';
        Pfk=reshape(PK1{targid}(k,:),TargetPos.fn(targid),TargetPos.fn(targid));
        Pfk=Pfk(1:2,1:2);
        XYellipse=sigmaEllipse(Xfk,Pfk,2);
        
        IN = inpolygon(SensorModel.Move.XYgridpoints(:,1),SensorModel.Move.XYgridpoints(:,2),XYellipse(:,1),XYellipse(:,2));
        TargSensNodes_index=horzcat(TargSensNodes_index,IN(:)');
    end
end
TargSensNodes_index=unique(TargSensNodes_index);
TargSensNodesPos=SensorModel.Move.XYgridpoints(TargSensNodes_index,:);

for b=1:1:length(Tn)
    for l=1:SensorModel.Static.Ns
        SensorModel.Static.Task{Tn(b)}=vertcat(SensorModel.Static.Task{Tn(b)},[l*ones(TargetPos.Nrar,1),1:TargetPos.Ntar]);  % task [sensor,targ] then make the measurement
    end
end

for movtarg=1:SensorModel.Move.Ns
    UAVnodepos0=SensorModel.Move.Grid2Abs(SensorModel.Move.Gridpos{movtarg}(curk-1,1),SensorModel.Move.Gridpos{movtarg}(curk-1,2),SensorModel);
    UAVnodepos0_index=SensorModel.Move.Grid2IndexVec(SensorModel.Move.Gridpos{movtarg}(curk-1,1),SensorModel.Move.Gridpos{movtarg}(curk-1,2),SensorModel);
    
    %     [~,UAVnode0ind]=min(sqrt(sum((repmat(xsenspos0,size(obj.SimProps.XY,1),1)-obj.SimProps.XY).^2,2)));
    %     UAVnode0=obj.SimProps.XYindex(UAVnode0ind);
    
    StateTree=GenrateTreefromGrid(SensorModel.Move.XYgridpoints,SensorModel.Move.XYADJindex,   length(Tn)  ,UAVnodepos0_index);
    
    P=Generatepaths(StateTree,UAVnodepos0_index,1);
    
    TargSensNodes=setdiff(TargSensNodes,UAVnode0);
    Pp=PrunePaths(P,TargSensNodes);
    
    if isempty(Pp) % the nodes are not reachable
        
        %if the number of steps cannot reach the nearest ifo node,
        %then form a triangle [UAV, nearest node, farthest node],
        % prune trajectories to be in this triangle
        TargSensNodesPos=obj.SimProps.XY(TargSensNodes,:);
        dis=repmat(xsenspos0,size(TargSensNodesPos,1),1)-TargSensNodesPos;
        [~,mindis]=min(sqrt(sum(dis.^2,2)));
        [~,maxdis]=max(sqrt(sum(dis.^2,2)));
        
        nearestnode=TargSensNodes(mindis);
        farthestnode=TargSensNodes(maxdis);
        
        pruntriagregion=obj.SimProps.XY([UAVnode0ind,nearestnode,farthestnode],:);
        
        IN = inpolygon(obj.SimProps.XY(:,1),obj.SimProps.XY(:,2),pruntriagregion(:,1),pruntriagregion(:,2));
        pruntriagnodes=obj.SimProps.XYindex(IN)';
        
        Pp=PrunePaths(P,pruntriagnodes);
    end
    
    
    % Debug plotting
    %                 keyboard
    
    %                 obj.plot_Dynamic_scenario(Tvec(1));
    %                 hold on
    %                 for jj=1:1:size(Pp,1)
    %                     for ij=1:1:length(Tvec)
    %                    plot(obj.SimProps.XY(Pp(jj,ij),1),obj.SimProps.XY(Pp(jj,ij),2),'ro' )
    %                     end
    %                 end
    %                 for jj=1:1:length(TargSensNodes)
    %                     plot(obj.SimProps.XY(TargSensNodes(jj),1),obj.SimProps.XY(TargSensNodes(jj),2),'bo','MarkerSize',8 )
    %                 end
    %                 plot3(obj.SimProps.XY(:,1),obj.SimProps.XY(:,2),obj.SimProps.XYindex,'.')
    %
    %
    
    
    % METHOD 1:   doing sequential estimation while taking all the static
    % sensors
    
    
    Mi_paths=zeros(size(Pp,1),1);
    
    for b=1:1:length(Tn)
        SensorModel.Move.Task{Tn(b)}=vertcat(SensorModel.Move.Task{Tn(b)},[movtarg*ones(TargetPos.Ntar,1),1:TargetPos.Ntar]);  % task [sensor,targ] then make the measurement
    end
    
    tic
    MIpaths=zeros(size(Pp,1));
    for i=1:1:size(Pp,1)
        TargetPos_copy=TargetPos;
        for b=1:1:length(Tn)
            SensorModel.Move.Gridpos{movtarg}(Tn(b),:)=[SensorModel.Move.IndexVec2Grid(Pp(i,b),SensorModel),1];%[SensorModel.Move.XYgridpoints(Pp(i,b),:),1];
            if b==1
                [TargetPos_copy,XK1,PK1]=TargetPriorProp(ntars,curk-1,Tn(b),TargetPos_copy,Time,further);
            else
                [TargetPos_copy,XK1,PK1]=TargetPriorProp(ntars,Tn(b-1),Tn(b),TargetPos_copy,Time,further);
            end
            [TargetPos_copy,XKF,PKF,XKU,PKU]=TargetPostUpdate(ntars,Tn(b),type,SensorModel,TargetPos_copy,Time,further);
            
            for tt=1:1:TargetPos.Ntar
                if det(PKF{tt})/det(PKU{tt})<1.01  %Pf/Pu
                    MIpaths(i)=MIpaths(i)+0;
                else
                    MIpaths(i)=MIpaths(i)+0.5*log(det(PKF{tt})/det(PKU{tt}));
                end
            end
            
        end
    end
    toc
    %                 keyboard
    [~,ind]=max(Mi_paths);

    for b=1:1:length(Tn)
        SensorModel.Move.Gridpos{movtarg}(Tn(b),:)=[SensorModel.Move.IndexVec2Grid(Pp(ind,b),SensorModel),1];%[SensorModel.Move.XYgridpoints(Pp(i,b),:),1];
    end
    
    
end
end


