function SensorModel=dubin_FDP(curk,futk,Time,SensorModel,TargetPos,further)
% Implement forward DP
% curk-1 is considered to be the latest updated
% find the positions for curk to futk


Tn=curk-1:1:futk ;
Tvec=Time.Tvec(Tn);

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

for movsens=1:1:SensorModel.Move.Ns
    disp(['FDP for moving sensor = ',num2str(movsens)])
    %                 keyboard
    %
    xsens_gridpos0=SensorModel.Move.Gridpos{movsens}(curk-1,:);
    
    UAVnode_abspos0=SensorModel.Move.Grid2Abs(xsens_gridpos0(1),xsens_gridpos0(2),SensorModel);
    UAVnode_indexpos0=SensorModel.Move.Grid2IndexVec(xsens_gridpos0(1),xsens_gridpos0(2),SensorModel);
    
    
    StateTree=GenrateTreefromGrid(SensorModel.Move.XYgridpoints,SensorModel.Move.XYADJindex,  length( Tn)  ,UAVnode_indexpos0);
    
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
    
    % StateTree{Tk}{i,1} is node i
    % StateTree{Tk}{i,2} is node i, left parent nodes
    % StateTree{Tk}{i,3} is node i, right child nodes
    
    % Save the optimum gain
    % StateTree{Tk}{i,4} is node i, Optimum(gain) mean
    % StateTree{Tk}{i,5} is node i, Optimum(gain) Cov
    % StateTree{Tk}{i,6} is node i, Optimum(gain) Info
    % StateTree{Tk}{i,7} is node i, previous parent optimum node: to
    % follwo back and get the whole trajectory
    
    % info is calculated from the P0 i.e. at start
    % at the end nodes, find the max info na dwork your way
    % back to the starting node to get the
    
    TargetPos_copy=TargetPos;
    
    
    for k=1:1:length(Tn)
        if k>1
            %                         keyboard
            SensorModel.Move.Task{Tn(k)}=vertcat(SensorModel.Move.Task{Tn(k)},[movsens*ones(TargetPos.Ntar,1),[1:TargetPos.Ntar]']); % task [sensor,targ]
        end
        
        for i=1:1:size(StateTree{k},1)
            if k==1
                StateTree{k}{i,4}=TargetPos.xf;
                StateTree{k}{i,5}=TargetPos.Pf;
                StateTree{k}{i,6}=0;
                StateTree{k}{i,7}=[];
            else
                
                
                
                currnode=StateTree{k}{i,1};
                MIcurr=-1e5;
                for parentnode=StateTree{k}{i,2}
                    
                    for s=1:1:size(StateTree{k-1},1)
                        if StateTree{k-1}{s,1}==parentnode
                            break
                        end
                    end
                    TargetPos_copy.xf=StateTree{k-1}{s,4};
                    TargetPos_copy.Pf=StateTree{k-1}{s,5};
                    
                    MIparent=StateTree{k-1}{s,6};
                    
                    % Do propagation from parent node at Tvec()
                    [TargetPos_copy,~,~]=TargetPriorProp(1:TargetPos.Ntar,Tn(k-1),Tn(k),TargetPos_copy,Time,{});
                    %                                 [curr_mu,curr_Pk]=DynProp_FDP(    obj     ,parent_mu,      parent_cov,     1:obj.TargProps.NumbTargets,Tvec(k-1)       );  % from k-1 to k
                    
                    % force the UAV sensor to current node and do measurement update
                    SensorModel.Move.Gridpos{movsens}(Tn(k),:)=[SensorModel.Move.IndexVec2Grid(currnode,SensorModel),1];%
                    
                    % task [sensor,targ] then make the measurement
                    [TargetPos_copy,XKF,PKF,XKU,PKU,flag]=TargetPostUpdate(1:TargetPos.Ntar,Tn(k),[],SensorModel,TargetPos_copy,Time,{'pseudo'});
                    %                                 obj.SensProps.SensorPos_allk{movtarg}(Tvec(k),:) =obj.SimProps.XY(currnode,:);
                    %                                 [curr_mu,curr_Pk]=MeasUpdate_FDP(obj,curr_mu,curr_Pk,   1:obj.TargProps.NumbTargets,  sensids,   Tvec(k));
                    %
                    MI=MI_FDP(PKF,PKU,flag)+MIparent;
                    %                                 MI=MI_FDP(P0,PKU,flag);
                    
                    if MI>MIcurr
                        StateTree{k}{i,4}=TargetPos_copy.xf;
                        StateTree{k}{i,5}=TargetPos_copy.Pf;
                        StateTree{k}{i,6}=MI;
                        StateTree{k}{i,7}=parentnode;
                        MIcurr=MI;
                    end
                end
                
            end
        end
    end
    %                 keyboard
    % Now get the Trajectory for UAV in terms of indices
    Trajvec=zeros(1,length(Tn));
    for k=length(Tn):-1:2
        % for the final step find the maximum MI node
        if k==length(Tn)
            m=-1e5;
            pp=0;
            for i=1:1:size(StateTree{k},1)
                if StateTree{k}{i,6}>m
                    pp=i;
                    m=StateTree{k}{i,6};
                end
            end
            Trajvec(k)=StateTree{k}{pp,1};
            Trajvec(k-1)=StateTree{k}{pp,7};
            continue
        end
        
        for i=1:1:size(StateTree{k},1)
            if StateTree{k}{i,1}==Trajvec(k)
                Trajvec(k-1)=StateTree{k}{i,7};
            end
        end
        %
        %                     for ij=1:1:length(Trajvec)
        %                        SensorModel.Move.Gridpos{movsens}(Tn(ij),:) =[SensorModel.Move.XYgridpos (Trajvec(ij) ),1];
        %                     end
        
    end
    
    
    
    %                 Debug plotting
    %                 keyboard
    %                 figure(3)
    %                 plot_scene_UAVnaval_onlygrid(TargetPos,SensorModel,curk-1)
    %                 hold on
    %                 for ij=1:1:length(Trajvec)
    %                     plot(SensorModel.Move.XYgridpoints(Trajvec(ij),1),SensorModel.Move.XYgridpoints(Trajvec(ij),2),'ro' )
    %                 end
    
    
    
    for k=1:1:length(Tn)
        try
            SensorModel.Move.Gridpos{movsens}(Tn(k),:)=[SensorModel.Move.XYgridpos(Trajvec(k),:),1];
        catch
            keyboard
        end
    end
    
end

end

