function StateTree=dubinFDP_InnerloopParallel(i,k,movsens,P0,Tn,StateTree,TargetPos,SensorModel,Time)

TargetPos_copy=TargetPos;

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
            
            
            % Do propagation from parent node at Tvec()
            [TargetPos_copy,~,~]=TargetPriorProp(1:TargetPos.Ntar,Tn(k-1),Tn(k),TargetPos_copy,Time,{});
            %                                 [curr_mu,curr_Pk]=DynProp_FDP(    obj     ,parent_mu,      parent_cov,     1:obj.TargProps.NumbTargets,Tvec(k-1)       );  % from k-1 to k
            
            % force the UAV sensor to current node and do measurement update
            SensorModel.Move.Gridpos{movsens}(Tn(k),:)=[SensorModel.Move.IndexVec2Grid(currnode,SensorModel),1];%
            
            % task [sensor,targ] then make the measurement
            [TargetPos_copy,XKF,PKF,XKU,PKU]=TargetPostUpdate(1:TargetPos.Ntar,Tn(k),[],SensorModel,TargetPos_copy,Time,{'pseudo'});
            %                                 obj.SensProps.SensorPos_allk{movtarg}(Tvec(k),:) =obj.SimProps.XY(currnode,:);
            %                                 [curr_mu,curr_Pk]=MeasUpdate_FDP(obj,curr_mu,curr_Pk,   1:obj.TargProps.NumbTargets,  sensids,   Tvec(k));
            
            MI=MI_FDP(P0,PKU);
            
            if MI>MIcurr
                StateTree{k}{i,4}=TargetPos_copy.xf;
                StateTree{k}{i,5}=TargetPos_copy.Pf;
                StateTree{k}{i,6}=MI;
                StateTree{k}{i,7}=parentnode;
                MIcurr=MI;
            end
        end
        
    end
