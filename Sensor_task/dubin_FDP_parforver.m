function SensorModel=dubin_FDP_parforver(curk,futk,Time,SensorModel,TargetPos,further)
            % Implement forward DP
            % curk-1 is considered to be the latest updated
            % find the positions for curk to futk

            
            Tn=curk-1:1:futk ;
            Tvec=Time.Tvec(Tn);
            % GOING TO FORCE TASK ALL SENSORS TO SEE EVERYTHING
%             obj=BruteForceTasking(obj,Tvec);
            
%              % First Copy the sensor properties from Tk all the way upto Tk1
%             for k=Tvec(1:end-1)
%                 obj=CopySensorProps2NextTime(obj,k,k+1,1:obj.SensProps.NumbSensors);
%             end
            
%             % get the moving and stationary sensor ids
%             [SensIDstat,SensIDmove]=GetMovingSensorIDs(obj,Tk);
% %             
%             
%             % Now doing it for each moving sensor in order
%             sensids=SensIDstat;

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
                    
                    ST=cell(size(StateTree{k},1),1);
                    parfor i=1:size(StateTree{k},1)
                        ST{i}=dubinFDP_InnerloopParallel(i,k,movsens,P0,Tn,StateTree,TargetPos,SensorModel,Time);
                    end
                    
                    for i=1:1:size(StateTree{k},1)
                        StateTree{k}{i,4}=ST{i}{k}{i,4};
                        StateTree{k}{i,5}=ST{i}{k}{i,5};
                        StateTree{k}{i,6}=ST{i}{k}{i,6};
                        StateTree{k}{i,7}=ST{i}{k}{i,7};
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
        
function MI=MI_FDP(P0,PKU)

MI=0;
for ntar=1:1:length(P0)
   

                if det(P0{ntar})/det(PKU{ntar})<1.01  %Pf/Pu
                    MI=MI+0;
                else
                    MI=MI+0.5*log(det(P0{ntar})/det(PKU{ntar}));
                end

    
end





end