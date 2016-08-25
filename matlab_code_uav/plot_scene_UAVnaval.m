function plot_scene_UAVnaval(TargetPos,Trajectories,SensorModel,curk)
hold on
Xc=0:TargetPos.Trafficdx:TargetPos.Traffic_xlim;
Yc=0:TargetPos.Trafficdx:TargetPos.Traffic_ylim;

plot(TargetPos.Pu(:,1),TargetPos.Pu(:,2),'k',TargetPos.Pd(:,1),TargetPos.Pd(:,2),'k',TargetPos.Pc(:,1),TargetPos.Pc(:,2),'k')
fill(TargetPos.Pu(:,1),TargetPos.Pu(:,2),'r')
fill(TargetPos.Pd(:,1),TargetPos.Pd(:,2),'r')
fill(TargetPos.Pc(:,1),TargetPos.Pc(:,2),'b')
alpha 0.2
plot_1sig_contour_multi(TargetPos,curk)
%     for i=1:1:length(Xc)
%         plot([Xc(i),Xc(i)],[0,TargetPos.Traffic_xlim],'r','linewidth',1)
%         plot([0,TargetPos.Traffic_ylim],[Yc(i),Yc(i)],'r','linewidth',1)
%     end
%plot sensors and their FOVs
if strcmp(SensorModel.Static.Mode,'off')
    % do nothing
elseif strcmp(SensorModel.Static.Mode,'passive') % just plot the radar locations
    for ij=1:1:SensorModel.Static.Ns
        plot(SensorModel.Static.pos{ij}(1),SensorModel.Static.pos{ij}(2),'ks','MarkerSize',6,'linewidth',3)
    end
elseif strcmp(SensorModel.Static.Mode,'active')
    for ij=1:1:SensorModel.Static.Ns
        plot(SensorModel.Static.pos{ij}(1),SensorModel.Static.pos{ij}(2),'ks','MarkerSize',6,'linewidth',3)
        plot_2Dcone(SensorModel.Static.pos{ij}(1),SensorModel.Static.pos{ij}(2),SensorModel.Static.FOV{ij}(1),SensorModel.Static.FOV{ij}(2),SensorModel.Static.FOV{ij}(3),'g','fill')
    end
end

%     keyboard
if strcmp(SensorModel.Move.Mode,'off')
    
    % do nothing
% elseif strcmp(SensorModel.Move.Mode,'passive') 
%     pos=SensorModel.Move.Grid2Abs(SensorModel.Move.Gridpos{ss}(curk,1),SensorModel.Move.Gridpos{ss}(curk,2),SensorModel);
%     plot_triangle(pos(1),pos(2),6,1.5,0,SensorModel.Move.trajcol{ss},'no')
elseif strcmp(SensorModel.Move.Mode,'active') || strcmp(SensorModel.Move.Mode,'passive') 
%     keyboard
    
    for sss=1:1:SensorModel.Move.Ns
        plot(Trajectories{sss}(:,1),Trajectories{sss}(:,2),SensorModel.Move.trajcol{sss})
        pos=SensorModel.Move.Grid2Abs(SensorModel.Move.Gridpos{sss}(curk,1),SensorModel.Move.Gridpos{sss}(curk,2),SensorModel);
        dir=SensorModel.Move.Gridpos{sss}(curk,3);
        if dir==1
            th=-pi/2;
        elseif dir==2
            th=0;
        elseif dir==3
            th=pi/2;
        elseif dir==4
            th=pi;
        else
            th=Nan;
        end
        plot_triangle(pos(1),pos(2),6,1.5,th,SensorModel.Move.trajcol{sss},'no')
        plot_ellipse(pos(1),pos(2),SensorModel.Move.FOV{sss}(2),SensorModel.Move.FOV{sss}(2),0,SensorModel.Move.trajcol{sss},'fill')
    end
end


axis([0 TargetPos.Traffic_xlim 0 TargetPos.Traffic_ylim])

