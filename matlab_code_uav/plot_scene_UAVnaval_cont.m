function plot_scene_UAVnaval_cont(TargetPos,Trajectories,SensorModel,nk,ck,curk)
%nk is the number steps to devide curk-1:curk interval
hold on
Xc=0:TargetPos.Trafficdx:TargetPos.Traffic_xlim;
Yc=0:TargetPos.Trafficdx:TargetPos.Traffic_ylim;

%plot the regions trafficability
plot(TargetPos.Pu(:,1),TargetPos.Pu(:,2),'k',TargetPos.Pd(:,1),TargetPos.Pd(:,2),'k',TargetPos.Pc(:,1),TargetPos.Pc(:,2),'k')
fill(TargetPos.Pu(:,1),TargetPos.Pu(:,2),'r')
fill(TargetPos.Pd(:,1),TargetPos.Pd(:,2),'r')
fill(TargetPos.Pc(:,1),TargetPos.Pc(:,2),'b')
alpha 0.2

%     for i=1:1:length(Xc)
%         plot([Xc(i),Xc(i)],[0,TargetPos.Traffic_xlim],'r','linewidth',1)
%         plot([0,TargetPos.Traffic_ylim],[Yc(i),Yc(i)],'r','linewidth',1)
%     end
%plot sensors and their FOVs
for ij=1:1:SensorModel.Static.Ns
    plot(SensorModel.Static.pos{ij}(1),SensorModel.Static.pos{ij}(2),'ks','MarkerSize',6,'linewidth',3)
    plot_2Dcone(SensorModel.Static.pos{ij}(1),SensorModel.Static.pos{ij}(2),SensorModel.Static.FOV{ij}(1),SensorModel.Static.FOV{ij}(2),SensorModel.Static.FOV{ij}(3),'g','fill')
end


%     keyboard

% for ik=1:1:nk
    
    if ck>nk/2 
        plot_1sig_contour_multi(TargetPos,curk)
    else
        if curk==1
            plot_1sig_contour_multi(TargetPos,curk)
        else
            plot_1sig_contour_multi(TargetPos,curk-1)
        end
    end
    for ss=1:1:SensorModel.Move.Ns
        ind=find(Trajectories{ss}(:,4)==curk);
        X=Trajectories{ss}(ind,1);
        Y=Trajectories{ss}(ind,2);
        phi=Trajectories{ss}(ind,3);
        
        Xr=interp1(1:1:length(X),X,linspace(1,length(X),nk));
        Yr=interp1(1:1:length(Y),Y,linspace(1,length(Y),nk));
        phir=interp1(1:1:length(phi),phi,linspace(1,length(phi),nk));
        
        plot(Xr(1:ck),Yr(1:ck),SensorModel.Move.trajcol{ss})
        plot_triangle(Xr(ck),Yr(ck),6,1.5,phir(ck)-pi/2,SensorModel.Move.trajcol{ss},'no')
        plot_ellipse(Xr(ck),Yr(ck),SensorModel.Move.FOV{ss}(2),SensorModel.Move.FOV{ss}(2),0,SensorModel.Move.trajcol{ss},'fill')
        
        if curk==1
        indp=find(Trajectories{ss}(:,4)==1);
        else
        indp=find(Trajectories{ss}(:,4)==curk-1);    
        end
        
        plot(Trajectories{ss}(1:indp(end),1),Trajectories{ss}(1:indp(end),2),SensorModel.Move.trajcol{ss})
%         pos=SensorModel.Move.Grid2Abs(SensorModel.Move.Gridpos{ss}(curk,1),SensorModel.Move.Gridpos{ss}(curk,2),SensorModel);
%         dir=SensorModel.Move.Gridpos{ss}(curk,3);
%         if dir==1
%             th=-pi/2;
%         elseif dir==2
%             th=0;
%         elseif dir==3
%             th=pi/2;
%         elseif dir==4
%             th=pi;
%         else
%             th=Nan;
%         end
%         plot_triangle(pos(1),pos(2),6,1.5,th,SensorModel.Move.trajcol{ss},'no')
%         plot_ellipse(pos(1),pos(2),SensorModel.Move.FOV{ss}(2),SensorModel.Move.FOV{ss}(2),0,SensorModel.Move.trajcol{ss},'fill')
    end
% end
axis([0 TargetPos.Traffic_xlim 0 TargetPos.Traffic_ylim])

