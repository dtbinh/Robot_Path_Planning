function [M,J]=dubin_DP(t0,dt,tf,c0,Nu,Ns,Nt,d,xlim,ylim,F,targ_dyn,prior_GMM,nz,h,g,R,Q,para,temp,sts,thresh)
% output mu{sens} and J{sens} are the control policy for each sensor
% Ns is no. of sensors and Nt is no. of targets and Nd is the number of
% directions
% sens0 is initial location info about sensors i.e x0 y0 and phi0
%
% calculates path from t0 to tf with intervals dt
% F is the grid dynamics (has to take into account boundaries)
% u is the set of admissible control
% V is the corresponding cost of each control action u (same size)
%  have to be compatible with each other and 
%    are a function of (x,y,phi). i.e. define this function separetly 
        %[xf,yf,phif,Vf]=F(x0,y0,phi0,u)

% prior_GMM is a cell of Gaussian mixes
% targ_models{targ_no} contains all the dyn for each target as a cell
  T=t0:dt:tf;      
[Xgrid,Ygrid]=meshgrid(0:d:xlim,0:d:ylim); % this grid is for sensor motion

% generate info map over time
% the information is in I{t}(xi,yi) i.e. a cell for each time step

[s1,s2]=meshgrid(0:d:xlim);
II=cell(length(T),1);
for ji=1:1:length(T)
II{ji,1}=zeros(size(s2));
end
tic
for k=1:1:length(T)
%     if strcmp('stationary',simple_targ)==1 && k==2
%         break;
%     end
%         
for i=1:1:size(s1,1)
    SS=zeros(1,size(s1,2));
    parfor j=1:size(s1,2)
        m=[];
        
        for nt=1:1:Nt
            prGMM=struct('mu',0,'P',0,'w',0);
%             prGMM.mu=0;
%             prGMM.P=0;
%             prGMM.w=0;
        prGMM.mu=prior_GMM.mu{nt};
        prGMM.P=prior_GMM.P{nt};
        prGMM.w=prior_GMM.w{nt};
        
        if sts~=0
        prGMM=marg_GMM(prGMM,sts);
        end
        [mmu,PPu]=GMM_mean_cov(prGMM);
        if max(eig(sqrtm(PPu)))<=thresh
%            disp('target information ignored')
            continue 
        end
%         prGMM
%         keyboard
        
        measGMM=meas_GMM(prGMM,nz{nt},h{nt},g{nt},R{nt},[s1(i,j),s2(i,j)]);
        m=horzcat(m,Expt_GMM(measGMM,@(z)dist_eval_meas(prGMM,z,h{nt},g{nt},R{nt},[s1(i,j),s2(i,j)],1)));
        end
        SS(j)=sum(m);
        
    end
    II{k}(i,:)=SS;
end
% figure
% mesh(s1,s2,II{k})
% pause(0.1)
% keyboard
for nt=1:1:Nt
        prGMM1.mu=prior_GMM.mu{nt};
        prGMM1.P=prior_GMM.P{nt};
        prGMM1.w=prior_GMM.w{nt};
prGMM1=prior_prop_GMM(prGMM1,targ_dyn{nt},para{nt},Q{nt});
        prior_GMM.mu{nt}=prGMM1.mu;
        prior_GMM.P{nt}=prGMM1.P;
        prior_GMM.w{nt}=prGMM1.w;
        
        prior_GMM.truth{nt}=vertcat(prior_GMM.truth{nt},targ_dyn{nt}(prior_GMM.truth{nt}(end,:)',para{nt})');
end
% plot_UAV_GMM2(d,xlim,prior_GMM,sts,60,Nt)
% keyboard
end
% if strcmp('stationary',simple_targ)==1
%         for k=2:1:length(T)
%            II{k}=II{1}; 
%         end
% %         figure(2)
% % mesh(s1,s2,II{1})
% % xlabel('s1')
% % ylabel('s2')
% end

%% Run DP sequentially centralized 
% setup the priority order
Ps=1:1:Ns;
% initialize DP and run
J=cell(Ns,length(T));
M=cell(Ns,length(T));
for i=1:Ns
    for j=1:length(T)
J{i,j}=zeros(size(Xgrid,1),size(Xgrid,2),4);
J{i,j}(:,:,1)=II{end,1}';
J{i,j}(:,:,2)=II{end,1}';
J{i,j}(:,:,3)=II{end,1}';
J{i,j}(:,:,4)=II{end,1}';
M{i,j}=zeros(size(Xgrid,1),size(Xgrid,2),4);
    end
end
Ic=zeros(length(T),1); % info collected by previous sensors
Xpath=cell(1,1);
Xpath{1}=[-1,-1,-1];
A=1e3;% information weighting
for ss=1:1:Ns
    for k=(length(T)-1):-1:1
        for i=1:1:size(Xgrid,1)
            for j=1:1:size(Xgrid,2)
                for dirn=1:1:4
                    JJ=1e15;
                    umin=0;
                    for U=1:1:Nu
                        [Xf,V]=F(i,j,dirn,U,xlim,d,temp);
                        if isnan(Xf)==1
                        continue
                        end
                        Jc=-A*(Ic(k)+II{k}(j,i))+1e2*collision_penalty(Xpath,i,j,size(Xpath,1),k)+1*V/10+J{Ps(ss),k+1}(Xf(1),Xf(2),Xf(3));
                        if Jc<JJ
                        JJ=Jc;
                        umin=U;
                        end
                    end      
                J{Ps(ss),k}(i,j,dirn)=JJ;
                M{Ps(ss),k}(i,j,dirn)=umin;
                end
            end            
        end
     
    end % time step end
    % calculate the information along the optimal path
    
    Xpath{ss}=[0 0 0];
    for k=1:1:(length(t0:dt:tf)-1)
        Ic(k)=Ic(k)+II{k}(c0{ss}(1,2),c0{ss}(1,1));
        Xpath{ss}(k,:)=c0{ss};
        u=M{ss,k}(c0{ss}(1,1),c0{ss}(1,2),c0{ss}(1,3));
        J{ss,k}(c0{ss}(1,1),c0{ss}(1,2),c0{ss}(1,3));
        [Xf,V]=UAV_grid_motion(c0{ss}(1,1),c0{ss}(1,2),c0{ss}(1,3),u,xlim,d,temp);
        c0{ss}(1,1)=Xf(1);c0{ss}(1,2)=Xf(2);c0{ss}(1,3)=Xf(3); 
        
    end
%     keyboard
end % end for each sensor
toc
% disp('123')
%    keyboard
