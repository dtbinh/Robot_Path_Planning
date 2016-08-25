close all

prior_GMM.mu{2}=[70,200,0,2,-0.03;75,210,0,2,-0.03];
prior_GMM.P{2}=[reshape(diag([1,1,0.1,0.1,(0.0001)^2]),1,25);reshape(diag([1,1,0.1,0.1,(0.0001)^2]),1,25)];
prior_GMM.w{2}=[0.5;0.5];

prrGmm.mu=prior_GMM.mu{2};
prrGmm.P=prior_GMM.P{2};
prrGmm.w=prior_GMM.w{2};
obj=GMM2obj(prrGmm);
 YY= random(obj);
 targ_pos(2,:)=YY(1:2);
 mov2=YY;
 figure(1)

 for i=1:1:20
     mov2=vertcat(mov2,targ_dyn{2}(mov2(end,:)',para)');
    prGMM1.mu=prior_GMM.mu{2};prGMM1.P=prior_GMM.P{2};prGMM1.w=prior_GMM.w{2};   
    prGMM1=prior_prop_GMM(prGMM1,targ_dyn{2},para{2},Q{2});
    prior_GMM.mu{2}=prGMM1.mu;prior_GMM.P{2}=prGMM1.P;prior_GMM.w{2}=prGMM1.w;
    plot_UAV_GMM(c0,d,xlim,prior_GMM,targ_pos,sts,r)
     hold on
    plot(mov2(:,1),mov2(:,2),'ko-')
     axis([0,400,0,400])
    pause(1)
    hold off
 end

