%%  Comparing UKF and CUT8


FixedRadFiltering('ut')
FixedRadFiltering('cut8')

%% without correction
load('Sim_ut_FixRadar_Comp_1')
TargetPos.Dyn
err=cell(TargetPos.Ntar,1);
for ntar=1:1:TargetPos.Ntar
    err{ntar}=TargetPos.xf{ntar}(:,[1,2])-TargetPos.truth{ntar}; 
    err{ntar}=sqrt( sum(err{ntar}.^2,2)    );
end
errUKF=err;

load('Sim_cut8_FixRadar_Comp_1')
err=cell(TargetPos.Ntar,1);
for ntar=1:1:TargetPos.Ntar
    err{ntar}=TargetPos.xf{ntar}(:,[1,2])-TargetPos.truth{ntar}; 
    err{ntar}=sqrt( sum(err{ntar}.^2,2)    );
end
errCUT8=err;

ntar=1;

figure
 plot(Time.Tvec,errUKF{ntar},Time.Tvec,errCUT8{ntar})
 legend('UKF','CUT8')

 %% with constant correction
load('Sim_ut_FixRadar_Comp_2')
TargetPos.Dyn
err=cell(TargetPos.Ntar,1);
for ntar=1:1:TargetPos.Ntar
    err{ntar}=TargetPos.xf{ntar}(:,[1,2])-TargetPos.truth{ntar}; 
    err{ntar}=sqrt( sum(err{ntar}.^2,2)    );
end
errUKF=err;

load('Sim_cut8_FixRadar_Comp_2')
err=cell(TargetPos.Ntar,1);
for ntar=1:1:TargetPos.Ntar
    err{ntar}=TargetPos.xf{ntar}(:,[1,2])-TargetPos.truth{ntar}; 
    err{ntar}=sqrt( sum(err{ntar}.^2,2)    );
end
errCUT8=err;

ntar=2;
figure
 plot(Time.Tvec,errUKF{ntar},Time.Tvec,errCUT8{ntar})
 legend('UKF','CUT8') 
 
  %% with adaptive correction
  nn=3;
load(strcat('Sim_ut_FixRadar_Comp_',num2str(nn)))
TargetPos.Dyn
err=cell(TargetPos.Ntar,1);
for ntar=1:1:TargetPos.Ntar
    err{ntar}=TargetPos.xf{ntar}(:,[1,2])-TargetPos.truth{ntar}; 
    err{ntar}=sqrt( sum(err{ntar}.^2,2)    );
end
errUKF=err;

load(strcat('Sim_cut8_FixRadar_Comp_',num2str(nn)))
err=cell(TargetPos.Ntar,1);
for ntar=1:1:TargetPos.Ntar
    err{ntar}=TargetPos.xf{ntar}(:,[1,2])-TargetPos.truth{ntar}; 
    err{ntar}=sqrt( sum(err{ntar}.^2,2)    );
end
errCUT8=err;

ntar=2;
figure
 plot(Time.Tvec,errUKF{ntar},Time.Tvec,errCUT8{ntar})
 legend('UKF','CUT8') 