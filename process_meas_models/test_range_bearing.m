
meth='ut';

    switch lower(meth)
        case 'ckf'
            qd_pts=@cubature_KF_points;
        case 'ut'
            qd_pts=@(m,P)UT_sigmapoints(m,P,2);
        case 'cut4'
            qd_pts=@conjugate_dir_gausspts;
        case 'cut6'
            qd_pts=@conjugate_dir_gausspts_till_6moment_scheme2;
        case 'cut8'
            qd_pts=@conjugate_dir_gausspts_till_8moment;
        case 'gh'
            qd_pts=@(m,P)GH_pts(m,P,para);
        otherwise
            error('smthg is wrong: DONT ask me what')
    end
    
    f=@(x)[sqrt(x(1)^2+x(2)^2),atan2(x(2),x(1))];
    
    m=[5;5];
    P=[2,-1;-1,5];
    eig(P)

%%    
    [X,w]=GH_pts(m,P,2);
    Xgh2=zeros(size(X));
    wgh2=w;
    for i=1:1:length(w)
        Xgh2(i,:)=f(X(i,:));
    end
    [y1,M1gh2]=Cal_moments_samples(Xgh2,wgh2,1,'raw');
    [y2,M2gh2]=Cal_moments_samples(Xgh2,wgh2,2,'raw');
    [y3,M3gh2]=Cal_moments_samples(Xgh2,wgh2,3,'raw');
    [y4,M4gh2]=Cal_moments_samples(Xgh2,wgh2,4,'raw');
    
    
%%
    [X,w]=GH_pts(m,P,3);
    Xgh3=zeros(size(X));
    wgh3=w;
    for i=1:1:length(w)
        Xgh3(i,:)=f(X(i,:));
    end
    [y1,M1gh3]=Cal_moments_samples(Xgh3,wgh3,1,'raw');
    [y2,M2gh3]=Cal_moments_samples(Xgh3,wgh3,2,'raw');
    [y3,M3gh3]=Cal_moments_samples(Xgh3,wgh3,3,'raw');
    [y4,M4gh3]=Cal_moments_samples(Xgh3,wgh3,4,'raw');


%%
    [X,w]=GH_pts(m,P,4);
    Xgh4=zeros(size(X));
    wgh4=w;
    for i=1:1:length(w)
        Xgh4(i,:)=f(X(i,:));
    end
    [y1,M1gh4]=Cal_moments_samples(Xgh4,wgh4,1,'raw');
    [y2,M2gh4]=Cal_moments_samples(Xgh4,wgh4,2,'raw');
    [y3,M3gh4]=Cal_moments_samples(Xgh4,wgh4,3,'raw');
    [y4,M4gh4]=Cal_moments_samples(Xgh4,wgh4,4,'raw');

%%
    [X,w]=GH_pts(m,P,5);
    Xgh5=zeros(size(X));
    wgh5=w;
    for i=1:1:length(w)
        Xgh5(i,:)=f(X(i,:));
    end
    [y1,M1gh5]=Cal_moments_samples(Xgh5,wgh5,1,'raw');
    [y2,M2gh5]=Cal_moments_samples(Xgh5,wgh5,2,'raw');
    [y3,M3gh5]=Cal_moments_samples(Xgh5,wgh5,3,'raw');
    [y4,M4gh5]=Cal_moments_samples(Xgh5,wgh5,4,'raw');

%%
    [X,w]=GH_pts(m,P,6);
    Xgh6=zeros(size(X));
    wgh6=w;
    for i=1:1:length(w)
        Xgh6(i,:)=f(X(i,:));
    end
    [y1,M1gh6]=Cal_moments_samples(Xgh6,wgh6,1,'raw');
    [y2,M2gh6]=Cal_moments_samples(Xgh6,wgh6,2,'raw');
    [y3,M3gh6]=Cal_moments_samples(Xgh6,wgh6,3,'raw');
    [y4,M4gh6]=Cal_moments_samples(Xgh6,wgh6,4,'raw');
    
%%
    [X,w]=GH_pts(m,P,7);
    Xgh7=zeros(size(X));
    wgh7=w;
    for i=1:1:length(w)
        Xgh7(i,:)=f(X(i,:));
    end
    [y1,M1gh7]=Cal_moments_samples(Xgh7,wgh7,1,'raw');
    [y2,M2gh7]=Cal_moments_samples(Xgh7,wgh7,2,'raw');
    [y3,M3gh7]=Cal_moments_samples(Xgh7,wgh7,3,'raw');
    [y4,M4gh7]=Cal_moments_samples(Xgh7,wgh7,4,'raw');
    
    
%%
    [X,w]= UT_sigmapoints(m,P,2);
    Xut=zeros(size(X));
    wut=w;
    for i=1:1:length(w)
        Xut(i,:)=f(X(i,:));
    end
    [y1,M1ut]=Cal_moments_samples(Xut,wut,1,'raw');
    [y2,M2ut]=Cal_moments_samples(Xut,wut,2,'raw');
    [y3,M3ut]=Cal_moments_samples(Xut,wut,3,'raw');
    [y4,M4ut]=Cal_moments_samples(Xut,wut,4,'raw');
    
    
%%
    [X,w]=conjugate_dir_gausspts(m,P);
    Xcut4=zeros(size(X));
    wcut4=w;
    for i=1:1:length(w)
        Xcut4(i,:)=f(X(i,:));
    end
    [y1,M1cut4]=Cal_moments_samples(Xcut4,wcut4,1,'raw');
    [y2,M2cut4]=Cal_moments_samples(Xcut4,wcut4,2,'raw');
    [y3,M3cut4]=Cal_moments_samples(Xcut4,wcut4,3,'raw');
    [y4,M4cut4]=Cal_moments_samples(Xcut4,wcut4,4,'raw');


%%
    [X,w]=conjugate_dir_gausspts_till_6moment_scheme2(m,P);
    Xcut6=zeros(size(X));
    wcut6=w;
    for i=1:1:length(w)
        Xcut6(i,:)=f(X(i,:));
    end
    [y1,M1cut6]=Cal_moments_samples(Xcut6,wcut6,1,'raw');
    [y2,M2cut6]=Cal_moments_samples(Xcut6,wcut6,2,'raw');
    [y3,M3cut6]=Cal_moments_samples(Xcut6,wcut6,3,'raw');
    [y4,M4cut6]=Cal_moments_samples(Xcut6,wcut6,4,'raw');


%%
    [X,w]=conjugate_dir_gausspts_till_8moment(m,P);
    Xcut8=zeros(size(X));
    wcut8=w;
    for i=1:1:length(w)
        Xcut8(i,:)=f(X(i,:));
    end
    [y1,M1cut8]=Cal_moments_samples(Xcut8,wcut8,1,'raw');
    [y2,M2cut8]=Cal_moments_samples(Xcut8,wcut8,2,'raw');
    [y3,M3cut8]=Cal_moments_samples(Xcut8,wcut8,3,'raw');
    [y4,M4cut8]=Cal_moments_samples(Xcut8,wcut8,4,'raw');        
    
    
    
%% Comparing


Xg=[M1gh2,M1gh3,M1gh4,M1gh5,M1gh6,M1gh7]
Xc=[M1ut,M1cut4,M1cut6,M1cut8]
n=1;
plot(1:size(Xg,2),Xg(n,:),'s-',1:size(Xc,2),Xc(n,:),'o-','MarkerSize',8)
legend('gh2-7','ut-cut4-cut6-cut8')

Xg=[M2gh2,M2gh3,M2gh4,M2gh5,M2gh6,M2gh7]
Xc=[M2ut,M2cut4,M2cut6,M2cut8]
n=2;
plot(1:size(Xg,2),Xg(n,:),'s-',1:size(Xc,2),Xc(n,:),'o-','MarkerSize',8)
legend('gh2-7','ut-cut4-cut6-cut8')

Xg=[M3gh2,M3gh3,M3gh4,M3gh5,M3gh6,M3gh7]    
Xc=[M3ut,M3cut4,M3cut6,M3cut8]
n=4;
plot(1:size(Xg,2),Xg(n,:),'s-',1:size(Xc,2),Xc(n,:),'o-','MarkerSize',8)
legend('gh2-7','ut-cut4-cut6-cut8')


Xg=[M4gh2,M4gh3,M4gh4,M4gh5,M4gh6,M4gh7] 


%%



X=[M2ut,M2cut4,M2cut6,M2cut8]
plot(1:size(X,2),X(4,:))