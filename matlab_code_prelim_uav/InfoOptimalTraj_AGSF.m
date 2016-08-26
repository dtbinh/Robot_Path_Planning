function [mu,Hnm,PDFk,I]=InfoOptimalTraj_AGSF(X,Nn,xc,yc,Nk,Nt,GMM,model,sens)
 Nt=Nt-Nk; %remaing steps

si=model.si;
A=model.A;
C=model.C;

model.frequpdt=0;
model.timeupdt=1;

% GMM.wtup=0;


%% evolve the pdf 
PDFk=cell(1,Nt);
PDFk{1}=GMM;
Hnm(1)=entropy_gmm(GMM);
for k=2:1:Nt
 [GMM,model,sens]=GMM_Filt_MOVETarg(GMM,model,sens,1);
PDFk{k}=GMM; 
 Hnm(k)=entropy_gmm(GMM);
%  plot_pdf(GMM,xc,yc,model)
%  pause(1)
end

%% find the optimal Mutual Information
I=zeros(Nt,size(X,1),3); %first col is information and next to cols are [phi,alpha]
tic
for k=1:1:Nt
parfor i=1:size(X,1)
    disp(strcat('Time step = ',num2str(k),' and node =',num2str(i)))
    I(k,i,:)=max_Info(X(i,:),PDFk{k},model);
end
end
toc
%% Finding the optimal traj using Dynamic Programming
disp('Now computing the optimal trajectory')
J=zeros(Nt+1,Nn^2,4); % optima cost to go function for all nodes
mu=zeros(Nt+1,Nn^2,4); % optimal policy for all nodes

 for i=1:1:Nn^2
     for j=1:1:4
    J(Nt+1,i,j)=si(X(i,:));
     end
 end
 
 for k=Nt:-1:1
   
    for ci=1:1:Nn^2
        for ch=1:1:4
        JJ=0;
        u=0;
        p=1;
        for j=-1:1:1
            ns=gridmove_heading(X,xc,yc,Nn,[ci,ch],j);
            if isnan(ns(1))~=1
                u(p)=j;
                JJ(p)=-100*I(k,ci,1)+C*u(p)^2+J(k+1,ns(1),ns(2));%+(X(ci,:)-GMM.mean')*A*(X(ci,:)-GMM.mean')'
                p=p+1;
            end
        end
                [J(k,ci,ch),ind]=min(JJ);
                mu(k,ci,ch)=u(min(ind));
        end
    end
end