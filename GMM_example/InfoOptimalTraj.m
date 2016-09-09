function mu=InfoOptimalTraj(X,Nn,xc,yc,Nk,Nt,GMM,model)
 Nt=Nt-Nk; %remaing steps

si=model.si;
A=model.A;
C=model.C;
%% find the optimal Mutual Information
I=zeros(size(X,1),3); %first col is information and next to cols are [phi,alpha]
parfor i=1:1:size(X,1)
    i
     tic
    I(i,:)=max_Info(X(i,:),GMM,model);
     toc
end
%% Finding the optimal traj using Dynamic Programming
disp('Now computing the optimal trajectory')
J=zeros(Nt,Nn^2,4); % optima cost to go function for all nodes
 mu=zeros(Nt,Nn^2,4); % optimal policy for all nodes

 for i=1:1:Nn^2
     for j=1:1:4
    J(Nt,i,j)=si(X(i,:));
     end
 end
 
 for k=Nt-1:-1:1
   
    for ci=1:1:Nn^2
        for ch=1:1:4
        JJ=0;
        u=0;
        p=1;
        for j=-1:1:1
            ns=gridmove_heading(X,xc,yc,Nn,[ci,ch],j);
            if isnan(ns(1))~=1
                u(p)=j;
                JJ(p)=-1*I(ci,1)+C*u(p)^2+J(k+1,ns(1),ns(2));%+(X(ci,:)-GMM.mean')*A*(X(ci,:)-GMM.mean')'
                p=p+1;
            end
        end
                [J(k,ci,ch),ind]=min(JJ);
                mu(k,ci,ch)=u(min(ind));
        end
    end
end