function s=det_leibniz(A,m)
% m is the dimension of A
t0=clock;
r=0;
if m==2
    s=A(1,1)*A(2,2)-A(1,2)*A(2,1);

elseif m==1
    s=A(1,1);
    
elseif m<=6 && m>2  % just use leibniz rule
    
    P = perms(1:m);

    II=eye(m);
    
    s=0;
    for n=1:1:size(P,1)
        p=1;
        for i=1:1:m
            p=p*A(i,P(n,i));
            
        end

%         p=prod(diag(A(1:1:m,P(n,1:1:m))));
        
        
        s=s+det(II(P(n,:),:) )*p;
        
        
    end
    
else
    r=1;
    s=0;
    for j=1:1:m
        B=A;
        B(1,:)=[];
        B(:,j)=[];
       s=s+(-1)^(j+1)*A(1,j)*det_leibniz(B,m-1);
    end
    
end
% 
if r==1
% disp(clock-t0)
end