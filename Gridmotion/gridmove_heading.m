% grid movement of UAV way points.
function ns=gridmove_heading(X,xc,yc,Nn,ps,u)

%present state [ci,ha]
% control u=-1:-turn right, 0:-straight, 1:-turn left
ci=ps(1);
ha=ps(2);

if u==0 % go straight
    ha=ha+u;
elseif u==-1 %turn right
    ha=ha+u;
    if ha==0
        ha=4;
    end
elseif u==1 % turn left
    ha=ha+u;
    if ha==5
        ha=1;
    end    
end

if ha==1  %heading angle towards east
    if X(ci,1)==xc
       ns=[NaN,NaN];
        return
    end
    ni=ci+1;
elseif ha==3 %heading angle is towards west
        if X(ci,1)==0
        ns=[NaN,NaN];
        return
        end
    ni=ci-1;
elseif ha==2 %heading north
        if X(ci,2)==yc
        ns=[NaN,NaN];
        return
        end
    ni=ci+Nn;
elseif ha==4  %heading south
        if X(ci,2)==0
        ns=[NaN,NaN];
        return
        end
    ni=ci-Nn;
end

ns=[ni,ha];