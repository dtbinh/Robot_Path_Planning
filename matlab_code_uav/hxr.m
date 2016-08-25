function [y,G,tp]= hxr(xtarg,xsenspos,FOV)
%dirn is the direction of looking/view/FOV
alpha=FOV(1);
Rmax=FOV(2);
dirn=FOV(3);

r=norm([xtarg(1),xtarg(2)]-[xsenspos(1),xsenspos(2)]);
th=atan2(xtarg(2)-xsenspos(2),xtarg(1)-xsenspos(1));
%alpha is (-pi,pi)
%dirn is (-pi,pi)
%th is (-pi,pi)

diff=dirn-th;
if diff>pi
    diff=diff-2*pi;
end
if diff<-pi
    diff=diff+2*pi;
end
   y=[r];
   tp=[0];% 0 =='r' , 1 == 'th'
if r>Rmax || abs(diff)>alpha
G=NaN;
else
G=1; 
end




