function plot_sens_view(xc,bnd,phi,l,col)
hold on
plot([xc(1),xc(1)+l*cos(phi)],[xc(2),xc(2)+l*sin(phi)],'k--','Linewidth',2)
fill([xc(1),xc(1)+l*cos(phi+bnd),xc(1)+l*cos(phi-bnd),xc(1)],[xc(2),xc(2)+l*sin(phi+bnd),xc(2)+l*sin(phi-bnd),xc(2)],col,'FaceAlpha',0.1)
hold off
end


