function xk1=diffdyn(xk,m0,cf)
xk=xk(:);
m0=m0(:);
dr=xk-m0;
if norm(dr)==0
    dr=0;
else
    dr=dr/norm(dr);
end
xk1=xk+dr*cf+0.2;