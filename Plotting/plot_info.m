function plot_info(Nn,X,II)
x=zeros(Nn,Nn);
y=zeros(Nn,Nn);
Z=zeros(Nn,Nn);
for j=1:1:Nn
    for i=1:1:Nn
        x(i,j)=X(i+(j-1)*Nn,1);
        y(i,j)=X(i+(j-1)*Nn,2);
        Z(i,j)=II(1,i+(j-1)*Nn,1);     
    end
end
mesh(x,y,Z)