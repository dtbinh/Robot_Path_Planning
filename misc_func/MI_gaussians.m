function MI=MI_gaussians(priorpdf,Pcc,m)
zpdf.mu=zpdf.mu(:);
 [z,w]=GH_points(3,zpdf.mu,zpdf.P);
 KL=0;
 n=length(priorpdf.mu);
 for i=1:1:length(w)
     
     %getting posterior pdf
     K=Pcc/zpdf.P;
     postpdf.mu=mk+K*(z(i,:)'-zpdf.mu);
     postpdf.P=Pk-K*Pz*K';
    KL=KL+w(i)*0.5*(trace(inv(priorpdf.P)*postpdf.P)+(priorpdf.mu-postpdf.mu)'*inv(priorpdf.P)(priorpdf.mu-postpdf.mu)-n-log(det(postpdf.P)/det(priorpdf.P))); 
 end







end