mov(1:Nt-1) = struct('cdata', [],...
                        'colormap', []);
                    close all
for k=1:1:Nt-1
figure
    plot_pdf(new_PDFk{k},xc,yc,model)
    hold on
    plot(x_mc(1:k,1),x_mc(1:k,2),'k*-')



    for hg=1:1:size(newpath,1)
  
           plot(Xtraj(newpath(hg,2):k,1),Xtraj(newpath(hg,2):k,2),newcol{newpath(hg,1)},'linewidth',2,'MarkerSize',6)

    end

    plot_sens_view(X(sensor_config(k,1),:),sensor_config(k,4),sensor_config(k,3),model.rmax,'b')
    axis([-10,xc+10,-10,yc+10])
    xlabel('x-axis')
    ylabel('y-axis')
     plot_prop_paper
     pause(1)
 saveas(h,strcat('MOVsensSTATtarg_LimView_',num2str(k)),'jpg')
    mov(k) = getframe(gcf);
           pause(1) 
           close
end
% movie(mov)

%  movie2avi(mov,'MOVsensSTATtarg_LimView_no_measupdt','compression','none','fps',1)



figure
plot(1:1:length(Hnm),Hnm,'bo-',1:1:length(Hc),Hc,'rs-','linewidth',2.5,'MarkerSize',6)
xlabel('Time step')
ylabel('Joint Entropy')
legend('No meas','With meas')
   plot_prop_paper