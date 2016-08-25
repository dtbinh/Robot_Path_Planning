function [Y,X,XY,XYindex,XYADJindex,XYgridpos]=SetupGrid(Xlim,dx)
            
            
            [X,Y]=meshgrid(Xlim(1):dx:Xlim(2));
%             XY=[Y(:),X(:)];
            XY=zeros(size(X,1)*size(X,2),2);
            XYindex=zeros(size(X,1)*size(X,2),1);
            XYgridpos=zeros(size(X,1)*size(X,2),2);
            k=1;
            for i=1:1:size(X,1)
                for j=1:1:size(X,2)
                    XY(k,:)=[Y(i,j),X(i,j)];
                    XYindex(k)=k;
                    XYgridpos(k,:)=[i,j];
                    k=k+1;
                end
            end
%             XYindex=[1:1:size(XY,1)]';
            % [XY,XYindex]
            ng=length(XYindex);
            XYADJindex=cell(ng,1);
            
            for i=1:1:length(XYindex)
                xy=XY(i,:);
                p=[];
                
                r=xy+[dx,0];
                if r(1)>Xlim(2) || r(2) >Xlim(2) || r(1)<Xlim(1) || r(2)<Xlim(1)
                    
                else
                    [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
                    p=horzcat(p,ind(1));
                end
                
                r=xy+[0,dx];
                if r(1)>Xlim(2) || r(2) >Xlim(2) || r(1)<Xlim(1) || r(2)<Xlim(1)
                    
                else
                    [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
                    p=horzcat(p,ind(1));
                end
                
                r=xy+[-dx,0];
                if r(1)>Xlim(2) || r(2) >Xlim(2) || r(1)<Xlim(1) || r(2)<Xlim(1)
                    
                else
                    [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
                    p=horzcat(p,ind(1));
                end
                
                r=xy+[0,-dx];
                if r(1)>Xlim(2) || r(2) >Xlim(2) || r(1)<Xlim(1) || r(2)<Xlim(1)
                    
                else
                    [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
                    p=horzcat(p,ind(1));
                end
                
                XYADJindex{i}=sort(p);
            end
        end