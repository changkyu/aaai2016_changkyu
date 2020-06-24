function show_topology( NV, iter, ETA, ETAij, color )
%SHOW_TOPOLOGY Summary of this function goes here
%   Detailed explanation goes here
r = 10;
    for i=1:NV    
        ci = cos((i-1)*(2*pi)/NV + pi/2);
        si = sin((i-1)*(2*pi)/NV + pi/2);

        xi = r*ci;
        yi = r*si;

        for j=1:NV
            if( i~=j )
                cj = cos((j-1)*(2*pi)/NV + pi/2);
                sj = sin((j-1)*(2*pi)/NV + pi/2);

                xj = r*cj;
                yj = r*sj;

                linewidth = (ETAij(iter,j,i)-5)/1.5*0.8;
                if( ETAij(iter,j,i) >= ETA )

                        rr = sqrt((xj - xi)^2 + (yj - yi)^2);
                        cij = (yi - yj)/rr;
                        sij = (xj - xi)/rr;
                        xx = rr * cij * 0.01 *linewidth;
                        yy = rr * sij * 0.01 *linewidth;
                    
                    if( i < j )
                        line([(xi) (xj)] + xx,[(yi) (yj)] + yy,'color',color(i,:), 'LineWidth', linewidth, 'LineStyle','-');
                    else
                        line([(xi) (xj)] + xx,[(yi) (yj)] + yy,'color',color(i,:), 'LineWidth', linewidth, 'LineStyle','-');
                    end
                end
            end        
        end

        %plot(xi+cos(angles), yi+sin(angles),color(i,:), 'LineWidth', 2);
    
    end
    
    for i=1:NV    
        ci = cos((i-1)*(2*pi)/NV + pi/2);
        si = sin((i-1)*(2*pi)/NV + pi/2);

        xi = r*ci;
        yi = r*si;

        h = plot(xi,yi,'o');
        set(h(1),'LineWidth',4,'MarkerEdgeColor',color(i,:),'MarkerFaceColor',color(i,:))
        text(xi+5*ci,yi+5*si,num2str(i),'color',color(i,:), 'FontSize',12,'FontWeight','bold');
    
    end
    
    title(['iter = ' num2str(iter)],'FontSize', 12);
    %xlabel(['iter = ' num2str(iter)]);
    set(gca,'Visible','off');
    set(gca,'XTick',[])
    set(gca,'XColor','w')
    set(gca,'YTick',[])
    set(gca,'YColor','w')    
    
    axis([-15 16 -15 16]);

end

