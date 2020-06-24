colors = ...
[ ...
    0.9290    0.6940    0.1250; % ADMM
    0         0.4470    0.7410; % ADMM (AP)
    0.8500    0.3250    0.0980; % ADMM (NAP)
];

colors_nodes = ...
[...
    1         0         0;
    0.9290    0.6940    0.1250;
    0.4660    0.6740    0.1880;
    0         0.4470    0.7410;
    0.4940    0.1840    0.5560;
];

marks = {'s', 'o','x'};

idx_run = 1;
ETA = 10;
NV = 5;

name_expr = 'Caltech';
param_expr = 4;
idx_obj = param_expr;
name_Network = 'Complete';

filepath_save_mat = sprintf('draw_caltech_obj%d.mat',idx_obj);
load(filepath_save_mat);
im = imread(sprintf('%d.JPG',idx_obj));

name_target = 'DPPCA-NAP';
name_xlabel = 'iteration';
name_ylabel = 'penalty constraint (eta)';

close gcf;
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/2-640 scrsz(4)/2-360 1280 720]);
winsize = get(gcf, 'Position');

writerObj = VideoWriter(sprintf('~/tmp2/caltech_obj%d.avi',idx_obj));
writerObj.FrameRate = 6; % second pause
writerObj.Quality = 100;
open(writerObj);

set(gca, 'NextPlot', 'replacechildren');
set(gcf, 'Renderer', 'zbuffer');

h_axes = axes;
pt = get(h_axes,'Position');

str_title = sprintf('%s %s \n(ETA:%d Obj:%d %s)',name_expr,name_target,ETA,param_expr,name_Network);                    
title(str_title);

eITER = max(cm1{idx_obj,1,idx_run}.eITER,cm1{idx_obj,1,idx_run}.eITER);
max_ssa = 180/pi*max(max(cm1{idx_obj,1,idx_run}.ssaArray), max(cm3{idx_obj,1,idx_run}.ssaArray));

min_iter = min(cm1{idx_obj,1,idx_run}.eITER, cm3{idx_obj,1,idx_run}.eITER);

for idx_iter = 1:(min_iter-1)
% idx_iter_cm1 = min(cm1{idx_obj,1,1}.eITER,idx_iter);
% idx_iter_cm3 = min(cm3{idx_obj,1,1}.eITER,idx_iter);

idx_iter_cm1 = min(cm1{idx_obj,1,idx_run}.eITER,min_iter) - 1;
idx_iter_cm3 = min(cm3{idx_obj,1,idx_run}.eITER,min_iter) - 1;

if( sum(idx_iter == [1 2 3 4 5 6 12 35 min_iter-1]) > 0) %1 12 21 35 cm3{idx_obj,1,1}.eITER
   range_x = [60:10:180 -180:10:60];
else 
    range_x = 60;
end

for x=range_x

clf;

h_ssa = axes('Position',[pt(1) pt(2) pt(3)/2 pt(4)/2]);
hold on;

% semilogy(h_ssa, 1:min(cm1{idx_obj,1,1}.eITER,idx_frame), ...
%   180/pi*cm1{idx_obj,1,1}.ssaArray(1:min(cm1{idx_obj,1,1}.eITER,idx_frame)), ...
%   [marks{1} '-'], 'color', colors(1,:), 'LineWidth', 2 );
semilogy(h_ssa, 1:idx_iter_cm1, ...
  180/pi*cm1{idx_obj,1,idx_run}.ssaArray(1:idx_iter_cm1), ...
  [marks{1} '-'], 'color', colors(1,:), 'LineWidth', 2 );

semilogy(h_ssa, 1:idx_iter_cm3, ...
  180/pi*cm3{idx_obj,1,idx_run}.ssaArray(1:idx_iter_cm3), ...
  [marks{3} '-'], 'color', colors(3,:), 'LineWidth', 2 );

semilogy(h_ssa, [idx_iter idx_iter], [0 max_ssa], '-k');

xlim(h_ssa, [0, min_iter]);
ylim(h_ssa, [0, max_ssa]);

xlabel('iteration','FontSize',12);
ylabel('subspace error','FontSize',12);

legend( {'DPPCA','DPPCA-NAP'}, 'Location','NE');

h_img = axes('Position',[pt(1) pt(2) pt(3)/6 pt(4)/6]);
imshow(im);

for a=1:NV
%h_eta = subplot(12,2,12+a*2);
h_eta = axes('Position',[pt(1)+pt(3)/2+pt(3)/10 pt(2)+(NV-a)*pt(4)/10 pt(3)/2-pt(3)/10 pt(4)/10-pt(4)/40]);
xlim(h_eta,[0.5, min_iter-0.5]);
ylim(h_eta,[ETA*0.5-1 max(max(max(cm3{idx_obj,1,idx_run}.ETA_ij_history))) + 1]);

hold on;
for b=1:NV
    if a==b 
        continue;
    end                  
    if cm3{idx_obj,1,idx_run}.ETA_ij_history(1:idx_iter_cm3,a,b)==0;
        continue;
    end
     [XX, YY] = stairs( (1:idx_iter_cm3) - 0.5, ...
                        cm3{idx_obj,1,idx_run}.ETA_ij_history((1:idx_iter_cm3)+1,a,b));   
     plot(XX, YY,            'color', colors_nodes(b,:), 'LineWidth', 1.5 );

end

text(idx_iter_cm3+2,8,['(' num2str(a) ')'],'color',colors_nodes(a,:), 'FontSize',12,'FontWeight','bold');
                    

if a == ceil(NV/2)
    ylabel('eta','FontSize',12);    
end

if a == NV
    xlabel('iteration','FontSize',12);
end

if a < NV
    set(gca,'XTick',[]);
else

end              
end

ah = axes('visible','off','Position',[pt(1)+pt(3)/2+pt(3)/10 pt(2) pt(3)/2-pt(3)/10 pt(4)/2]);
xlim([0 min_iter]);
ylim([0 1]);
hold on;
line('XData',[idx_iter idx_iter] - 0.5,'YData',[0 1], 'parent',ah,'linewidth',1,'LineStyle','-','color',[0 0 0]);

h_network = axes('visible','off','Position',[pt(1) + pt(3)*4/10, pt(2) + pt(4)*1/2 + pt(4)/20, pt(3)/5, pt(3)/4]);%pt(4)/4 - pt(4)/20  ]);
hold on;
show_topology( NV, idx_iter+1, ETA, cm3{idx_obj,1,idx_run}.ETA_ij_history, colors_nodes );

for idx_node = 1:5
    
if( idx_node == 2)
    h_recon = axes('Position',[pt(1) + pt(3)*1/10, pt(2) + pt(4)*1/2 + pt(4)/10 + pt(4)/4, pt(3)/5, pt(3)/5] );%pt(4)/4 - pt(4)/10  ]);
elseif( idx_node == 1)
    h_recon = axes('Position',[pt(1) + pt(3)*4/10, pt(2) + pt(4)*1/2 + pt(4)/10 + pt(4)/4, pt(3)/5, pt(3)/5] );%pt(4)/4 - pt(4)/10 ]);
elseif( idx_node == 5)
    h_recon = axes('Position',[pt(1) + pt(3)*7/10, pt(2) + pt(4)*1/2 + pt(4)/10 + pt(4)/4, pt(3)/5, pt(3)/5] );%pt(4)/4 - pt(4)/10 ]);
elseif( idx_node == 3)
    h_recon = axes('Position',[pt(1) + pt(3)*1/10, pt(2) + pt(4)*1/2 + pt(4)/10,           pt(3)/5, pt(3)/5] );%pt(4)/4 - pt(4)/10 ]);
elseif( idx_node == 4)
    h_recon = axes('Position',[pt(1) + pt(3)*7/10, pt(2) + pt(4)*1/2 + pt(4)/10,           pt(3)/5, pt(3)/5] );%pt(4)/4 - pt(4)/10 ]);
end
 
% xyz_min = min(min(cm3{idx_obj,1,1}.W_history(:,idx_node,:,:),[],1),[],3);
% xyz_max = max(max(cm3{idx_obj,1,1}.W_history(:,idx_node,:,:),[],1),[],3);
% 
% xyz_min = squeeze(xyz_min);
% xyz_max = squeeze(xyz_max);
% axis([xyz_min(3),xyz_max(3), ...
%     xyz_min(1),xyz_max(1), ...
%     xyz_min(2),xyz_max(2) ...
% ]);

%axis equal;
axis([-250 250 -250 250 -250 250]);


% hold on;
plot3(  squeeze(cm3{idx_obj,1,idx_run}.W_history(idx_iter,idx_node,:,3)), ...
        squeeze(cm3{idx_obj,1,idx_run}.W_history(idx_iter,idx_node,:,1)),...
        -squeeze(cm3{idx_obj,1,idx_run}.W_history(idx_iter,idx_node,:,2)), '.', 'color',colors_nodes(idx_node,:));

% plot3(  squeeze(GT{idx_obj,1,1}(:,3)), ...
%         squeeze(GT{idx_obj,1,1}(:,1)),...
%         -squeeze(GT{idx_obj,1,1}(:,2)), '.', 'color',colors_nodes(idx_node,:));

% plot3(  squeeze(cm1{idx_obj,1,1}.W{1}(:,3)), ...
%         squeeze(cm1{idx_obj,1,1}.W{1}(:,1)),...
%         -squeeze(cm1{idx_obj,1,1}.W{1}(:,2)), '.', 'color',colors_nodes(idx_node,:));

set(gca,'XTick',[]); set(gca,'XColor','w')
set(gca,'YTick',[]); set(gca,'YColor','w')
set(gca,'ZTick',[]); set(gca,'ZColor','w')
axis off;
view([x,15]);    

end

drawnow;
writeVideo(writerObj, getframe(gcf));

end

end
% close video
close(writerObj);
