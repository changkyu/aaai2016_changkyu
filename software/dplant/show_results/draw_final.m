dir_summary = '~/tmp2/results_final2/summary';

fontsz_title = 14;
fontsz_x     = 12;
fontsz_y     = 12;

marks = {'s','o','x'};
colors = ...
[ ...
    0.9290    0.6940    0.1250; % ADMM
    0    0.4470    0.7410;% ADMM (AP)
    0.8500    0.3250    0.0980;% ADMM (NAP)
];
    


idx_ex = 4;
range_models = [ 4 6 5];

experiments_desc = load_experiments_desc;
expr_desc = experiments_desc{idx_ex};
name_expr = expr_desc.name;
load(fullfile(dir_summary,sprintf('%s_summary.mat',name_expr)));

idx_ETA = 1;
idx_Network =1;
idx_NV=1;
idx_run = 1;

for idx_ex_param = [1 3];
    idx_ex_param2 = 1;
    
iter = zeros(length(range_models),1);
ssa  =  cell(length(range_models),1);
for i=1:length(range_models)                        
    iter(i) = squeeze(results{range_models(i),idx_ex_param,idx_ex_param2}.iter(idx_ETA, :, idx_Network, :));
    ssa{i}  = results{range_models(i),idx_ex_param,idx_ex_param2}.ssa{idx_ETA, idx_NV, idx_Network, idx_run};
                    
end         

figure();
hold on;

for i=1:length(range_models)
    plot( 1:iter(i), ...
            180/pi*ssa{i}(1:iter(i)), ...
            [marks{i} '-'], 'color', colors(i,:), 'LineWidth', 2 );
end

set(haxes,'YScale','log');

xlim([0, max(iter_median(:,idx_NV))]);
ylim([0, inf]);

name_target = 'Subspace Error';
name_xlabel = 'iteration';
name_ylabel = 'subspace error';
str_title = sprintf('%s %s \n(ETA:%d Nodes:%d %s)',names_expr{idx_ex},name_target,ETA,NV,names_Networks{idx_Network});
title( str_title, 'FontSize', fontsz_title);
xlabel(name_xlabel,'FontSize',fontsz_x);
ylabel(name_ylabel,'FontSize',fontsz_y);
legend( {names_Models{range_models}}, 'Location','SW');

str_title = sprintf('%s %s (ETA:%d Nodes:%d %s)',names_expr{idx_ex},name_target,ETA,NV,names_Networks{idx_Network});
file_img = strrep(str_title, ' ', '_');
file_img = strrep(file_img, ':', '');
file_img = strrep(file_img, '(', '');                    
file_img = strrep(file_img, ')', '');                                        
saveas(gcf, sprintf('../../img/%s',file_img), 'epsc');
saveas(gcf, sprintf('../../img/%s',file_img), 'png');

end