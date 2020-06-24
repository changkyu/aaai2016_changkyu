%% Do TV using some ADMM variants
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
reset(s,0);

addpath('altmany-export_fig');
addpath('subtightplot');

method_name = {'1_admm', '2_fadmmwr', '3_admmvp', '4_admmap_2node'};
MUs = [0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01];

%% Prepare data
% gen_org_images();

%% Run
for idm = 2 : 2%length(method_name)
    run_tv_one_method(method_name{idm}, MUs);
end

%% Show objective curve
% show_tv_obj_curve(method_name, MUs);
[best_relerr, worst_relerr, total_stats] = show_tv_obj_curve_relerr(method_name, MUs);
% show_tv_obj_curve_subopt(method_name, MUs);

%% Show
% show_tv_org_images();
% for idm = 1 : length(method_name)
%     show_tv_one_method(method_name{idm}, MUs, best_relerr, worst_relerr, idm);
% end
