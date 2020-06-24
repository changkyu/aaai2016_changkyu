% Envirnonment Setting

dir_root     = '../../';
dir_software = [dir_root 'software/'];
dir_dataset  = '~/tmp2/git_repo/project/code/dataset/';
dir_result   = '~/tmp2/results_1';

%% Path
if ~exist(dir_result,'dir')
    mkdir(dir_result);
end

%% Add Path of 3rdparty
addpath(fullfile(dir_software,'3rdparty/yoon15iccv/'));

%% Add Path of Generate_Dataset
addpath('./generate_dataset/');
addpath('./save_parfor/');
addpath('./show_results/');

