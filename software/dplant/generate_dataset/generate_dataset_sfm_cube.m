function [X, n_X, dim_X, dim_Z, GT] = generate_dataset_sfm_cube( noise, idr )

env_setting;
dir_cube = fullfile(dir_dataset,'cube');
if ~exist(dir_cube, 'dir')
    synthetic_cube_data_gen_ex(dir_dataset);
end

load(fullfile(dir_cube, sprintf('%.5f_%03d.mat', noise, idr)));

X = measurement_matrix';
[n_X dim_X] = size(X);
dim_Z = size(GT,2);

end
