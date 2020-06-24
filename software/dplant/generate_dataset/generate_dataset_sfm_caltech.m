function [X, n_X, dim_X, dim_Z, GT] = generate_dataset_sfm_caltech( idx )

env_setting;
dir_caltech = fullfile(dir_dataset,'caltech');
if ~exist(dir_caltech, 'dir')
    error(['There is no such a directory: ' dir_caltech]);
end

load( fullfile( dir_caltech,sprintf('%03d.mat',idx)) );

X = measurement_matrix';
[n_X dim_X] = size(X);
dim_Z = 3;
[~, GT, ~, ~, ~, ~] = sfm_affine(X');
GT = GT';
end

