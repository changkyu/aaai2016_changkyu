function [X, n_X, dim_X, dim_Z, GT] = generate_dataset_sfm_hopkins( idx )

env_setting;
dir_hopkins = fullfile(dir_dataset,'hopkins');
if ~exist(dir_hopkins, 'dir')
    error(['There is no such a directory: ' dir_hopkins]);
end

load( fullfile( dir_hopkins,sprintf('%03d.mat',idx)) );

X = measurement_matrix';
[n_X dim_X] = size(X);
dim_Z = 3;
[~, GT, ~, ~, ~, ~] = sfm_affine(X');
GT = GT';

end


