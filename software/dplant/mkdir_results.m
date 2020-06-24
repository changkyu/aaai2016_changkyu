function mkdir_results(models_desc, name_experiment)
% mkdir for results directories

env_setting;

for idx_model=1:numel(models_desc)
    model_desc = models_desc{idx_model};
    dir_model = fullfile(dir_result, model_desc.name);
    if ~exist(dir_model,'dir')
        mkdir(dir_model);
    end
    
    dir_experiment = fullfile( dir_model, name_experiment );
    if ~exist(dir_experiment,'dir')
        mkdir(dir_experiment);
    end    
end

end