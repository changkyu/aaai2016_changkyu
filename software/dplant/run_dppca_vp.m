rng(0,'v5uniform');

nmatlabs = 4;
if verLessThan('matlab', '8.2')
    if matlabpool('size') == 0 
       matlabpool(nmatlabs) ;
    end
else
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        parpool(nmatlabs);
    end
end

clearvars;
env_setting;
range_models = [8];
range_experiments = [4];
experiments;
