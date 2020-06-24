nmatlabs = 5;
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
range_models = 2;
experiments;

clearvars;
env_setting;
range_models = 3;
experiments;
