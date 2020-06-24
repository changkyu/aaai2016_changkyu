function run_tv_one_method(method_name, MUs)

if ~exist('results/res_org.mat', 'file'), 
    error('Please run gen_org_images.m first!'); 
else
    load('results/res_org');
end

result = cell(length(IMGs), length(SIGMAs), length(MUs));

fprintf('### %s ###\n', method_name);

for ii = 1 : 1%length(IMGs)
    imgs_org{ii} = imread(IMGs{ii});

    for is = 1 : 1%length(SIGMAs)
        imgs_noise{ii,is} = double(imnoise(imgs_org{ii}, 'gaussian', 0, (SIGMAs(is)/256)^2));

        for im = 9 : 9%length(MUs)
            if strcmp(method_name, '1_admm')
                [res_im, res_his] = ...
                    tv_admm(imgs_noise{ii,is}, MUs(im) / 2, 1, double(imgs_org{ii}));
            elseif strcmp(method_name, '2_fadmmwr')
                [res_im, res_his] = ...
                    tv_fadmmrs(imgs_noise{ii,is}, MUs(im) / 2, 1, 0.999, double(imgs_org{ii}));
            elseif strcmp(method_name, '3_admmvp')
                [res_im, res_his] = ...
                    tv_admmvp(imgs_noise{ii,is}, MUs(im) / 2, 1, 10, 2, 2, 50, double(imgs_org{ii}));
            elseif strcmp(method_name, '4_admmap_2node')
                [res_im, res_his] = ...
                    tv_admmap_2node(imgs_noise{ii,is}, MUs(im) / 2, 1, 10, 300, double(imgs_org{ii}));
            elseif strcmp(method_name, '4_admmap_3node')
                [res_im, res_his] = ...
                    tv_admmap_3node(imgs_noise{ii,is}, MUs(im) / 2, 1, 50, double(imgs_org{ii}));
            end
            
            result{ii,is,im}.image = res_im;
            result{ii,is,im}.history = res_his;
            result{ii,is,im}.psnr = psnr(res_im, double(imgs_org{ii}));                

            fprintf('DONE: %s, sigma = %d, mu = %.2f (iter = %d, time = %.3f sec, method = %s)\n', ...
                IMGs{ii}, SIGMAs(is), MUs(im), res_his.eITER, res_his.eTIME, method_name);
        end
    end
end

% Save result
save(['results/res_' method_name '.mat'], 'result', 'MUs');

end