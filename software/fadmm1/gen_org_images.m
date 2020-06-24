function gen_org_images()

IMGs = {'data/cameraman.tif', 'data/barbara.png', 'data/boat.png'};
SIGMAs = [20, 50];

imgs_org = cell(length(IMGs), 1);
imgs_noise = cell(length(IMGs), length(SIGMAs));

for ii = 1 : length(IMGs)
    imgs_org{ii} = imread(IMGs{ii});

    for is = 1 : length(SIGMAs)
        imgs_noise{ii,is} = double(imnoise(imgs_org{ii}, 'gaussian', 0, (SIGMAs(is)/256)^2));
    end
end

save('results/res_org.mat', 'imgs_org', 'imgs_noise', 'IMGs', 'SIGMAs');

end
