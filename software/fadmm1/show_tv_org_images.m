function show_tv_org_images()

load('results/res_org');
IMG_NAMEs = {'cameraman', 'barbara', 'boat'};

%% Show corrupted image
h = figure;
rows = length(IMGs);
cols = length(SIGMAs) + 1;
for ii = 1 : rows
    for is = 1 : cols
        if is == 1
            subplot(rows, cols, (ii - 1) * cols + is);
            imagesc(imgs_org{ii});
            title(IMG_NAMEs{ii});
        else
            subplot(rows, cols, (ii - 1) * cols + is);
            imagesc(imgs_noise{ii,is-1});
            title(['\sigma = ' sprintf('%d', SIGMAs(is-1))]);
        end
        colormap(gray); axis equal; 
        xlim([0 size(imgs_org{ii}, 2)]);
        ylim([0 size(imgs_org{ii}, 1)]);
        % hide ticks
        set(gca,'xtick',[]);
        set(gca,'xticklabel',[]);
        set(gca,'ytick',[]);
        set(gca,'yticklabel',[]);
    end
end
tightfig;
set(h, 'Position', [0 0 512 512]);
set(h, 'Color', 'white'); % white bckgr
export_fig(h, 'results/org_all.eps', '-eps');
close gcf;

end