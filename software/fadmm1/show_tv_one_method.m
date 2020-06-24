function show_tv_one_method(method_name, MUs, best_relerr, worst_relerr, idm)

load('results/res_org.mat');
if exist(sprintf('results/res_%s.mat', method_name), 'file')
    load(sprintf('results/res_%s.mat', method_name));
else
    error([sprintf('results/res_%s.mat', method_name) ' does not exist!']);
end

%% Show result of this method
h = figure;
rows = length(IMGs) * length(SIGMAs);
cols = length([2 4 6 8 10]);%length(MUs);
for ii = 1 : length(IMGs)
    for is = 1 : length(SIGMAs)
        for im = [2 4 6 8 10]%1 : length(MUs)
            subplot(rows, cols, (ii - 1) * (cols*2) + (is - 1) * cols + im/2);
            
            i_best = best_relerr{ii, is, idm}.idx;
            i_worst = worst_relerr{ii, is, idm}.idx;            
            
            imagesc(result{ii,is,im}.image);
            colormap(gray); axis equal; 
            xlim([0 size(result{ii,is,im}.image, 2)]);
            ylim([0 size(result{ii,is,im}.image, 1)]);
            
            if im == i_best
                title({['\mu = ' sprintf('%.2f', MUs(im)) ...
                    sprintf(', err = %.4f', result{ii,is,im}.history.relerr(result{ii,is,im}.history.eITER))], ...
                    ['t = ' sprintf('%d', result{ii,is,im}.history.eITER) ...
                    sprintf(' (%.3f s)', result{ii,is,im}.history.eTIME)]}, ...
                    'Color', 'b', 'FontWeight','bold');
            elseif im == i_worst
                title({['\mu = ' sprintf('%.2f', MUs(im)) ...
                    sprintf(', err = %.4f', result{ii,is,im}.history.relerr(result{ii,is,im}.history.eITER))], ...
                    ['t = ' sprintf('%d', result{ii,is,im}.history.eITER) ...
                    sprintf(' (%.3f s)', result{ii,is,im}.history.eTIME)]}, ...
                    'Color', 'r');
            else
                title({['\mu = ' sprintf('%.2f', MUs(im)) ...
                    sprintf(', err = %.4f', result{ii,is,im}.history.relerr(result{ii,is,im}.history.eITER))], ...
                    ['t = ' sprintf('%d', result{ii,is,im}.history.eITER) ...
                    sprintf(' (%.3f s)', result{ii,is,im}.history.eTIME)]});
            end
            % hide ticks
            set(gca,'xtick',[]);
            set(gca,'xticklabel',[]);
            set(gca,'ytick',[]);
            if im == 1
                ylabel(['\sigma = ' sprintf('%d', SIGMAs(is))]);
            end
        end
    end
end
tightfig;
set(h, 'Position', [0 0 980 970]);
set(gcf, 'Color', 'white'); % white bckgr
export_fig( gcf, ...   % figure handle
    sprintf('results/res_%s.eps', method_name),... % name of output file without extension
    '-eps');           % file format
close gcf;

end
