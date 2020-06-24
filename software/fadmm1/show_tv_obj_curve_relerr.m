function [best_relerr, worst_relerr, total_stats] = show_tv_obj_curve_relerr(method_names, MUs)

method_names_for_legend = {'ADMM', 'F-ADMM-RS', 'ADMM-VP', 'Ours (ADMM-AP)'};
markers = {'s', '^', '+', 'o'}.'; % for some reason, we MUST transpose...

load('results/res_org.mat');
best_relerr = cell(length(IMGs), length(SIGMAs), length(method_names));
worst_relerr = cell(length(IMGs), length(SIGMAs), length(method_names));

total_stats.time = zeros(length(IMGs), length(SIGMAs), length(method_names));
total_stats.iter = zeros(length(IMGs), length(SIGMAs), length(method_names));
total_stats.relerr = zeros(length(IMGs), length(SIGMAs), length(method_names));
total_stats.relerrF = zeros(length(IMGs), length(SIGMAs), length(method_names), length(MUs));
total_stats.relerrO = zeros(length(IMGs), length(SIGMAs), length(method_names), length(MUs));
total_stats.relerrI = zeros(length(IMGs), length(SIGMAs), length(method_names), length(MUs));

colorOrder = get(gca, 'ColorOrder');
close gcf;
colorOrder = colorOrder([3 4 5 1],:);  

for ii = 1 : length(IMGs)
    for is = 1 : length(SIGMAs)
        for im = 1 : length(MUs) %[2 4 6 8 10]
            h = figure;
            
            max_iter = 0;
            max_time = 0;
            toPlotY = [];
            toPlotX = [];
            toPlotCx = [];
            toPlotCy = [];
            toPlotCobjx = cell(length(method_names),1);
            toPlotCobjzx = cell(length(method_names),1);
            toPlotCobjzy = cell(length(method_names),1);
            
            for imt = 1:length(method_names)    
                if exist(sprintf('results/res_%s.mat', method_names{imt}), 'file')
                    load(sprintf('results/res_%s.mat', method_names{imt}));
                else
                    error([sprintf('results/res_%s.mat', method_names{imt}) ' does not exist!']);
                end
                
                %----------------------------------------------------------
                % accumulate total stats
                total_stats.time(ii,is,imt) = total_stats.time(ii,is,imt) + result{ii,is,im}.history.eTIME;
                total_stats.iter(ii,is,imt) = total_stats.iter(ii,is,imt) + result{ii,is,im}.history.eITER;
                total_stats.relerr(ii,is,imt) = total_stats.relerr(ii,is,imt) + ...
                    result{ii,is,im}.history.relerr(result{ii,is,im}.history.eITER);
                total_stats.relerrF(ii,is,imt,im) = ...
                    norm(result{ii,is,im}.image - double(imgs_org{ii}), 'fro') / norm(double(imgs_org{ii}), 'fro');
                total_stats.relerrO(ii,is,imt,im) = ...
                    norm(imgs_noise{ii,is} - double(imgs_org{ii}), 'fro') / norm(double(imgs_org{ii}), 'fro');
                total_stats.relerrI(ii,is,imt,im) = ...
                    (total_stats.relerrO(ii,is,imt,im) - total_stats.relerrF(ii,is,imt,im)) / ...
                    result{ii,is,im}.history.eITER;
                %----------------------------------------------------------

                % calculate time-normalized iterations (time / iteration)
                % psnr_gain = (best_psnr - worst_psnr)
                % psnr_at_iter = (psnr_earned_or_lost_in_that_patch) / iterations
  
                max_iter = max(max_iter, result{ii,is,im}.history.eITER);
                max_time = max(max_time, result{ii,is,im}.history.eTIME);
                
                % pad nan to non-existing values
                toPlotYtmp = nan(imt,max_iter);
                for timt = 1 : imt-1
                    toPlotYtmp(timt,1:length(toPlotY(timt,:))) = toPlotY(timt,:);
                end
                toPlotYtmp(imt,1:result{ii,is,im}.history.eITER) = ...
                    result{ii,is,im}.history.relerr;
                toPlotY = toPlotYtmp;
                
                % pad nan to non-existing values
                toPlotXtmp = nan(imt,max_iter);
                for timt = 1 : imt-1
                    toPlotXtmp(timt,1:length(toPlotX(timt,:))) = toPlotX(timt,:);
                end
                for idi = 1 : result{ii,is,im}.history.eITER
                    toPlotXtmp(imt,idi) = idi * result{ii,is,im}.history.eTIME / result{ii,is,im}.history.eITER;
                end
                toPlotX = toPlotXtmp;
                
                % get consensus
                toPlotCx = [toPlotCx; result{ii,is,im}.history.consensus(result{ii,is,im}.history.eITER).xzx];
                toPlotCy = [toPlotCy; result{ii,is,im}.history.consensus(result{ii,is,im}.history.eITER).xzy];
                toPlotCobjx{imt} = [result{ii,is,im}.history.consensus(:).objx];
                toPlotCobjzx{imt} = [result{ii,is,im}.history.consensus(:).objzx];
                toPlotCobjzy{imt} = [result{ii,is,im}.history.consensus(:).objzy];
                
                % find which mu's the best for each method
                if isfield(best_relerr{ii, is, imt}, 'relerr')
                    if best_relerr{ii, is, imt}.relerr > toPlotY(imt, result{ii,is,im}.history.eITER)
                        best_relerr{ii, is, imt}.relerr = toPlotY(imt, result{ii,is,im}.history.eITER);
                        best_relerr{ii, is, imt}.idx = im;
                    end
                else
                    best_relerr{ii, is, imt}.relerr = toPlotY(imt, result{ii,is,im}.history.eITER);
                    best_relerr{ii, is, imt}.idx = im;
                end
                
                % also find the worst...
                if isfield(worst_relerr{ii, is, imt}, 'relerr')
                    if worst_relerr{ii, is, imt}.relerr < toPlotY(imt, result{ii,is,im}.history.eITER)
                        worst_relerr{ii, is, imt}.relerr = toPlotY(imt, result{ii,is,im}.history.eITER);
                        worst_relerr{ii, is, imt}.idx = im;
                    end
                else
                    worst_relerr{ii, is, imt}.relerr = toPlotY(imt, result{ii,is,im}.history.eITER);
                    worst_relerr{ii, is, imt}.idx = im;
                end
                
                % add up total time...
                if isfield(best_relerr{ii, is, imt}, 'totalTime')
                    best_relerr{ii, is, imt}.totalTime = ...
                        best_relerr{ii, is, imt}.totalTime + result{ii,is,im}.history.eTIME;
                else
                    best_relerr{ii, is, imt}.totalTime = ...
                        result{ii,is,im}.history.eTIME;
                end
            end
            
            %% Plot relative error against iterations and time
            subplot(1,2,1);
            p = semilogy(toPlotY(:,1:max_iter)','LineWidth',2);
            for imt = 1 : length(method_names)
                p(imt).Color = colorOrder(imt,:);
            end
            set(p, {'Marker'}, markers);

            %legend(method_names_for_legend, 'Location', 'northeast');
            xlabel('Iterations');
            ylabel('Relative error: || x - x* || / || x* ||');
            
            subplot(1,2,2);
            p = semilogy(toPlotX(:,1:max_iter)', toPlotY(:,1:max_iter)','LineWidth',2);
            for imt = 1 : length(method_names)
                p(imt).Color = colorOrder(imt,:);
            end
            set(p, {'Marker'}, markers);
            
            legend(method_names_for_legend, 'Location', 'northeast');
            xlabel('Seconds');
            ylabel('Relative error: || x - x* || / || x* ||');
            
            %% Plot consensus error as bar
            
            consensusV = axes('Position', [0.3 0.75 0.15 0.15]);
            barGraph = bar(diag(toPlotCx), 'stacked');
            for imt = 1 : length(method_names)
                % change color so that it corresponds to semilogy plots
                set(barGraph(imt), 'FaceColor', colorOrder(imt,:));
            end
            set(consensusV, 'Box', 'off');
            ylabel('|| Dx * x - zx ||');
            set(gca,'xtick',[]);
            
            consensusV = axes('Position', [0.3 0.5 0.15 0.15]);
            barGraph = bar(diag(toPlotCy), 'stacked');
            for imt = 1 : length(method_names)
                % change color so that it corresponds to semilogy plots
                set(barGraph(imt), 'FaceColor', colorOrder(imt,:));
            end            
            set(consensusV, 'Box', 'off');
            ylabel('|| Dy * x - zy ||');
            set(gca,'xtick',[]);            
            
            %% Final save
            set(h, 'Position', [0 0 1280 540])
            set(gcf, 'Color', 'white'); % white bckgr
            saveas( gcf, ...   % figure handle
                sprintf('results/relerr/plot_i%d_s%d_m%d.eps', ii,is,im),... % name of output file without extension
                'epsc');           % file format
            saveas( gcf, ...   % figure handle
                sprintf('results/relerr/plot_i%d_s%d_m%d.png', ii,is,im),... % name of output file without extension
                'png');           % file format
            close gcf;         
            
            fprintf('DONE ii = %d, is = %d, im = %d\n', ii, is, im);
        end
    end
end

subplot2 = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);

for ii = 1 : length(IMGs)
    for is = 1 : length(SIGMAs)
        h = figure;            

        for imt = 1 : length(method_names)    
            if exist(sprintf('results/res_%s.mat', method_names{imt}), 'file')
                load(sprintf('results/res_%s.mat', method_names{imt}));
            else
                error([sprintf('results/res_%s.mat', method_names_for_legend{imt}) ' does not exist!']);
            end
            
            i_best = best_relerr{ii, is, imt}.idx;
            i_worst = worst_relerr{ii, is, imt}.idx;

            subplot2(2, length(method_names), imt);
            imagesc(result{ii,is,i_best}.image);
            colormap gray;
            set(gca,'ytick',[]);            
            set(gca,'yticklabel',[]);
            set(gca,'xtick',[]);
            set(gca,'xticklabel',[]);
            
%             title(sprintf('+%.3f / %.3f s', ...
%                 best_relerr{ii, is, imt}.relerr - worst_relerr{ii, is, imt}.relerr, ...
%                 best_relerr{ii, is, imt}.totalTime));

            subplot2(2, length(method_names), imt + length(method_names));
            imagesc(result{ii,is,i_worst}.image);
            colormap gray;
            set(gca,'ytick',[]);            
            set(gca,'yticklabel',[]);
            set(gca,'xtick',[]);
            xlabel(method_names_for_legend{imt});
            
            fprintf('DONE ii = %d, is = %d, im = %d\n', ii, is, im);
            
            fprintf('%s / Best: %f / Worst: %f / Diff: %f / Time : %f / Diff (to proposed best): %f / Diff (to proposed worst): %f\n', ...
                method_names{imt}, ...
                best_relerr{ii, is, imt}.relerr, ...
                worst_relerr{ii, is, imt}.relerr, ...
                best_relerr{ii, is, imt}.relerr - worst_relerr{ii, is, imt}.relerr, ...
                best_relerr{ii, is, imt}.totalTime, ...
                best_relerr{ii, is, imt}.relerr - best_relerr{ii, is, length(method_names)}.relerr, ...
                worst_relerr{ii, is, imt}.relerr - worst_relerr{ii, is, length(method_names)}.relerr);
        end
        
        set(h, 'Position', [0 0 1280 720]);
        set(gcf, 'Color', 'white'); % white bckgr
        saveas( gcf, ...   % figure handle
            sprintf('results/gain_i%d_s%d.eps', ii,is),...
            'epsc');           % file format
        saveas( gcf, ...   % figure handle
            sprintf('results/gain_i%d_s%d.png', ii,is),...
            'png');           % file format
        close gcf;     
    end
end

end
