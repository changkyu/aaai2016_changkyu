function show_tv_obj_curve(method_names, MUs)

method_names_for_legend = {'ADMM', 'F-ADMM-RS', 'ADMM-VP', 'Ours (ADMM-AP)'};
markers = {'s', 'x', '+', 'o'}.'; % for some reason, we MUST transpose...

load('results/res_org.mat');
best_psnr = cell(length(IMGs), length(SIGMAs), length(method_names));
worst_psnr = cell(length(IMGs), length(SIGMAs), length(method_names));

for ii = 1 : length(IMGs)
    for is = 1 : length(SIGMAs)
        for im = 1 : length(MUs)
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
                    result{ii,is,im}.history.objval;
                toPlotY = toPlotYtmp;                
                
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
                if isfield(best_psnr{ii, is, imt}, 'psnr')
                    if best_psnr{ii, is, imt}.psnr < result{ii,is,im}.psnr
                        best_psnr{ii, is, imt}.psnr = result{ii,is,im}.psnr;
                        best_psnr{ii, is, imt}.idx = im;
                    end
                else
                    best_psnr{ii, is, imt}.psnr = result{ii,is,im}.psnr;
                    best_psnr{ii, is, imt}.idx = im;
                end
                
                % also find the worst...
                if isfield(worst_psnr{ii, is, imt}, 'psnr')
                    if worst_psnr{ii, is, imt}.psnr > result{ii,is,im}.psnr
                        worst_psnr{ii, is, imt}.psnr = result{ii,is,im}.psnr;
                        worst_psnr{ii, is, imt}.idx = im;
                    end
                else
                    worst_psnr{ii, is, imt}.psnr = result{ii,is,im}.psnr;
                    worst_psnr{ii, is, imt}.idx = im;
                end
                
                % add up total time...
                if isfield(best_psnr{ii, is, imt}, 'totalTime')
                    best_psnr{ii, is, imt}.totalTime = ...
                        best_psnr{ii, is, imt}.totalTime + result{ii,is,im}.history.eTIME;
                else
                    best_psnr{ii, is, imt}.totalTime = ...
                        result{ii,is,im}.history.eTIME;
                end
            end
            
            %% Plot total objective value against iterations and time
            subplot(1,2,1);
            p = semilogy(toPlotY(:,1:max_iter)');
            set(p, {'Marker'}, markers);

            legend(method_names_for_legend, 'Location', 'northwest');
            xlabel('Iterations');
            ylabel('Objective value');
            
            subplot(1,2,2);
            p = semilogy(toPlotX(:,1:max_iter)', toPlotY(:,1:max_iter)');
            set(p, {'Marker'}, markers);
            
            %legend(method_names_for_legend, 'Location', 'northwest');
            xlabel('Seconds');
            ylabel('Objective value');
            
            %% Plot objective partitions
            toPlotC = zeros(length(method_names), max_iter);
            for imt = 1 : length(method_names)
                toPlotC(imt,:) = [toPlotCobjx{imt}, nan(1,max_iter - length(toPlotCobjx{imt}))];
            end

            consensusV = axes('Position', [0.34 0.75 0.12 0.12]);
            semilogy(toPlotC');
            set(consensusV, 'Box', 'off');
            ylabel('\mu * 0.5 * || x - b ||');
            set(gca,'xtick',[]);
            %xlabel('iterations');
            %calc_diffs(toPlotC)

            toPlotC = zeros(length(method_names), max_iter);
            for imt = 1 : length(method_names)
                toPlotC(imt,:) = [toPlotCobjzx{imt}, nan(1,max_iter - length(toPlotCobjzx{imt}))];
            end
            
            consensusV = axes('Position', [0.34 0.58 0.12 0.12]);
            semilogy(toPlotC');
            set(consensusV, 'Box', 'off');
            ylabel('|| zx ||');
            set(gca,'xtick',[]);
            %xlabel('iterations');
            %calc_diffs(toPlotC)
            
            toPlotC = zeros(length(method_names), max_iter);
            for imt = 1 : length(method_names)
                toPlotC(imt,:) = [toPlotCobjzy{imt}, nan(1,max_iter - length(toPlotCobjzy{imt}))];
            end
            
            consensusV = axes('Position', [0.34 0.38 0.12 0.12]);
            semilogy(toPlotC');
            set(consensusV, 'Box', 'off');
            ylabel('|| zy ||');
            set(gca,'xtick',[]);
            xlabel('iterations');        
            %calc_diffs(toPlotC)
            
            %% Plot consensus error as bar
            colorOrder = get(gca, 'ColorOrder');
            
            consensusV = axes('Position', [0.75 0.75 0.15 0.15]);
            barGraph = bar(diag(toPlotCx), 'stacked');
            for imt = 1 : length(method_names)
                % change color so that it corresponds to semilogy plots
                set(barGraph(imt), 'FaceColor', colorOrder(imt,:));
            end
            set(consensusV, 'Box', 'off');
            ylabel('|| Dx * x - zx ||');
            set(gca,'xtick',[]);
            
            consensusV = axes('Position', [0.75 0.58 0.15 0.15]);
            barGraph = bar(diag(toPlotCy), 'stacked');
            for imt = 1 : length(method_names)
                % change color so that it corresponds to semilogy plots
                set(barGraph(imt), 'FaceColor', colorOrder(imt,:));
            end            
            set(consensusV, 'Box', 'off');
            ylabel('|| Dy * x - zy ||');
            set(gca,'xtick',[]);
            
            %% Final save
            set(h, 'Position', [0 0 1280 720])
            set(gcf, 'Color', 'white'); % white bckgr
            tightfig;
            saveas( gcf, ...   % figure handle
                sprintf('results/plot_i%d_s%d_m%d.eps', ii,is,im),... % name of output file without extension
                'epsc');           % file format
            close gcf;         
            
            fprintf('DONE ii = %d, is = %d, im = %d\n', ii, is, im);
        end
    end
end

for ii = 1 : length(IMGs)
    for is = 1 : length(SIGMAs)
        h = figure;            

        for imt = 1:length(method_names)    
            if exist(sprintf('results/res_%s.mat', method_names{imt}), 'file')
                load(sprintf('results/res_%s.mat', method_names{imt}));
            else
                error([sprintf('results/res_%s.mat', method_names_for_legend{imt}) ' does not exist!']);
            end
            
            i_best = best_psnr{ii, is, imt}.idx;
            i_worst = worst_psnr{ii, is, imt}.idx;

            subplot(2, length(method_names), imt);
            imagesc(result{ii,is,i_best}.image);
            colormap gray;
            set(gca,'ytick',[]);            
            set(gca,'yticklabel',[]);
            set(gca,'xtick',[]);
            set(gca,'xticklabel',[]);
            
            title(sprintf('+%.3f / %.3f s', ...
                best_psnr{ii, is, imt}.psnr - worst_psnr{ii, is, imt}.psnr, ...
                best_psnr{ii, is, imt}.totalTime));

            subplot(2, length(method_names), imt + length(method_names));
            imagesc(result{ii,is,i_worst}.image);
            colormap gray;
            set(gca,'ytick',[]);            
            set(gca,'yticklabel',[]);
            set(gca,'xtick',[]);
            xlabel(method_names_for_legend{imt});
            
            fprintf('DONE ii = %d, is = %d, im = %d\n', ii, is, im);
            
            fprintf('%s / Best: %f / Worst: %f / Diff: %f / Time : %f / Diff (to proposed best): %f / Diff (to proposed worst): %f\n', ...
                method_names{imt}, ...
                best_psnr{ii, is, imt}.psnr, ...
                worst_psnr{ii, is, imt}.psnr, ...
                best_psnr{ii, is, imt}.psnr - worst_psnr{ii, is, imt}.psnr, ...
                best_psnr{ii, is, imt}.totalTime, ...
                best_psnr{ii, is, imt}.psnr - best_psnr{ii, is, length(method_names)}.psnr, ...
                worst_psnr{ii, is, imt}.psnr - worst_psnr{ii, is, length(method_names)}.psnr);
        end
        
        set(h, 'Position', [0 0 1280 720]);
        set(gcf, 'Color', 'white'); % white bckgr
        tightfig;
        saveas( gcf, ...   % figure handle
            sprintf('results/gain_i%d_s%d.eps', ii,is),...
            'epsc');           % file format
        close gcf;     
    end
end

end

%% DEBUG function
function diffs = calc_diffs(toPlotC)

[nmethods, ~] = size(toPlotC);

diffs = zeros(nmethods, 1);

for imt = 1 : nmethods
    indices = find(~isnan(toPlotC(imt,:)));
    diffs(imt) = toPlotC(imt, indices(end));
end

diffs = diffs - diffs(nmethods);

end
