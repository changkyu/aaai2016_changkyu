function show_tv_obj_curve_subopt(method_names, MUs)

method_names_for_legend = {'ADMM', 'F-ADMM-RS', 'ADMM-VP', 'Proposed'};
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
  
                toPlotY = abs([toPlotY; result{ii,is,im}.history.subopt]);
                max_iter = max(max_iter, result{ii,is,im}.history.eITER);
                max_time = max(max_time, result{ii,is,im}.history.eTIME);
                
                toPlotX = abs([toPlotX; zeros(size(result{ii,is,im}.history.subopt))]);
                for idi = 1 : result{ii,is,im}.history.eITER
                    toPlotX(imt,idi) = idi * result{ii,is,im}.history.eTIME / result{ii,is,im}.history.eITER;
                end
                toPlotX(imt,idi+1:end) = Inf;
                
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
            ylabel('Log absolute suboptimality: log(abs(F(x) - F(x*)))');
            
            subplot(1,2,2);
            p = semilogy(toPlotX(:,1:max_iter)', toPlotY(:,1:max_iter)');
            set(p, {'Marker'}, markers);
            
            %legend(method_names_for_legend, 'Location', 'northwest');
            xlabel('Seconds');
            ylabel('Log absolute suboptimality: log(abs(F(x) - F(x*)))');
            
            %% Final save
            set(h, 'Position', [0 0 1280 720])
            set(gcf, 'Color', 'white'); % white bckgr
            saveas( gcf, ...   % figure handle
                sprintf('results/subopt/plot_i%d_s%d_m%d.eps', ii,is,im),... % name of output file without extension
                'epsc');           % file format
            close gcf;         
            
            fprintf('DONE ii = %d, is = %d, im = %d\n', ii, is, im);
        end
    end
end

end
