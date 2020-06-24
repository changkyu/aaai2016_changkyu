%--------------------------------------------------------------------------
% Synthetic data result plot
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.02.05 (last modified on 2015/04/02)
%--------------------------------------------------------------------------
range = 10^3; % this value should be smaller or equal to max #iterations
fontsize = 18;
show_in_log_scale = false;

% load options if forgotten
load('data/synth/options.mat');

%% Figure 1. Different ETA
figure;
h = axes;
hold on;
hl = xlabel('iterations'); 
set(hl, 'FontSize',fontsize);
hl = ylabel('Objective function value'); 
set(hl, 'FontSize',fontsize);

min_val = inf;
max_val = -inf;
load('results/synth/N05_G002_E16.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,6); min_val]);
max_val = max([cm_to_show.objArray(1:range,6); max_val]);
load('results/synth/N05_G002_E12.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,6); min_val]);
max_val = max([cm_to_show.objArray(1:range,6); max_val]);
load('results/synth/N05_G002_E10.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,6); min_val]);
max_val = max([cm_to_show.objArray(1:range,6); max_val]);
load('results/synth/N05_G002_E08.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,6); min_val]);
max_val = max([cm_to_show.objArray(1:range,6); max_val]);
load('results/synth/N01_G001_E10.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,1); min_val]);
max_val = max([cm_to_show.objArray(1:range,1); max_val]);

load('results/synth/N05_G002_E16.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,6) = cm_to_show.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,6), ':r','LineWidth',2.5,'MarkerSize',7);

load('results/synth/N05_G002_E12.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,6) = cm_to_show.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,6), '-.g','LineWidth',2,'MarkerSize',7);

load('results/synth/N05_G002_E10.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,6) = cm_to_show.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,6), '--b','LineWidth',2,'MarkerSize',7);

load('results/synth/N05_G002_E08.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,6) = cm_to_show.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,6), '-.m','LineWidth',2,'MarkerSize',5);

load('results/synth/N01_G001_E10.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,1) = cm_to_show.objArray(1:range,1) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,1), '-k','LineWidth',2,'MarkerSize',7);

set(h,'XScale','log');
if show_in_log_scale
%     ylim([10^(3.52) 10^(3.525)]);
    ylim([min_val-100, max_val+100]);
    set(h,'YScale','log');
else
    ylim([min_val-10, max_val+10]);
end
set(h,'FontSize',fontsize);

legend('\eta=16','\eta=12','\eta=10','\eta=8','Centralized');

%% Figure 2. Different number of nodes
figure;
h = axes;
hold on;
hl = xlabel('iterations'); 
set(hl, 'FontSize',fontsize);
hl = ylabel('Objective function value'); 
set(hl, 'FontSize',fontsize);

min_val = inf;
max_val = -inf;
load('results/synth/N10_G002_E10.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,11); min_val]);
max_val = max([cm_to_show.objArray(1:range,11); max_val]);
load('results/synth/N08_G002_E10.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,9); min_val]);
max_val = max([cm_to_show.objArray(1:range,9); max_val]);
load('results/synth/N05_G002_E10.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,6); min_val]);
max_val = max([cm_to_show.objArray(1:range,6); max_val]);
load('results/synth/N02_G002_E10.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,3); min_val]);
max_val = max([cm_to_show.objArray(1:range,3); max_val]);
load('results/synth/N01_G001_E10.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,1); min_val]);
max_val = max([cm_to_show.objArray(1:range,1); max_val]);

load('results/synth/N10_G002_E10.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,11) = cm_to_show.objArray(1:range,11) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,11), ':r','LineWidth',2.5,'MarkerSize',5);

load('results/synth/N08_G002_E10.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,9) = cm_to_show.objArray(1:range,9) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,9), '-.g','LineWidth',2,'MarkerSize',7);

load('results/synth/N05_G002_E10.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,6) = cm_to_show.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,6), '--b','LineWidth',2,'MarkerSize',7);

load('results/synth/N02_G002_E10.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,3) = cm_to_show.objArray(1:range,3) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,3), '-.m','LineWidth',2,'MarkerSize',7);
               
load('results/synth/N01_G001_E10.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,1) = cm_to_show.objArray(1:range,1) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,1), '-k','LineWidth',2,'MarkerSize',7);

set(h,'XScale','log');
if show_in_log_scale
%     ylim([10^(3.52) 10^(3.525)]);
    ylim([min_val-100, max_val+100]);
    set(h,'YScale','log');
else
    ylim([min_val-10, max_val+10]);
end
set(h,'FontSize',fontsize);

legend('N=10','N=8','N=5','N=2','Centralized');

%% Figure 3. Different topology
figure;
h = axes;
hold on;
hl = xlabel('iterations'); 
set(hl, 'FontSize',fontsize);
hl = ylabel('Objective function value'); 
set(hl, 'FontSize',fontsize);

min_val = inf;
max_val = -inf;
load('results/synth/N05_G002_E10.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,6); min_val]);
max_val = max([cm_to_show.objArray(1:range,6); max_val]);
load('results/synth/N05_G003_E10.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,6); min_val]);
max_val = max([cm_to_show.objArray(1:range,6); max_val]);
load('results/synth/N05_G004_E10.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,6); min_val]);
max_val = max([cm_to_show.objArray(1:range,6); max_val]);
load('results/synth/N01_G001_E10.mat'); cm_to_show = cm4;
min_val = min([cm_to_show.objArray(1:range,1); min_val]);
max_val = max([cm_to_show.objArray(1:range,1); max_val]);

load('results/synth/N05_G002_E10.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,6) = cm_to_show.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,6), ':r','LineWidth',2.5,'MarkerSize',7);

load('results/synth/N05_G003_E10.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,6) = cm_to_show.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,6), '-.g','LineWidth',2,'MarkerSize',7);

load('results/synth/N05_G004_E10.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,6) = cm_to_show.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,6), '--b','LineWidth',2,'MarkerSize',7);

load('results/synth/N01_G001_E10.mat'); cm_to_show = cm4;
if show_in_log_scale
%    cm_to_show.objArray(1:range,1) = cm_to_show.objArray(1:range,1) - repmat(min_val, [range 1]);
end
plot(cm_to_show.objArray(1:range,1), '-k','LineWidth',2,'MarkerSize',7);

set(h,'XScale','log');
if show_in_log_scale
%     ylim([10^(3.52) 10^(3.525)]);
	ylim([min_val-100, max_val+100]);
    set(h,'YScale','log');
else
    ylim([min_val-10, max_val+10]);
end
set(h,'FontSize',fontsize);

legend('Ring','Star','Chain','Centralized');

%% Figure 4. Different rate of missing values (MAR)
figure;
h = axes;
hold on;
hl = xlabel('Missing Rates (%)'); 
set(hl, 'FontSize',fontsize);
hl = ylabel('Average Root Mean Squared error'); 
set(hl, 'FontSize',fontsize);

RMS_all = [];
for idk = 1 : length(RATEarr)
    RMS_cs1 = [];
    RMS_ds1 = [];
%     RMS_cs2 = [];
    RMS_ds2 = [];
%     RMS_cs3 = [];
%     RMS_cs4 = [];
    
    for idr = 1 : REPEATS
        load(sprintf('data/synth/data_R%d_i%02d.mat', idk, idr));
        load(sprintf('results/synth/R%d_i%02d_mar.mat', idk, idr));

        cm_c = cm1;
        RMS_c1 = calc_ppca_rms( Xmar, cm_c.W, cm_c.EZ, cm_c.MU );
        
        cm_d = cm3;
        RMS_d1 = calc_ppca_rms( Xmar, cm_d.W, cm_d.EZ, cm_d.MU );
        
%         cm_c = cm2;
%         RMS_c2 = calc_ppca_rms( Xmar, cm_c.mW, cm_c.mZ, cm_c.mMU );

        cm_d = cm4;
        RMS_d2 = calc_ppca_rms( Xmar, cm_d.W, cm_d.EZ, cm_d.MU );
%         
%         cm_c = cm5;
%         RMS_c3 = calc_ppca_rms( Xmar, cm_c.W, cm_c.EZ, cm_c.MU );
% 
%         cm_c = cm6;
%         RMS_c4 = calc_ppca_rms( Xmar, cm_c.mW, cm_c.mZ, cm_c.mMU );
        
        RMS_cs1 = [RMS_cs1 RMS_c1];
        RMS_ds1 = [RMS_ds1 RMS_d1];
%         RMS_cs2 = [RMS_cs2 RMS_c2];
        RMS_ds2 = [RMS_ds2 RMS_d2];
%         RMS_cs3 = [RMS_cs3 RMS_c3];
%         RMS_cs4 = [RMS_cs4 RMS_c4];
    end
    
    RMS_c1 = mean(RMS_cs1);
    RMS_d1 = mean(RMS_ds1);
%     RMS_c2 = mean(RMS_cs2);
    RMS_d2 = mean(RMS_ds2);
%     RMS_c3 = mean(RMS_cs3);
%     RMS_c4 = mean(RMS_cs4);

    RMS_all = [RMS_all; RMS_c1, RMS_d1, RMS_d2];
%     RMS_all = [RMS_all; RMS_c1, RMS_d1, RMS_c2, RMS_d2, RMS_c3, RMS_c4];
end
bar(RMS_all);

%ylim([0.4, 0.48]);
set(gca, 'XLim', [.5 5.5]);
set(gca, 'XTick', [1:5]);
set(gca, 'XTickLabel', ...
    {sprintf('%.2f', RATEarr{1}*100); ...
    sprintf('%.2f', RATEarr{2}*100); ...
    sprintf('%.2f', RATEarr{3}*100); ...
    sprintf('%.2f', RATEarr{4}*100); ...
    sprintf('%.2f', RATEarr{5}*100)});

set(h,'FontSize',fontsize);
legend('PPCA', 'D-PPCA', 'D-PPCA-ANT');
% legend('PPCA', 'D-PPCA', 'BPCA', 'D-BPCA', 'PPCA(IR)', 'VBPCA(IR)');
set(gca,'YGrid','on');

%% Figure 5. Different rate of missing values (MNAR)
figure;
h = axes;
hold on;
hl = xlabel('Missing Rates (%)'); 
set(hl, 'FontSize',fontsize);
hl = ylabel('Average Root Mean Squared error'); 
set(hl, 'FontSize',fontsize);

RMS_all = [];
for idk = 1 : length(RATEarr)
    RMS_cs1 = [];
    RMS_ds1 = [];
%     RMS_cs2 = [];
    RMS_ds2 = [];
%     RMS_cs3 = [];
%     RMS_cs4 = [];    
    
    for idr = 1 : REPEATS
        load(sprintf('data/synth/data_R%d_i%02d.mat', idk, idr));
        load(sprintf('results/synth/R%d_i%02d_mnar.mat', idk, idr));

        cm_c = cm1;
        RMS_c1 = calc_ppca_rms( Xmnar, cm_c.W, cm_c.EZ, cm_c.MU );
        
        cm_d = cm3;
        RMS_d1 = calc_ppca_rms( Xmnar, cm_d.W, cm_d.EZ, cm_d.MU );
        
%         cm_c = cm2;
%         RMS_c2 = calc_ppca_rms( Xmnar, cm_c.mW, cm_c.mZ, cm_c.mMU );

        cm_d = cm4;
        RMS_d2 = calc_ppca_rms( Xmnar, cm_d.W, cm_d.EZ, cm_d.MU );
% 
%         cm_c = cm5;
%         RMS_c3 = calc_ppca_rms( Xmar, cm_c.W, cm_c.EZ, cm_c.MU );
% 
%         cm_c = cm6;
%         RMS_c4 = calc_ppca_rms( Xmar, cm_c.mW, cm_c.mZ, cm_c.mMU );        

        RMS_cs1 = [RMS_cs1 RMS_c1];
        RMS_ds1 = [RMS_ds1 RMS_d1];
%         RMS_cs2 = [RMS_cs2 RMS_c2];
        RMS_ds2 = [RMS_ds2 RMS_d2];
%         RMS_cs3 = [RMS_cs3 RMS_c3];
%         RMS_cs4 = [RMS_cs4 RMS_c4];        
    end
    
    RMS_c1 = mean(RMS_cs1);
    RMS_d1 = mean(RMS_ds1);
%     RMS_c2 = mean(RMS_cs2);
    RMS_d2 = mean(RMS_ds2);
%     RMS_c3 = mean(RMS_cs3);
%     RMS_c4 = mean(RMS_cs4);    
    
    RMS_all = [RMS_all; RMS_c1, RMS_d1, RMS_d2];
%     RMS_all = [RMS_all; RMS_c1, RMS_d1, RMS_c2, RMS_d2, RMS_c3, RMS_c4];
end
bar(RMS_all);

%ylim([0.4, 0.48]);
set(gca, 'XLim', [.5 5.5]);
set(gca, 'XTick', [1:5]);
set(gca, 'XTickLabel', ...
    {sprintf('%.2f', RATEarr{1}*100); ...
    sprintf('%.2f', RATEarr{2}*100); ...
    sprintf('%.2f', RATEarr{3}*100); ...
    sprintf('%.2f', RATEarr{4}*100); ...
    sprintf('%.2f', RATEarr{5}*100)});

set(h,'FontSize',fontsize);
legend('PPCA', 'D-PPCA', 'D-PPCA-ANT');
% legend('PPCA', 'D-PPCA', 'BPCA', 'D-BPCA', 'PPCA(IR)', 'VBPCA(IR)');
set(gca,'YGrid','on');
