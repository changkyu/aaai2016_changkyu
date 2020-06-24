%--------------------------------------------------------------------------
% Synthetic data result plot
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.02.05 (last modified on 2012/03/14)
%--------------------------------------------------------------------------
fontsize = 18;
range = 10^4;

show_in_log_scale = true;

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
load('result/synth_dppca/dppca_N05_G002_E16.mat');
min_val = min([cm.objArray(1:range,6); min_val]);
max_val = max([cm.objArray(1:range,6); max_val]);
load('result/synth_dppca/dppca_N05_G002_E12.mat');
min_val = min([cm.objArray(1:range,6); min_val]);
max_val = max([cm.objArray(1:range,6); max_val]);
load('result/synth_dppca/dppca_N05_G002_E10.mat');
min_val = min([cm.objArray(1:range,6); min_val]);
max_val = max([cm.objArray(1:range,6); max_val]);
load('result/synth_dppca/dppca_N05_G002_E08.mat');
min_val = min([cm.objArray(1:range,6); min_val]);
max_val = max([cm.objArray(1:range,6); max_val]);
load('result/synth_dppca/dppca_N01_G001_E10.mat');
min_val = min([cm.objArray(1:range,1); min_val]);
max_val = max([cm.objArray(1:range,1); max_val]);

load('result/synth_dppca/dppca_N05_G002_E16.mat');
if show_in_log_scale
%    cm.objArray(1:range,6) = cm.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,6), ':r','LineWidth',2.5,'MarkerSize',7);

load('result/synth_dppca/dppca_N05_G002_E12.mat');
if show_in_log_scale
%    cm.objArray(1:range,6) = cm.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,6), '-.g','LineWidth',2,'MarkerSize',7);

load('result/synth_dppca/dppca_N05_G002_E10.mat');
if show_in_log_scale
%    cm.objArray(1:range,6) = cm.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,6), '--b','LineWidth',2,'MarkerSize',7);

load('result/synth_dppca/dppca_N05_G002_E08.mat');
if show_in_log_scale
%    cm.objArray(1:range,6) = cm.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,6), '-.m','LineWidth',2,'MarkerSize',5);

load('result/synth_dppca/dppca_N01_G001_E10.mat');
if show_in_log_scale
%    cm.objArray(1:range,1) = cm.objArray(1:range,1) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,1), '-k','LineWidth',2,'MarkerSize',7);

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
load('result/synth_dppca/dppca_N10_G002_E10.mat');
min_val = min([cm.objArray(1:range,11); min_val]);
max_val = max([cm.objArray(1:range,11); max_val]);
load('result/synth_dppca/dppca_N08_G002_E10.mat');
min_val = min([cm.objArray(1:range,9); min_val]);
max_val = max([cm.objArray(1:range,9); max_val]);
load('result/synth_dppca/dppca_N05_G002_E10.mat');
min_val = min([cm.objArray(1:range,6); min_val]);
max_val = max([cm.objArray(1:range,6); max_val]);
load('result/synth_dppca/dppca_N02_G002_E10.mat');
min_val = min([cm.objArray(1:range,3); min_val]);
max_val = max([cm.objArray(1:range,3); max_val]);
load('result/synth_dppca/dppca_N01_G001_E10.mat');
min_val = min([cm.objArray(1:range,1); min_val]);
max_val = max([cm.objArray(1:range,1); max_val]);

load('result/synth_dppca/dppca_N10_G002_E10.mat');
if show_in_log_scale
%    cm.objArray(1:range,11) = cm.objArray(1:range,11) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,11), ':r','LineWidth',2.5,'MarkerSize',5);

load('result/synth_dppca/dppca_N08_G002_E10.mat');
if show_in_log_scale
%    cm.objArray(1:range,9) = cm.objArray(1:range,9) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,9), '-.g','LineWidth',2,'MarkerSize',7);

load('result/synth_dppca/dppca_N05_G002_E10.mat');
if show_in_log_scale
%    cm.objArray(1:range,6) = cm.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,6), '--b','LineWidth',2,'MarkerSize',7);

load('result/synth_dppca/dppca_N02_G002_E10.mat');
if show_in_log_scale
%    cm.objArray(1:range,3) = cm.objArray(1:range,3) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,3), '-.m','LineWidth',2,'MarkerSize',7);
               
load('result/synth_dppca/dppca_N01_G001_E10.mat');
if show_in_log_scale
%    cm.objArray(1:range,1) = cm.objArray(1:range,1) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,1), '-k','LineWidth',2,'MarkerSize',7);

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
load('result/synth_dppca/dppca_N05_G002_E10.mat');
min_val = min([cm.objArray(1:range,6); min_val]);
max_val = max([cm.objArray(1:range,6); max_val]);
load('result/synth_dppca/dppca_N05_G003_E10.mat');
min_val = min([cm.objArray(1:range,6); min_val]);
max_val = max([cm.objArray(1:range,6); max_val]);
load('result/synth_dppca/dppca_N05_G004_E10.mat');
min_val = min([cm.objArray(1:range,6); min_val]);
max_val = max([cm.objArray(1:range,6); max_val]);
load('result/synth_dppca/dppca_N01_G001_E10.mat');
min_val = min([cm.objArray(1:range,1); min_val]);
max_val = max([cm.objArray(1:range,1); max_val]);

load('result/synth_dppca/dppca_N05_G002_E10.mat');
if show_in_log_scale
%    cm.objArray(1:range,6) = cm.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,6), ':r','LineWidth',2.5,'MarkerSize',7);

load('result/synth_dppca/dppca_N05_G003_E10.mat');
if show_in_log_scale
%    cm.objArray(1:range,6) = cm.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,6), '-.g','LineWidth',2,'MarkerSize',7);

load('result/synth_dppca/dppca_N05_G004_E10.mat');
if show_in_log_scale
%    cm.objArray(1:range,6) = cm.objArray(1:range,6) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,6), '--b','LineWidth',2,'MarkerSize',7);

load('result/synth_dppca/dppca_N01_G001_E10.mat');
if show_in_log_scale
%    cm.objArray(1:range,1) = cm.objArray(1:range,1) - repmat(min_val, [range 1]);
end
plot(cm.objArray(1:range,1), '-k','LineWidth',2,'MarkerSize',7);

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
