%--------------------------------------------------------------------------
% Synthetic data result plot (Missing Values)
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.06.01 (last modified on 2012/06/01)
%--------------------------------------------------------------------------
fontsize = 18;
range = 10^3;

%% Figure 4. Different rate of missing values (MAR)
figure;
h = axes;
hold on;
hl = xlabel('Missing Rates (%)'); 
set(hl, 'FontSize',fontsize);
hl = ylabel('Average Root Mean Squared error'); 
set(hl, 'FontSize',fontsize);

RMS_all = [];
for idk = 1:5
    load(sprintf('result/synth_dppca_m/data_R%d.mat', idk));
    X = cm;
    MissIDX = cm_init;
    
    RMS_cs = [];
    RMS_ds = [];
    
    for idr = 1:5%20
        load(sprintf('result/synth_dppca_m/cppca_R%d_i%02d.mat', idk, idr));
        cm_c = cm;
        Xbar_c = miss_value_reconst_c(cm_c.W, cm_c.EZ, cm_c.MU);

        load(sprintf('result/synth_dppca_m/dppca_R%d_i%02d.mat', idk, idr));
        cm_d = cm;
        Xbar_d = miss_value_reconst_d(cm_d.W, cm_d.EZ, cm_d.MU);

        RMS_c = sqrt(mse(X - Xbar_c));
        RMS_d = sqrt(mse(X - Xbar_d));
        
        RMS_cs = [RMS_cs RMS_c];
        RMS_ds = [RMS_ds RMS_d];
    end
    
    RMS_c = mean(RMS_cs);
    RMS_d = mean(RMS_ds);
    
    RMS_all = [RMS_all; RMS_c, RMS_d];
end
bar(RMS_all);

ylim([0.4, 0.48]);
set(gca, 'XLim', [.5 5.5]);
set(gca, 'XTick', [1:5]);
set(gca, 'XTickLabel',{'1'; '5'; '10'; '20'; '30'});

set(h,'FontSize',fontsize);
legend('PPCA','D-PPCA');
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
for idk = 1:5
    load(sprintf('result/synth_dppca_m/data_NR%d.mat', idk));
    X = cm;
    MissIDX = cm_init;
    
    RMS_cs = [];
    RMS_ds = [];
    
    for idr = 1:5%20
        load(sprintf('result/synth_dppca_m/cppca_NR%d_i%02d.mat', idk, idr));
        cm_c = cm;
        Xbar_c = miss_value_reconst_c(cm_c.W, cm_c.EZ, cm_c.MU);

        load(sprintf('result/synth_dppca_m/dppca_NR%d_i%02d.mat', idk, idr));
        cm_d = cm;
        Xbar_d = miss_value_reconst_d(cm_d.W, cm_d.EZ, cm_d.MU);

        RMS_c = sqrt(mse(X - Xbar_c));
        RMS_d = sqrt(mse(X - Xbar_d));
        
        RMS_cs = [RMS_cs RMS_c];
        RMS_ds = [RMS_ds RMS_d];
    end
    
    RMS_c = mean(RMS_cs);
    RMS_d = mean(RMS_ds);
    
    RMS_all = [RMS_all; RMS_c, RMS_d];
end
bar(RMS_all);

ylim([0.4, 0.48]);
set(gca, 'XLim', [.5 5.5]);
set(gca, 'XTick', [1:5]);
set(gca, 'XTickLabel',{'1'; '5'; '10'; '20'; '30'});

set(h,'FontSize',fontsize);
legend('PPCA','D-PPCA');
set(gca,'YGrid','on');
