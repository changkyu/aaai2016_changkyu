function [cm3, cm4] = expr_run_dppca(X, M, V, E, ETA, THRESHd, objfreq_d, m_init_d, ZeroMean)

%% ------------------------------------------------------------------------
disp('*** Distributed Setting ***');

disp('* D-PPCA *');
cm3 = dppca(X, M, V, E, 'Eta', ETA, ...
    'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...
    'ZeroMean', ZeroMean, 'MaxIter', 10000, 'VaryETA', false);

% disp('* D-BPCA *');
% cm4 = dbpca(X, M, V, E, 'Eta', ETA, ...
%     'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...
%     'ZeroMean', ZeroMean);

% disp('* D-PPCA-ANT *');
% cm4 = dppca_ant(X, M, V, E, 'Eta', ETA, ...
%     'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...
%     'ZeroMean', ZeroMean, 'MaxIter', 1000);

cm4 = cm3;

end
