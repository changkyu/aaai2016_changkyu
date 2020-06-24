function [cm1, cm2] = expr_run_cppca(X, M, THRESHc, objfreq_c, m_init_c, ZeroMean)

%% ------------------------------------------------------------------------
disp('*** Centralized Setting ***');

disp('* PPCA (EM) *');
cm1 = cppca_em(X, M, ...
    'InitModel', m_init_c, 'Threshold', THRESHc, 'ShowObjPer', objfreq_c, ...
    'ZeroMean', ZeroMean, 'MaxIter', 1000);

% disp('* BPCA *');
% cm2 = cbpca(X, M, ...
%     'InitModel', m_init_c, 'Threshold', THRESHc, 'ShowObjPer', objfreq_c, ...
%     'ZeroMean', ZeroMean);

cm2 = cm1;

end
