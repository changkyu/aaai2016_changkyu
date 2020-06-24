function angle = calc_ppca_max_ssa( W_old, W_new )

angle = -Inf;

if iscell(W_old)
    % distributed case
    assert(iscell(W_new), 'Old W is distributed but new W is centralized!');
        
    NV_old = length(W_old);
    NV = length(W_new);

    assert(NV_old == NV, 'Number of nodes for old and new W are different!');

    for j = 1 : NV
        cur_angle = subspace( W_old{j}, W_new{j} );
        if cur_angle > angle
            angle = cur_angle;
        end
    end
else
    % centralized case
    assert(~iscell(W_new), 'Old W is centralized but new W is distributed!');        

    angle = subspace( W_old, W_new );
end

end