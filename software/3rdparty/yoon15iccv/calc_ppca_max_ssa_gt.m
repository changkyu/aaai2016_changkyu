function angle = calc_ppca_max_ssa_gt( W, W_gt, NV )

assert(~iscell(W_gt), 'W_gt must be a matrix!');

angle = -Inf;

if iscell(W)
    % distributed case
    for j = 1 : NV
        cur_angle = subspace( W{j}, W_gt );
        if cur_angle > angle
            angle = cur_angle;
        end
    end
else
    % centralized case
    angle = subspace( W, W_gt );
end

end