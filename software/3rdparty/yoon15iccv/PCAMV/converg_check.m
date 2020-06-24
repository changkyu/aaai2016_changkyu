%  CONVERG_CHECK - Check convergence criteria
%
%  convmsg = CONVERG_CHECK( opts, lc, angleA, sd_iter ) checks
%  convergence criteria and generates a message that should be
%  displayed.

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function convmsg = converg_check( opts, lc, angleA, sd_iter )

convmsg = '';

if angleA < opts.minangle
    convmsg = ...
        sprintf( [ 'Convergence achieved (angle between subpaces'...
                   ' smaller than %.2e)\n' ], opts.minangle );
    
elseif opts.earlystop && lc.prms(end) > lc.prms(end-1)
    convmsg = sprintf( 'Early stopping.\n' );
    
elseif isfield( opts, 'rmsstop' ) && ~isempty(opts.rmsstop) ...
        && length(lc.rms)-1 > opts.rmsstop(1)
    
    numiter = opts.rmsstop(1);
    abs_tol = opts.rmsstop(2);
    rel_tol = [];
    if length(opts.rmsstop) > 2
        rel_tol = opts.rmsstop(3);
    end
    rms1 = lc.rms(end-numiter);
    rms2 = lc.rms(end);
    
    if abs(rms1-rms2) < abs_tol || ...
       ( length(opts.rmsstop) > 2 && abs( (rms1-rms2)/rms2 ) < rel_tol )
        
        convmsg = ...
            sprintf( 'Stop: RMS does not change much for %d iterations.\n',...
                     numiter );
    end

elseif isfield( opts, 'cfstop' ) && ~isempty(opts.cfstop) && ...
        length(lc.cost)-1 > opts.cfstop(1)
    
    numiter = opts.cfstop(1);
    abs_tol = opts.cfstop(2);
    rel_tol = [];
    if length(opts.cfstop) > 2
        rel_tol = opts.cfstop(3);
    end
    cost1 = lc.cost(end-numiter);
    cost2 = lc.cost(end);
    
    if abs(cost1-cost2) < abs_tol || ...
       ( length(opts.rmsstop) > 2 && abs( (cost1-cost2)/cost2 ) < rel_tol )
        
        convmsg = ...
            sprintf( 'Stop: Cost does not change much for %d iterations.\n',...
                     numiter );
    end
elseif nargin >4 && sd_iter == 40
    convmsg = ...
        sprintf(...
            [ 'Slowing-down stop. ' ...
              'You may continue by changing the gradient type.\n' ] );

end

