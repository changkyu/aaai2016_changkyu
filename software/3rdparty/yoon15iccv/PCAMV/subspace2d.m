%  SUBSPACE2D - Visualize principal components in 2D
%
%  SUBSPACE2D( S, Sv, Sr ) visualizes principal components in S and
%  their posterior variances Sv. Each column of S should be
%  two-dimensional, Sv should contain the posterior variances of
%  S. Function COVDG is useful for obtaining Sv from the results
%  returned by PCA_FULL. Optional parameter Sr can be given to plot
%  the known values of the components.
%
%  The estimated elements are plotted with dots and the corresponding
%  true values are the crosses connected with the dots. Moving the
%  mouse over the points provides estimated error bars.
%
%  See also COVDG

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function subspace2d( S, Sv, Sr )

npts = size(S,2);

p = 0:pi/100:2*pi;

figure
set( gca, 'DrawMode', 'fast' );
hold off
if 1
    plot( S(1,:), S(2,:), '.' )
elseif 0
    varSmax = sum( diag( S*S'/npts ) )/10;
    colormap(bone)
    cmap = colormap;
    for i = 1:npts
        varS = sum(diag(Sv{i}));
        h = plot( S(1,i), S(2,i), 'k.', 'Markersize', 10 );
        set( h, 'color', ...
                cmap( ceil( min(varS,varSmax-eps)...
                             /varSmax*size(cmap,1) ), : ) )
        hold on
    end
else
    p = 0:pi/100:2*pi;
    for i = 1:npts
        [Vc,Dc] = eig(Sv{i});
        lv = 0.5*Vc*sqrt(Dc)*[cos(p);sin(p)];
        plot( S(1,i)+lv(1,:), S(2,i)+lv(2,:), 'b' );
        hold on
    end
end
hold on
ud.S = S;
ud.Sv = Sv;

if nargin > 2
    plot( Sr(1,:), Sr(2,:), 'x' )
    for i = 1:npts
        plot( [ S(1,i) Sr(1,i) ], [ S(2,i) Sr(2,i) ], ':' )
    end
    ud.Sr = Sr;
else
    ud.Sr = [];
end

xl = get( gca, 'xlim' );
yl = get( gca, 'ylim' );
set( gca, 'xlim', xl, 'ylim', yl, 'tag', 'Sub2dAxes' );

if nargin > 1 & ~isempty(Sv)
    set( gcf, 'WindowButtonMotionFcn', @Sub2dWBMot, 'UserData', ud )
    set( gcf, 'WindowButtonUpFcn', @Sub2dWBUp, 'UserData', ud )
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Sub2dWBMot(src,eventdata)

%set( src, 'units', 'pixels' )
%mcl = get( src, 'CurrentPoint' );

h = findobj( src, 'tag', 'Sub2dAxes' );
set( h, 'units', 'pixels' )
mcl = get( h, 'CurrentPoint' );

x = mcl(1,1);
y = mcl(1,2);

xl = get( h, 'xlim' );
yl = get( h, 'ylim' );

uda = get( h, 'UserData' );
if isempty(uda)
    uda.fixed = 0;
    set( h, 'UserData', uda )
end

if x >= xl(1) & x <= xl(2) & y >= yl(1) & y <= yl(2) & ~uda.fixed

    axes(h)

    ud = get( src, 'UserData' );
    S = ud.S;
    Sv = ud.Sv;
    dist = S - repmat( [x;y], 1, size(S,2) );
    dist = sum(dist.^2,1);
    [tmp,I] = min(dist);
    hold on

    p = 0:pi/100:2*pi;
    if iscell(Sv)
        [Vc,Dc] = eig(Sv{I});
        lv = 3.035*Vc*sqrt(Dc)*[cos(p);sin(p)];
    else
        lv = 3.035*diag(sqrt(Sv(:,I)))*[cos(p);sin(p)];
    end

    hbnd = findobj( h, 'tag', 'Sub2dBnd' );
    hpt = findobj( h, 'tag', 'Sub2dPt' );
    if isempty(hbnd)
        hbnd = plot( S(1,I)+lv(1,:), S(2,I)+lv(2,:), 'r' );
        set( hbnd, 'tag', 'Sub2dBnd' );
        hpt = plot( S(1,I), S(2,I), 'ro' );
        set( hpt, 'tag', 'Sub2dPt' );
    else
        set( hbnd, 'XData', S(1,I)+lv(1,:), 'YData', S(2,I)+lv(2,:) )
        set( hpt, 'XData', S(1,I), 'YData', S(2,I) );
    end

end

set( h, 'units', 'normalized' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Sub2dWBUp(src,eventdata)

set( src, 'units', 'pixels' )
%mcl = get( src, 'CurrentPoint' );

h = findobj( src, 'tag', 'Sub2dAxes' );
set( h, 'units', 'pixels' )
mcl = get( h, 'CurrentPoint' );

x = mcl(1,1);
y = mcl(1,2);

xl = get( h, 'xlim' );
yl = get( h, 'ylim' );

if x >= xl(1) & x <= xl(2) & y >= yl(1) & y <= yl(2)

    uda = get( h, 'UserData' );
    uda.fixed = ~uda.fixed;
    set( h, 'UserData', uda )
    if ~uda.fixed
        Sub2dWBMot(src,[])
    end

end

set( h, 'units', 'normalized' )
