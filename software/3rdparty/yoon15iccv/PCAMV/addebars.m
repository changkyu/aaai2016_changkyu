%  ADDEBARS - Add error bars to current figure
%
%  ADDEBARS(Yb) adds error bars in Yb to the current figure which
%  should be created by TSPLOT. Yb should contain error bars for
%  corresponding elements in Y displayed with TSPLOT.
%
%  See also TSPLOT

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function addebars( Sbars )

hax = findobj( gcf, 'type', 'axes' );

for i = 1:length(hax)
    pos = get( hax(i), 'position');
    pos2(i) = pos(2);
end
[tmp,I] = sort( -pos2 );
hax = hax(I);

for i = 1:length(hax)
    axes(hax(i))
    Add_1_EBar( Sbars(i,:) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Add_1_EBar( Sbars )

h = findobj( gca, 'type', 'line' );
if length(h) > 1
    for i = 1:length(h)
        npts(i) = length( get( h(i), 'xdata' ) );
    end
    [tmp,i] = max(npts);
    h = h(i);
end
time = get( h, 'xdata' );
S = get( h, 'ydata' );
S(isnan(S)) = 0;

hold on
hp = fill( [ time fliplr(time) ],...
            [ S+Sbars fliplr(S-Sbars) ],...
            [ 0.6 0.6 1 ], 'edgecolor', 'none' );
set( gca, 'children', flipud(get( gca, 'children' )) )
