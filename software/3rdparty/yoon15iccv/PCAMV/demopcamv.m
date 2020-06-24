clear all
randn('state',4)
rand('state',4)

clc;
disp('This demo gives a simple example of using the toolbox for PCA')
disp('with missing data. ')
disp(' ')
disp('First, we generate synthetic data according to the model')
disp(' ')
disp('        X(i,j) = A(i,:)*S(:,j) + noise')
disp(' ')
disp('with two-dimensional S(:,j) and 5-dimensional X(:,j). Matrix A' )
disp('is sampled from a Gaussian distribution, S is a simple dynamical' )
disp('process and Gaussian noise is added.')
disp(' ')
disp('We remove 20% of the observations from the data matrix X.')
disp(' ')
disp('Press any key to see plots of S and X.')
pause;

% Generate the matrix of inputs x and targets t.

n1 = 5;
n2 = 500;
ndata = 2000;
nprobe = 250;
ncomp = 2;

Ar = orth(randn(n1,ncomp))*diag(ncomp:-1:1);
T = 1:n2;
Sr = [ exp(-T/150).*cos( 2*pi*T/50 )
       exp(-T/150).*sin( 2*pi*T/50 ) ];
% Normalizing to zero mean and unit variance
Sr = ( Sr - repmat( mean(Sr,2), 1, n2 ) );
Sr = Sr ./ repmat( sqrt( mean( Sr.^2, 2 ) ), 1, n2 );
Xr = Ar * Sr;
Xrnoise = Xr + 0.2 * randn(n1,n2);

h1 = tsplot(Sr);

rp = randperm( n1*n2 );

if 1
    % Missing values replaced with NaNs
    X = Xrnoise;
    Xprobe = Xrnoise;
    X( rp(1:end-ndata) ) = NaN;
    Xprobe( rp(nprobe+1:end) ) = NaN;
else
    % Sparse matrices with only observed values
    data = Xrnoise( rp(end-ndata+1:end) )';
    data( data == 0 ) = eps;
    pdata = Xrnoise( rp(1:nprobe) )';
    pdata( pdata == 0 ) = eps;
    
    [IX,JX] = ind2sub( [n1,n2], rp(end-ndata+1:end) );
    [pIX,pJX] = ind2sub( [n1,n2], rp(1:nprobe) );
    X = sparse( IX, JX, data, n1, n2 );
    Xprobe = sparse( pIX, pJX, pdata, n1, n2 );
end

h2 = tsplot(X);

disp(' ')
disp('Press any key to continue')
pause; clc;

close( [ get( h1(1), 'Parent' ), get( h2(1), 'Parent' ) ] )

disp('Now we train a variational Bayesian PCA model with two components.')
disp(' ')
disp('Press any key to continue')
pause;
echo on

opts = struct( 'maxiters', 30,...
               'algorithm', 'vb',...
               'xprobe', Xprobe,...
               'uniquesv', 0,...
               'cfstop', [ 100 0 0 ],...
               'minangle', 0 );
[ A, S, Mu, V, cv, hp, lc ] = pca_full( X, ncomp, opts );

echo off

disp('Learning is finished.')
disp(' ')
disp('Press any key to continue')
pause; clc;

disp('The next plot displays the found principal components (green)')
disp('and confidence regions (blue shade) as well as the original' )
disp('components (red).')
disp(' ')
disp('Press any key to see the plot')
pause;

%Sr1 = pinv(A) * Xr;
Sr1 = round( Sr*S'/n2 )*Sr;
Sv = cov_f2d(cv.S,cv.Isv);

hax = tsplot( S, 'g' );
Sbnd = sqrt(Sv)*3;
addebars( Sbnd );
addtsplot( Sr1, 'r' )

axes(hax(1))
legend( 'error bars', 'found', 'true' )

disp(' ')
disp('Press any key to continue')
pause; clc;

disp('The next plot displays the principal subspace in which the found')
disp('components are marked with dots and the original components are')
disp('shown with crosses. Get the estimated confidence regions by moving')
disp('the mouse over the dots.')
disp(' ')
disp('Press any key to see the plot')
pause;

subspace2d( S(:,1:200), cv.S, Sr1(:,1:200) )

disp(' ')
disp('Press any key to continue')
pause; clc;

disp('Finally we can plot the reconstructions of the missing data (green)')
disp('and the corresponding confidence intervals (blue shade). The original')
disp('noiseless data are shown in red.')
disp(' ')
disp('Press any key to end.')
pause;

Xrec = repmat(Mu,1,n2) + A*S;

for i = 1:size(X,1)
    for j = 1:size(X,2)
        Vr(i,j) = A(i,:) * cv.S{j} * A(i,:)' + ...
                  S(:,j)' * cv.A{i} * S(:,j) + ...
                  sum( sum( cv.S{j} .* cv.A{i} ) ) + cv.Mu(i);
    end
end
%Vr = Vr + V;

hax = tsplot( Xrec, 'g' );
Xbnd = sqrt(Vr)*3;
addebars( Xbnd );
addtsplot( Xr, 'r' )

axes(hax(1))
legend( 'error bars', 'reconstruction', 'true' )
