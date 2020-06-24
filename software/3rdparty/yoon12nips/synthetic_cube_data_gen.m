%% Synthetic cube data generation
%--------------------------------------------------------------------------
clear;

% Noise levels
noises = [1e-5, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1];
% We generate this number of independent trials
reps = 20; 
model = 'cube';

% create directory if not available
if ~exist(sprintf('data/%s', model), 'dir')
    mkdir(sprintf('data/%s',model));
end

for idr = 1:reps

    for k = 1:length(noises)
        % clear variables...
        clearvars -except noises model reps idr k;

        % Create data of viewing a 3D unit cube from a set of 
        % cameras distributed on a circle
        datasetname = model;

        % Unit cube
        x = [0  1  0  0  1  1  0  1];  x = x-.5;
        y = [0  0  1  0  1  0  1  1];  y = y-.5;
        z = [0  0  0  1  0  1  1  1];  z = z-.5;

        % This plot index will be used to draw the cube.  
        % e.g., plot3(x(plotIdx),y(plotIdx),z(plotIdx)).
        plotIdx = [1 2 5 3 1 4 6 8 7 4 6 2 5 8 7 3];
        
        GT = [x' y' z'];

        [m,n] = size(x);

        %% Create a set of stationary cameras on a circle
        Ncam = 5;  % number of cameras
        thetaCam = linspace(0,pi/2,Ncam);  % cameras are spaced according to this

        % Cameras in (x,y) place
        xc = cos(thetaCam);
        yc = sin(thetaCam);
        zc = zeros(size(thetaCam));

        % Rotate the camera arc outside of (x,y) plane
        axisCircle  = [1  0  0];
        angleCircle = pi/6;
        M  = makehgtform('axisrotate',axisCircle,angleCircle);
        Xc = M*[ xc(:) yc(:) zc(:) ones(Ncam,1) ]';

        % Transform to spherical for camera
        [Az,El] = cart2sph(Xc(1,:),Xc(2,:),Xc(3,:));

        % Due to weird difference between cart2sph and view
        % functions' azimuth/elevation convention, we need to
        % rotate 90 degree more for azimuth.
        %
        % References
        % http://www.mathworks.com/help/matlab/ref/cart2sph.html?searchHighlight=cart2sph
        % http://www.mathworks.com/help/matlab/ref/view.html?searchHighlight=view    
        Az = Az + pi/2;

        %% Create cube views (2D vertex coordinates) from each camera

        % Rotate the cube around z-axis, then compute views from each Ncam cameras
        % for each cube rotation.

        Ncube = 5; % # of frames
        thetaCube = linspace(0,pi/6,Ncube);

        for r=1:Ncube,

            fprintf('Generating %s data for noise level %f, frame# %d (Run: %d)...\n', ...
                model, noises(k), r, idr);

            M = makehgtform('axisrotate',[0 0 1],thetaCube(r));

            Xrc = M*[ x;y;z;ones(1,m*n)];

            % Plot the camera-object scene
            clf
            hold on;

            for i=1:Ncam

                A = viewmtx(Az(i)/pi*180,El(i)/pi*180);

                Xv = A*Xrc;

                % Image coordinates + noise
                X{i,r} = Xv(1:2,:) + rand(2,size(Xv,2))*noises(1,k);

                subplot('Position',[0.15*i,0.85,0.1,0.1]);
                line(X{i,r}(1,plotIdx),X{i,r}(2,plotIdx));
                title(sprintf('Camera %d',i));

                GTX{i,r} = Xv(1:2,:);
            end

            subplot('Position',[0.05,0.05,0.75,0.75]);
            plot3(Xrc(1,plotIdx),Xrc(2,plotIdx),Xrc(3,plotIdx),'linewidth',2);

            hold on;

            % Draw cameras
            plot3(5*Xc(1,:),5*Xc(2,:),5*Xc(3,:),'or','linewidth',2);
            Vc = .99*Xc - Xc;
            quiver3(5*Xc(1,:),5*Xc(2,:),5*Xc(3,:),5*Vc(1,:),5*Vc(2,:),5*Vc(3,:),0.5);
            view(3);
            grid on;
            axis equal

            drawnow
        end

        clf;

        % make measurement matrix
        X = X';
        [r, c] = size(X);
        X = reshape(X, [r * c, 1]);
        X = cell2mat(X);
        measurement_matrix = X;

        % GTX is ground truth coordinates without noise
        GTX = GTX';
        [r, c] = size(GTX);
        GTX = reshape(GTX, [r * c, 1]);
        GTX = cell2mat(GTX);        

        save(sprintf('data/%s/%03d_%03d.mat', model, k, idr), ...
            'measurement_matrix', 'plotIdx', 'GT', 'GTX', 'datasetname');

    end

    close gcf;
end

%% Synthetic Cube data generation (with Missing Values)
%--------------------------------------------------------------------------
clear;

noise = 1e-5;

% Create data of viewing a 3D unit cube from a set of cameras distributed
% on a circle
model = 'cube';
datasetname = model;

switch( model )

  case 'cube'
    % Unit cube
    x = [0  1  0  0  1  1  0  1];  x = x-.5;
    y = [0  0  1  0  1  0  1  1];  y = y-.5;
    z = [0  0  0  1  0  1  1  1];  z = z-.5;

    % This plot index will be used to draw the cube.  E.g.,
    % plot3(x(plotIdx),y(plotIdx),z(plotIdx)).
    plotIdx = [1 2 5 3 1 4 6 8 7 4 6 2 5 8 7 3];

end

GT = [x' y' z'];

[m,n] = size(x);

%% Create a set of stationary cameras on a circle
Ncam = 5;  % number of cameras
thetaCam = linspace(0,pi/2,Ncam);  % cameras are spaced according to this

% Cameras in (x,y) place
xc = cos(thetaCam);
yc = sin(thetaCam);
zc = zeros(size(thetaCam));

% Rotate the camera arc outside of (x,y) plane
axisCircle  = [1  0  0];
angleCircle = pi/6;
M  = makehgtform('axisrotate',axisCircle,angleCircle);
Xc = M*[ xc(:) yc(:) zc(:) ones(Ncam,1) ]';

% Transform to spherical for camera
[Az,El] = cart2sph(Xc(1,:),Xc(2,:),Xc(3,:));


%% Create cube views (2D vertex coordinates) from each camera

% Rotate the cube around z-axis, then compute views from each Ncam cameras
% for each cube rotation.

Ncube = 5; % # of frames
thetaCube = linspace(0,pi/6,Ncube);

for r=1:Ncube,

    fprintf('Generating %s data for frame# %d...\n',...
        model,  r);

    M = makehgtform('axisrotate',[0 0 1],thetaCube(r));

    Xrc = M*[ x;y;z;ones(1,m*n)];

    % Plot the camera-object scene
    clf
    hold on;

    for i=1:Ncam

        A = viewmtx(Az(i)/pi*180,El(i)/pi*180);

        Xv = A*Xrc;

        % Image coordinates + noise
        X{i,r} = Xv(1:2,:) + rand(2,size(Xv,2))*noise;

        subplot('Position',[0.15*i,0.85,0.1,0.1]);
        line(X{i,r}(1,plotIdx),X{i,r}(2,plotIdx));
        title(sprintf('Camera %d',i));
 
        GTX{i,r} = Xv(1:2,:);
    end

    subplot('Position',[0.05,0.05,0.75,0.75]);
    plot3(Xrc(1,plotIdx),Xrc(2,plotIdx),Xrc(3,plotIdx),'linewidth',2);

    hold on;

    % Draw cameras
    plot3(5*Xc(1,:),5*Xc(2,:),5*Xc(3,:),'or','linewidth',2);
    Vc = .99*Xc - Xc;
    quiver3(5*Xc(1,:),5*Xc(2,:),5*Xc(3,:),5*Vc(1,:),5*Vc(2,:),5*Vc(3,:),0.5);
    view(3);
    grid on;
    axis equal

    drawnow
end

% make measurement matrix
X = X';
[r, c] = size(X);
X = reshape(X, [r * c, 1]);
X = cell2mat(X);
measurement_matrix = X;

% GTX is ground truth coordinates without noise
GTX = GTX';
[r, c] = size(GTX);
GTX = reshape(GTX, [r * c, 1]);
GTX = cell2mat(GTX);        

% This is the missing index matrix manually obtained missing points (i.e. 
% occluded points from each camera viewpoint) Note that the first point is 
% always occluded, thus never observed across all 25 frames.
%          x                y
%          1 2 3 4 5 6 7 8  1 2 3 4 5 6 7 8
MissIDX = [0 1 0 0 1 1 0 1; 0 1 0 0 1 1 0 1; % Cam 1
           1 1 0 1 1 1 0 1; 1 1 0 1 1 1 0 1;
           1 1 0 1 1 1 0 1; 1 1 0 1 1 1 0 1;
           1 1 0 1 1 1 0 1; 1 1 0 1 1 1 0 1;
           1 1 0 1 1 1 0 1; 1 1 0 1 1 1 0 1; 
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1; % Cam 2
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           1 1 0 1 1 1 1 1; 1 1 0 1 1 1 1 1;
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1; % Cam 3
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1; % Cam 4
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           0 0 1 1 1 1 1 1; 0 0 1 1 1 1 1 1; % Cam 5
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1;
           0 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 1];

% flip the indicator (0 = non-missing, 1 = missing)
MissIDX = ~MissIDX;
       
% Note that for missing data exeperiments, independent runs or noise levels
% are not the main concern. We only control different set of missing values
% to see the exact effect of missing values.
save(sprintf('data/%s/data_missing.mat', model), ...
    'measurement_matrix', 'plotIdx', 'GT', 'GTX', 'datasetname', 'MissIDX');

close gcf;
