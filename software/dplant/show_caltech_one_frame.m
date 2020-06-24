addpath('../3rdparty/yoon15iccv');
addpath('../fadmm1');
addpath('../fadmm1/altmany-export_fig');

obj_names = {'BallSander','BoxStuff','Rooster','Standing','StorageBin'};

resize_factor = 1;

GTs = load('../../../../data/caltech_turntable/data.mat');
GTs = {GTs.BallSander, GTs.BoxStuff, GTs.Rooster, GTs.Standing, GTs.StorageBin};
CENTs = cell(length(GTs),1);
for idm = 1:length(GTs)
    CENTs{idm} = mean(GTs{idm}, 2);
    GTs{idm} = round(GTs{idm}*resize_factor);
end
PTs = {size(GTs{1},2), size(GTs{2},2), size(GTs{3},2), size(GTs{4},2), size(GTs{5},2),};

for idm = 1 : 5
    clf;

    img_frame = imread(sprintf('../../../../data/caltech_turntable/%s/img_0001.png', obj_names{idm}));
    [h, w, d] = size(img_frame);

    %% do SfM example
    [sfm_M, sfm_S, sfm_b, U3, W3, V3] = sfm_affine(GTs{idm});
    
    % Compute reprojected points
    Est3D = sfm_S';
    EstCam = sfm_M';

    % Make homogeneous coordinates
    Est3Dh = [Est3D ones(size(Est3D, 1), 1)]';
    EstCamh = permute(reshape(EstCam', [2, size(EstCam, 2)/2, 3]), [1 3 2]);
    EstCamh(3, 3, :) = 1; EstCamh(3, 4, :) = 0; 
    EstCamh(1:2, 4, :) = reshape(CENTs{idm}, 2, 1, size(EstCam, 2)/2);
    
    % Compute world-coordinate projected 3D points
    frames = size(GTs{idm},1)/2;
    Est3Dhp = zeros(3, PTs{idm}, frames);
    for idf = 1:frames
        Est3Dhp(:,:,idf) = EstCamh(:,:,idf) * Est3Dh;
    end    
    Est3Dhp(1:2,:,:) = Est3Dhp(1:2,:,:) * resize_factor;    

    %% visualize
    subplot(1,2,1);
    imshow(img_frame);
    hold on;
    scatter(GTs{idm}(1,:),GTs{idm}(2,:),'g+');
    
    subplot(1,2,2);
    idf = 15; % frame to show
    plot3(Est3Dhp(1,:,idf)-(w/2), -Est3Dhp(3,:,idf), -Est3Dhp(2,:,idf)+(h/2), ...
            'r+', 'MarkerSize', 7);
    xlabel('x');
    ylabel('z');
    zlabel('y');
    grid on;
    set(gca,'ztick',[]);            
    set(gca,'ytick',[]);            
    set(gca,'xtick',[]);    
    
    tightfig;
    set(gcf, 'Position', [0 0 1280 500]);
    set(gcf,'Color','white');
    export_fig(gcf, ...
        sprintf('../../../img/supp/%s.eps', obj_names{idm}), ...
        '-eps');
end

close gcf;