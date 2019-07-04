load experimental_measurements.mat

nx = 8;
ny = 12;
sb = (nx+1)*(ny+1);
nb = 4;

f_gt = 15/0.03;%15 mm, pixels are 30 um 
u0_gt = 160;
board_size = 25/0.03;%25 mm, pixels are 30 um

measurements = cell(1,nb);
for b = 1:nb
    measurements{b} = xs_SWIR(:,(b-1)*sb+(1:sb));
end
Xs_HS = zeros(3,sb);
idx = 1;
for y = 1:ny+1
for x = 1:nx+1
    Xs_HS(1,idx) = y * board_size;
    Xs_HS(2,idx) = x * board_size;
    idx = idx + 1;
end
end

%visualize to make sure we have correspondences right
% figure,axes,hold all
% Xs_HS = Xs_HS / 100;
% scatter(measurements{1}(1,:),measurements{1}(2,:))
% scatter(Xs_HS(1,:),Xs_HS(2,:))
% plot([Xs_HS(1,1),measurements{1}(1,1)],[Xs_HS(2,1),measurements{1}(2,1)],'r')
% plot([Xs_HS(1,9),measurements{1}(1,9)],[Xs_HS(2,9),measurements{1}(2,9)],'g')
% plot([Xs_HS(1,9*13),measurements{1}(1,9*13)],[Xs_HS(2,9*13),measurements{1}(2,9*13)],'c')
% plot([Xs_HS(1,9*12+1),measurements{1}(1,9*12+1)],[Xs_HS(2,9*12+1),measurements{1}(2,9*12+1)],'m')
% Xs_HS = Xs_HS * 100;

%% actually calibration the camera
    [camera_est,R_est,t_est] = calibrate_linear_camera_donne(measurements,Xs_HS,f_gt,u0_gt);
    [camera_ba,R_ba,t_ba,final_MSE] = refine_linear_camera(measurements,Xs_HS,camera_est,R_est,t_est,false);