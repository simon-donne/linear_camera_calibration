%calibrate_synthetic_example
%   Performs calibration of a random camera on a random scene, and
%   evaluates the result.
%
%   Original code by Simon Donné, January 2017
function [] = example_calibrate_synthetic(noiselevel)

    if nargin < 1
        noiselevel = 0.5;%pixels
    end

    %% the scene
    nx = 9;
    ny = 9;
    nb = 10;
    rb = 50;%pixels
    [scene,base_board,R_gt,t_gt] = create_synthetic_scene(nx,ny,rb,nb,0.5,0.0);
    
    
    %% the camera
    camera_gt = create_linear_camera(1000,500,50);
    f_gt = camera_gt(1,1);
    u0_gt = camera_gt(1,3);
    s_gt = camera_gt(2,2);
    
    %% the measurements
    measured = cell(1,nb);
    for b = 1:nb
        projected = linear_camera_projection(camera_gt,scene{b});
        measured{b} = projected + randn(size(projected))*noiselevel;
    end
    
    %% calibrate the camera
    [camera_est,R_est,t_est] = calibrate_linear_camera_donne(measured,base_board);
    f_est = camera_est(1,1);
    u0_est = camera_est(1,3);
    s_est = camera_est(2,2);
    
    [camera_ba,R_ba,t_ba,final_MSE] = refine_linear_camera(measured,base_board,camera_est,R_est,t_est);
    f_ba = camera_ba(1,1);
    u0_ba = camera_ba(1,3);
    s_ba = camera_ba(2,2);
    
    %% print the results
    fprintf('Focal distance groundtruth: %4.4f\t\tEstimated: %4.4f\t\tAfter BA: %4.4f\n',f_gt,f_est,f_ba);
    fprintf('Optical center groundtruth: %4.4f\t\tEstimated: %4.4f\t\tAfter BA: %4.4f\n',u0_gt,u0_est,u0_ba);
    fprintf('Scanning speed groundtruth: %4.4f\t\tEstimated: %4.4f\t\tAfter BA: %4.4f\n',s_gt,s_est,s_ba);
    fprintf('MSE after bundle adjustment: %4.4f\n',final_MSE);
    
end