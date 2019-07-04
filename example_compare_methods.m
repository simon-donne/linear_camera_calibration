%calibrate_synthetic_example
%
%   Generates the comparison plots from 
%       "Robustifying Plane-base Calibration for Linear Cameras"
%
%   Original code by Simon Donné, January 2017

    %% the number of iterations per simulation
    simulations = 500;

    %% the scene parameters, when the simulation does not override them
    nx = 9;
    ny = 9;
    nb = 10;
    rb = 50;%pixels
    
    %% the camera (always the same)
    camera_gt = create_linear_camera(1000,500,50);
    f_gt = camera_gt(1,1);
    u0_gt = camera_gt(1,3);
    s_gt = camera_gt(2,2);
    
    %% varying the noise level
    if ~exist('sigmas','var')
    sigmas = 0.2:0.2:2.0;
    nrsigmas = numel(sigmas);
    focal_errors_drareni_nl = zeros(nrsigmas,simulations);
    center_errors_drareni_nl = zeros(nrsigmas,simulations);
    focal_errors_donne_nl = zeros(nrsigmas,simulations);
    center_errors_donne_nl = zeros(nrsigmas,simulations);
    valids_drareni_nl = zeros(nrsigmas,simulations);
    for sigma = 1:nrsigmas
        noiselevel = sigmas(sigma);
        for simulation = 1:simulations
            clc
            fprintf('Simulation %d/%d of noiselevel %d/%d\n',simulation,simulations,sigma,nrsigmas);
            
            %% getting the measurements
            [scene,base_board,R_gt,t_gt] = create_synthetic_scene(nx,ny,rb,nb);
            measured = cell(1,nb);
            for b = 1:nb
                projected = linear_camera_projection(camera_gt,scene{b});
                measured{b} = projected + randn(size(projected))*noiselevel;
            end

            %% calibration with the existing method
            [camera_est,R_est,t_est] = calibrate_linear_camera_drareni(measured,base_board);
            [camera_ba,R_ba,t_ba,drareni_MSE] = refine_linear_camera(measured,base_board,camera_est,R_est,t_est);
            valids_drareni_nl(sigma,simulation) = drareni_MSE < 4*noiselevel^2;
            focal_errors_drareni_nl(sigma,simulation) = abs(camera_ba(1,1)-f_gt);
            center_errors_drareni_nl(sigma,simulation) = abs(camera_ba(1,3)-u0_gt);
            
            %% calibration with the proposed method
            [camera_est,R_est,t_est] = calibrate_linear_camera_donne(measured,base_board);
            [camera_ba,R_ba,t_ba] = refine_linear_camera(measured,base_board,camera_est,R_est,t_est);            
            focal_errors_donne_nl(sigma,simulation) = abs(camera_ba(1,1)-f_gt);
            center_errors_donne_nl(sigma,simulation) = abs(camera_ba(1,3)-u0_gt);
        end
    end
    end
    
    %it does seem like there are still painful outliers around, for drareni
    valids_drareni_nl = valids_drareni_nl & (focal_errors_drareni_nl < 100) & (center_errors_drareni_nl < 100);
    focal_errors_drareni_nl(valids_drareni_nl~=1) = 0;
    center_errors_drareni_nl(valids_drareni_nl~=1) = 0;
    
    figure,axes,hold all
    title('Varying the noiselevel of the measurements')
    xlabel('Sigma (pixels)')
    ylabel('Absolute error (pixels)')
    plot(sigmas,sum(focal_errors_drareni_nl,2)./sum(valids_drareni_nl,2),'go-')
    plot(sigmas,mean(focal_errors_donne_nl,2),'bo-')
    plot(sigmas,sum(center_errors_drareni_nl,2)./sum(valids_drareni_nl,2),'g+-')
    plot(sigmas,mean(center_errors_donne_nl,2),'b+-')
    legend({'Focal distance (Drareni)','Focal distance (Donne)','Optical Center (Drareni)','Optical center (Donne)'});
    hold off
    drawnow
    
    %% varying the number of planes
    if ~exist('plane_counts','var')
    plane_counts = 6:2:24;
    nrplanecounts = numel(plane_counts);
    focal_errors_drareni_pc = zeros(nrplanecounts,simulations);
    center_errors_drareni_pc = zeros(nrplanecounts,simulations);
    focal_errors_donne_pc = zeros(nrplanecounts,simulations);
    center_errors_donne_pc = zeros(nrplanecounts,simulations);
    valids_drareni_pc = zeros(nrplanecounts,simulations);
    for planes = 1:nrplanecounts
        planecount = plane_counts(planes);
        noiselevel = 0.5;
        
        for simulation = 1:simulations
            clc
            fprintf('Simulation %d/%d of planecount %d/%d\n',simulation,simulations,planes,nrplanecounts);
            
            %% getting the measurements
            [scene,base_board,R_gt,t_gt] = create_synthetic_scene(nx,ny,rb,planecount);
            measured = cell(1,planecount);
            for b = 1:planecount
                projected = linear_camera_projection(camera_gt,scene{b});
                measured{b} = projected + randn(size(projected))*noiselevel;
            end

            %% calibration with the existing method
            [camera_est,R_est,t_est] = calibrate_linear_camera_drareni(measured,base_board);
            [camera_ba,R_ba,t_ba,drareni_MSE] = refine_linear_camera(measured,base_board,camera_est,R_est,t_est);
            valids_drareni_pc(planes,simulation) = drareni_MSE < 4*noiselevel^2;
            focal_errors_drareni_pc(planes,simulation) = abs(camera_ba(1,1)-f_gt);
            center_errors_drareni_pc(planes,simulation) = abs(camera_ba(1,3)-u0_gt);
            
            %% calibration with the proposed method
            [camera_est,R_est,t_est] = calibrate_linear_camera_donne(measured,base_board);
            [camera_ba,R_ba,t_ba] = refine_linear_camera(measured,base_board,camera_est,R_est,t_est);
            
            focal_errors_donne_pc(planes,simulation) = abs(camera_ba(1,1)-f_gt);
            center_errors_donne_pc(planes,simulation) = abs(camera_ba(1,3)-u0_gt);
        end
    end
    end
    
    %it does seem like there are still painful outliers around, for drareni
    valids_drareni_pc = valids_drareni_pc & (focal_errors_drareni_pc < 100) & (center_errors_drareni_pc < 100);
    focal_errors_drareni_pc(valids_drareni_pc~=1) = 0;
    center_errors_drareni_pc(valids_drareni_pc~=1) = 0;
    
    figure,axes,hold all
    title('Varying the number of planes')
    xlabel('Number of planes')
    ylabel('Absolute error (pixels)')
    plot(plane_counts,sum(focal_errors_drareni_pc,2)./sum(valids_drareni_pc,2),'go-')
    plot(plane_counts,mean(focal_errors_donne_pc,2),'bo-')
    plot(plane_counts,sum(center_errors_drareni_pc,2)./sum(valids_drareni_pc,2),'g+-')
    plot(plane_counts,mean(center_errors_donne_pc,2),'b+-')
    legend({'Focal distance (Drareni)','Focal distance (Donne)','Optical Center (Drareni)','Optical center (Donne)'});
    hold off
    drawnow
    
    save compare_results.mat

%% output for gnuplot
% fh = fopen('synthetic_results_noise.txt','w');
% fprintf(fh,'Sigma\tFocal distance (Drareni)\tFocal distance (Donne)\tOptical center (Drareni)\tOptical center (Donne)\n');
% fdr = sum(focal_errors_drareni_nl,2)./sum(valids_drareni_nl,2);
% fdo = mean(focal_errors_donne_nl,2);
% odr = sum(center_errors_drareni_nl,2)./sum(valids_drareni_nl,2);
% odo = mean(center_errors_donne_nl,2);
% for i = 1:numel(sigmas)
%     fprintf(fh,'%f\t%f\t%f\t%f\t%f\n',sigmas(i),fdr(i),fdo(i),odr(i),odo(i));
% end
% fclose(fh);
% fh = fopen('synthetic_results_planes.txt','w');
% fprintf(fh,'Number of planes\tFocal distance (Drareni)\tFocal distance (Donne)\tOptical center (Drareni)\tOptical center (Donne)\n');
% fdr = sum(focal_errors_drareni_pc,2)./sum(valids_drareni_pc,2);
% fdo = mean(focal_errors_donne_pc,2);
% odr = sum(center_errors_drareni_pc,2)./sum(valids_drareni_pc,2);
% odo = mean(center_errors_donne_pc,2);
% for i = 1:numel(plane_counts)
%     fprintf(fh,'%f\t%f\t%f\t%f\t%f\n',plane_counts(i),fdr(i),fdo(i),odr(i),odo(i));
% end
% fclose(fh);