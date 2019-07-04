%calibrate_linear_camera
%   Calibrates a linear camera given a cell array of I measurements of a
%   given 3xN planar point set.
%   Usage:
%       K = calibrate_linear_camera(measurements,points)
%
%   Implementation of the method described in 
%       "Plane-base Calibration for Linear Cameras"
%       Jamil Draréni, Sébastien Roy, Peter Sturm
%
%   Original code by Simon Donné, January 2017
function [K,R,t] = calibrate_linear_camera_drareni(uv,ab)
    I = numel(uv);
    N = size(ab,2);
    
    %% 3.3 Homography estimation
    H = cell(1,I);
    a = ab(1,:)';
    b = ab(2,:)';
    for i = 1:I
        u = uv{i}(1,:)';
        v = uv{i}(2,:)';
        %the system by drareni et al.
        H_system = [ [a,b,ones(N,1),zeros(N,6),-u.*a, -u.*b, -u];
                     [zeros(N,3),a,b,ones(N,1),a.^2,b.^2,a.*b,-a.*v,-b.* v,-v];
                     [-a.*v,-b.*v,-v,a.*u,b.*u,u,a.^2.*u,b.^2.*u,a.*b.*u,zeros(N,3)]];
                 
        [U,S,V] = svd(H_system);
        h = V(:,12)/V(12,12);
        
        %extract the homography
        H{i} = [h(1:3)',zeros(1,3);
                h(4:9)';
                h(10:11)',1,zeros(1,3)];
    end
    
%     %bonus: evaluate the reprojection error of the homographies
%     MSE = zeros(I,1);
%     for i = 1:I
%         homogeneous_reprojection = H{i} * [a,b,ones(N,1),a.^2,b.^2,a.*b]';
%         reprojection = [homogeneous_reprojection(1,:)./homogeneous_reprojection(3,:);
%                         homogeneous_reprojection(2,:)./homogeneous_reprojection(3,:)];
%         MSE(i) = mean( sum( (reprojection-uv{i}(1:2,:)).^2, 1));
%     end
%     fprintf('Average homography reprojection MSE: %f\n',mean(MSE));
    
    %% 3.4 Intrinsic parameter estimation
    % focal length and optical center
    X_system = zeros(2*I,3+I);
    for i = 1:I
        M = [H{i}(1,1),H{i}(1,2);
            H{i}(2,4)/H{i}(3,1),H{i}(2,5)/H{i}(3,2);
             H{i}(3,1),H{i}(3,2)];
         
        X_system(2*i-1,1:3) = [M(1,1)*M(1,2), M(1,1)*M(3,2)+M(1,2)*M(3,1), M(3,1)*M(3,2)];
        X_system(2*i-1,3+i) = M(2,1)*M(2,2);
        X_system(2*i  ,1:3) = [M(1,1)^2 - M(1,2)^2, 2*(M(1,1)*M(3,1)-M(1,2)*M(3,2)), M(3,1)^2-M(3,2)^2];
        X_system(2*i  ,3+i) = M(2,1)^2 - M(2,2)^2;
    end
    [U,S,V] = svd(X_system,0);
    u0 = -V(2,3+I)/V(1,3+I);
    f = sqrt(V(3,3+I)/V(1,3+I)-u0^2);
    
    % scanning speed and camera distances
    X_system = zeros(3*I,1+I);
    X_rhs = repmat([1;1;0],[I,1]);
    for i = 1:I
        M = [H{i}(1,1),H{i}(1,2);
             H{i}(2,1)-H{i}(3,1)*H{i}(2,3),H{i}(2,2)-H{i}(3,2)*H{i}(2,3);
             H{i}(3,1),H{i}(3,2)];
        
        X_system(3*i-2,1  ) = M(2,1)^2;
        X_system(3*i-2,1+i) = ((u0^2+f^2)*M(3,1)^2 - 2*u0*M(1,1)*M(3,1) + M(1,1)^2)/f^2;
        X_system(3*i-1,1  ) = M(2,2)^2;
        X_system(3*i-1,1+i) = ((u0^2+f^2)*M(3,2)^2 - 2*u0*M(1,2)*M(3,2) + M(1,2)^2)/f^2;
        X_system(3*i  ,1  ) = M(2,1)*M(2,2);
        X_system(3*i  ,1+i) = ((u0^2+f^2)*M(3,1)*M(3,2) - u0*(M(1,2)*M(3,1)+M(1,1)*M(3,2)) + M(1,1)*M(1,2))/f^2;
    end
    res = X_system\X_rhs;
    s = sqrt(1/res(1));
    t3 = sqrt(res(2:end));
    
    %this gives us the entire intrinsic matrix:
    K = [[f,0,u0];[0,s,0];[0,0,1]];
    
    %% 3.5 Extrinsic parameter estimation
    R = cell(1,I);
    t = cell(1,I);
    for i = 1:I
        t{i} = [(H{i}(1,3)-u0)*t3(i)/f;H{i}(2,3)/s;t3(i)];
        R{i} = [(H{i}(1,1:2)-u0*H{i}(3,1:2))*t{i}(3)/f;H{i}(2,1:2)/s-H{i}(3,1:2)*t{i}(2);H{i}(3,1:2)*t{i}(3)];
        R{i} = [R{i},cross(R{i}(:,1),R{i}(:,2))];
        [U,~,V] = svd(R{i});
        R{i} = U*V';        
    end
end