%refine_linear_camera
%   Refines the intrinsic and extrinsic matrices of a linear camera
%   calibration using bundle adjustment.
%   Usage:
%       K = refine_linear_camera(measurements,points,K,R,t)
%
%   Original code by Simon Donné, January 2017
function [K,R,t,final_MSE] = refine_linear_camera(uv,ab,K,R,t,do_equivalencies)
    if nargin < 6
        do_equivalencies = 1;
    end

    I = numel(R);
    N = size(ab,2);
    a = ab(1,:)';
    b = ab(2,:)';
    
    lambda = 0.01;%damping factor for levenberg marquardt
                 %this is a pretty easy problem, so it is not crucial
    iterations = 100;
    MSE = zeros(iterations,1);
    
    for iteration = 1:iterations
        f = K(1,1);
        u0 = K(1,3);
        s = K(2,2);
        
        %% calculate the current reprojection errors and the jacobian
        Jacobian = zeros(2*I*N,3+6*I);
        errorvector = zeros(2*I*N,1);
        SE = zeros(I*N,1);
        
        ddf = 0;
        ddu0 = 0;
        dds = 0;
        
        for i = 1:I
            reprojected = K*(R{i}*ab+repmat(t{i},[1,N]));
            e1 = uv{i}(1,:) - reprojected(1,:)./reprojected(3,:);
            e2 = uv{i}(2,:) - reprojected(2,:);
            errorvector(2*N*(i-1) + (1:2*N)) = [e1,e2];
            SE(N*(i-1) + (1:N)) = e1.^2+e2.^2;
            
            L1 = a*R{i}(1,1) + b*R{i}(1,2) + t{i}(1);
            L2 = a*R{i}(3,1) + b*R{i}(3,2) + t{i}(3);
            L3 = a*R{i}(2,1) + b*R{i}(2,2) + t{i}(2);
            L4 = a*R{i}(1,1) + b*R{i}(1,2);
            L5 = a*R{i}(3,1) + b*R{i}(3,2);
            L6 = a*R{i}(2,1) + b*R{i}(2,2);
            
            de1df  = -L1./L2;
            de1du0 = -ones(N,1);
            de2ds  = -L3;
            
            ddf  = ddf + sum(e1.*de1df');
            ddu0 = ddu0+ sum(e1.*de1du0');
            dds  = dds + sum(e2.*de2ds');
            
            Jacobian( 2*N*(i-1) + (1:N), 1) = de1df;
            Jacobian( 2*N*(i-1) + (1:N), 2) = de1du0;
            Jacobian( 2*N*(i-1) + (1:N)+N, 3) = de2ds;
            
            de1dt1 = -f./L2;
            de1dt3 = f*L1./(L2.^2);
            de2dt2 = -s*ones(N,1);
            
            de1dw1 = L6.*L1*f./(L2.^2);
            de1dw2 = -f*(L1.*L4*f+L2.*L5)./(L2.^2);
            de1dw3 = L6./L2*f;
            de2dw1 = s*L5;
            de2dw3 = -s*L4;
            
            Jacobian( 2*N*(i-1) + (1:N), 3 + 6*(i-1)+1) = de1dt1;
            Jacobian( 2*N*(i-1) + (1:N), 3 + 6*(i-1)+3) = de1dt3;
            Jacobian( 2*N*(i-1) + (1:N)+N, 3 + 6*(i-1)+2) = de2dt2;
            
            Jacobian( 2*N*(i-1) + (1:N), 3 + 6*(i-1)+4) = de1dw1;
            Jacobian( 2*N*(i-1) + (1:N), 3 + 6*(i-1)+5) = de1dw2;
            Jacobian( 2*N*(i-1) + (1:N), 3 + 6*(i-1)+6) = de1dw3;
            Jacobian( 2*N*(i-1) + (1:N) + N, 3 + 6*(i-1)+4) = de2dw1;
            Jacobian( 2*N*(i-1) + (1:N) + N, 3 + 6*(i-1)+6) = de2dw3;
        end
        MSE(iteration) = mean(SE);
        
        %% Levenberg-Marquardt: calculate optimization step
        JTJ = Jacobian'*Jacobian;
        delta = -(JTJ+lambda*diag(diag(JTJ)))\(Jacobian'*errorvector);
        
        %% Levenberg-Marquardt: apply the optimization step
        %intrinsic parameters, update the intrinsic matrix

        if do_equivalencies
            f = f + delta(1);
            u0= u0+ delta(2);
        end
        s = s + delta(3);
        K = [[f,0,u0];[0,s,0];[0,0,1]];        
        %extrinsic parameters
        for i = 1:I
            t{i}(1:3) = t{i}(1:3) + delta(3 + 6*(i-1)+(1:3));
            phi   = delta(3 + 6*(i-1)+4);
            theta = delta(3 + 6*(i-1)+5);
            psi   = delta(3 + 6*(i-1)+6);
            Deltai = [[1,0,0];[0,cos(phi),-sin(phi)];[0,sin(phi),cos(phi)]] * ...
                     [[cos(theta),0,sin(theta)];[0,1,0];[-sin(theta),0,cos(theta)]] * ...
                     [[cos(psi),-sin(psi),0];[sin(psi),cos(psi),0];[0,0,1]];
            R{i} = Deltai*R{i};
        end
    end
    
    final_MSE = MSE(end);
%     plot(MSE);
%     keyboard
    fprintf('After bundle adjustment, the MSE is %4.4f\n',final_MSE);
end