function [X, Z] = ekf1(sensor, vic, varargin)
% EKF1 Extended Kalman Filter with Vicon velocity as inputs
%
% INPUTS:
%   sensor - struct stored in provided dataset, fields include
%          - is_ready: logical, indicates whether sensor data is valid
%          - t: sensor timestamp
%          - rpy, omg, acc: imu readings
%          - img: uint8, 240x376 grayscale image
%          - id: 1xn ids of detected tags
%          - p0, p1, p2, p3, p4: 2xn pixel position of center and
%                                four corners of detected tags
%            Y
%            ^ P3 == P2
%            | || P0 ||
%            | P4 == P1
%            o---------> X
%   vic    - struct for storing vicon linear velocity in world frame and
%            angular velocity in body frame, fields include
%          - t: vicon timestamp
%          - vel = [vx; vy; vz; wx; wy; wz]
%   varargin - any variables you wish to pass into the function, could be
%              a data structure to represent the map or camera parameters,
%              your decision. But for the purpose of testing, since we don't
%              know what inputs you will use, you have to specify them in
%              init_script by doing
%              ekf1_handle = ...
%                  @(sensor, vic) ekf1(sensor, vic, your input arguments);
%
% OUTPUTS:
% X - nx1 state of the quadrotor, n should be greater or equal to 6
%     the state should be in the following order
%     [x; y; z; roll; pitch; yaw; other states you use]
%     we will only take the first 6 rows of X
% OPTIONAL OUTPUTS:
% Z - mx1 measurement of your pose estimator, m shoulb be greater or equal to 6
%     the measurement should be in the following order
%     [x; y; z; roll; pitch; yaw; other measurement you use]
%     note that this output is optional, it's here in case you want to log your
%     measurement

persistent x_hat_old P old_t first_time phi theta psi
persistent xdot_func jacobian_x_simp_func jacobian_n_simp_func

if isempty(first_time)
    first_time=0;
    init_script2
    x_hat_old=zeros(6,1);
    old_t=0;
    P=eye(6);
    X=zeros(6,1);
    Z=zeros(6,1);
    phi=0;
    theta=0;
    psi=0;
    return
end

X = x_hat_old;
Z=zeros(6,1);

if (isempty(sensor) && isempty(vic))
    X=nan(6,1);
    Z=nan(6,1);
    return
end


%%%%%%%%%%%%%%%%%%%%%% PREDICTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make At, Ut, Ft Vt

if ~isempty(vic)
    dt=vic.t-old_t;
    old_t=vic.t;
    v_m1=vic.vel(1);
    v_m2=vic.vel(2);
    v_m3=vic.vel(3);
    omega1=vic.vel(4)
    omega2=vic.vel(5)
    omega3=vic.vel(6)
    phi
    theta
    v_m1
    v_m2
    v_m3
    
    
    xdot=xdot_func(0,0,0,0,0,0,omega1,omega2,omega3,phi,theta,v_m1,v_m2,v_m3);
    At=jacobian_x_simp_func(0,0,omega1,omega3,phi,theta);
    Ut=jacobian_n_simp_func(phi,theta);
    dt
    At
    Ft=eye(6)+dt*At;
    Vt=dt*Ut;
    
    Q=eye(6); %%%% MAYBE???
    
    % make prediction
    X=x_hat_old+xdot*dt;
    % predicted estimate covariance
    P=Ft*P*Ft'+Vt*Q*Vt';
end

%%%%%%%%%%%%%%%%%%%%%% UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(sensor)
    
    if ((sensor.is_ready) && ~isempty(sensor.id))
        
        
        [pos, eul,t] = estimate_pose(sensor, varargin);
        
        phi=eul(1);
        theta=eul(2);
        psi=eul(3);
        
        Z=[pos;eul];
        
        
        
        
        
        % make some variables
        Ct=eye(6);
        Wt=eye(6);
        R=eye(6);
        
        
        % some other calculations
        % innovation or measurement residual
        y_tilda=[pos;eul]-eye(6)*X;
        % innovation (or residual) covariance
        Sk=Ct*P*Ct'+Wt*R*Wt';
        % Optimal kalman gain
        K=P*Ct'/Sk;
        
        % update estimations
        X=X+K*y_tilda;
        P=P-K*Ct*P;
        
    end
    
end

% update old variables
x_hat_old=X;

end
