function [X, Z] = ekf2(sensor, varargin)
% EKF2 Extended Kalman Filter with IMU as inputs
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
%   varargin - any variables you wish to pass into the function, could be
%              a data structure to represent the map or camera parameters,
%              your decision. But for the purpose of testing, since we don't
%              know what inputs you will use, you have to specify them in
%              init_script by doing
%              ekf1_handle = ...
%                  @(sensor) ekf2(sensor, your input arguments);
%
% OUTPUTS:
% X - nx1 state of the quadrotor, n should be greater or equal to 9
%     the state should be in the following order
%     [x; y; z; vx; vy; vz; roll; pitch; yaw; other states you use]
%     we will only take the first 9 rows of X
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
    init_script3
    x_hat_old=zeros(9,1);
    old_t=0;
    P=eye(9);
    X=zeros(9,1);
    Z=zeros(9,1);
    phi=0;
    theta=0;
    psi=0;
    return
end

X = x_hat_old;
Z=zeros(9,1);

if (isempty(sensor) && isempty(vic))
    X=nan(9,1);
    Z=nan(9,1);
    return
end


%%%%%%%%%%%%%%%%%%%%%% PREDICTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make At, Ut, Ft Vt

if ~isempty(sensor)
    
    if ((sensor.is_ready) && ~isempty(sensor.id))
        if (sensor.id(1)~=0)
            
            dt=sensor.t-old_t;
            old_t=sensor.t;
            
            [pos, eul,t] = estimate_pose(sensor, varargin);
            [vel, ~] = estimate_vel(sensor, varargin);
            phi=eul(1);
            theta=eul(2);
            psi=eul(3);
            Z=[pos;eul;vel];
            
            
            
            v_m1=vel(1);
            v_m2=vel(2);
            v_m3=vel(3);
            omega1=sensor.omg(1);
            omega2=sensor.omg(2);
            omega3=sensor.omg(3);
            a_m1=sensor.acc(1);
            a_m2=sensor.acc(2);
            a_m3=sensor.acc(3);
            
            %         xdot=xdot_func(a_m1,a_m2,a_m3,0,0,-9.81,0,0,0,0,0,0,0,0,0,0,0,0,omega1,omega2,omega3,phi,psi,theta,v_m1,v_m2,v_m3,0,0,0,0,0,0)
            %         At=jacobian_x_simp_func(a_m1,a_m2,a_m3,0,0,0,0,0,omega1,omega3,phi,psi,theta,0,0,0,0,0)
            %         Ut=jacobian_n_simp_func(phi,psi,theta)
            
            xdot=xdot_func(a_m1,a_m2,a_m3,0,0,-9.81,0,0,0,0,0,0,omega1,omega2,omega3,phi,psi,theta,v_m1,v_m2,v_m3);
            At=jacobian_x_simp_func(a_m1,a_m2,a_m3,0,0,0,0,0,omega1,omega3,phi,psi,theta)
            Ut=jacobian_n_simp_func(phi,psi,theta)
            
            %xdot=xdot_func(0,0,0,0,0,0,omega1,omega2,omega3,phi,theta,v_m1,v_m2,v_m3);
            %At=jacobian_x_simp_func(0,0,omega1,omega3,phi,theta);
            %Ut=jacobian_n_simp_func(phi,theta);
            
            Ft=eye(9)+dt*At
            Vt=dt*Ut
            
            Q=eye(9); %%%% MAYBE???
            
            % make prediction
            X=x_hat_old+xdot*dt
            % predicted estimate covariance
            P=Ft*P*Ft'+Vt*Q*Vt';
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%% UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(sensor)
    
    if ((sensor.is_ready) && ~isempty(sensor.id))
        if (sensor.id(1)~=0)
            
            % make some variables
            Ct=eye(9);
            Wt=eye(9);
            R=eye(9);
            
            
            % some other calculations
            % innovation or measurement residual
            y_tilda=[pos;vel;eul]-eye(9)*X;
            % innovation (or residual) covariance
            Sk=Ct*P*Ct'+Wt*R*Wt';
            % Optimal kalman gain
            K=P*Ct'/Sk;
            
            % update estimations
            X=X+K*y_tilda;
            P=P-K*Ct*P;
            
        end
    end
end

% update old variables
x_hat_old=X;


end
