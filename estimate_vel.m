function [vel, omg] = estimate_vel(sensor, varargin)
%ESTIMATE_VEL 6DOF velocity estimator
%   sensor - struct stored in provided dataset, fields include
%          - is_ready: logical, indicates whether sensor data is valid
%          - t: timestamp
%          - rpy, omg, acc: imu readings, you should not use these in this phase
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
%              estimate_vel_handle = ...
%                  @(sensor) estimate_vel(sensor, your personal input arguments);
%   vel - 3x1 velocity of the quadrotor in world frame
%   omg - 3x1 angular velocity of the quadrotor

if isempty(sensor.id)
    vel = [];
    omg = [];
    return
end

persistent old_points old_timestamp pointTracker

vel = zeros(3,1);
omg = zeros(3,1);


K=[314.1779 0         199.4848; ...
    0         314.2218  113.7838; ...
    0         0         1];

% rename variable so it works with the variable names from the script. 
    data=sensor;

% if first time running
    if (isempty(pointTracker))
        image1=data.img;
        image1Points = detectSURFFeatures(image1);
        [image1Features, image1Points] = extractFeatures(image1, image1Points);

        pointTracker = vision.PointTracker;
        initialize(pointTracker,image1Points.Location,image1)
        disp('initializing')
        old_points=image1Points.Location;
        old_timestamp=data.t;
    end

% find the new coordinates of these points
    [points,point_validity,scores] = step(pointTracker,data.img);


% extract the valid points between images to use for estimating
% velocity
    valid_old_points=old_points(point_validity,:);
    valid_new_points=points(point_validity,:);
    old_points=points;

% calculate the Z used to help estimate velocity
    [pos, eul,t] = estimate_pose(data);
    R = RPYtoRot_ZXY(eul(1),eul(2),eul(3));
    R_c=(R\[sind(45), -sind(45), 0; -sind(45),-sind(45),0;0,0,-1])'; %here 

% establish the change in time
    delta_T=.02; %20 milliseconds, .02 seconds. does better than timestamp history.

% which indices to use? (go back and do ransac)
    point_indices=1:length(valid_old_points);

    [~,inlierpoints1,inlierpoints2] = estimateGeometricTransform(valid_new_points,valid_old_points,'projective');
    point_indices=1:length(inlierpoints1);
    
% calculate VOmega
    %[VOmega] = estimate_vel_3_pts(point_indices, valid_new_points, valid_old_points, delta_T, R_c, pos,t);
    [VOmega] = estimate_vel_3_pts(point_indices, inlierpoints1, inlierpoints2, delta_T, R_c, pos,t);

% rotate VOmega into global frame;
    %%%%%%% NEED TO COME BACK AND ROTATE THIS
    %vel=[sind(45), -sind(45), 0; -sind(45),-sind(45),0;0,0,-1]*VOmega(1:3);
    %omg=[sind(45), -sind(45), 0; -sind(45),-sind(45),0;0,0,-1]*VOmega(4:6);
    
    vel=R_c'*VOmega(1:3);
    omg=R_c*VOmega(4:6);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decide if need to recalculate new points for next round
if (sum(point_validity)/length(point_validity)) < .7
        image1=data.img;
        image1Points = detectSURFFeatures(image1);
        [image1Features, image1Points] = extractFeatures(image1, image1Points);
        setPoints(pointTracker,image1Points.Location);
        
        old_points=image1Points.Location;
        old_timestamp=data.t;

    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
