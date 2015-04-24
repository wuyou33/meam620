function [pos, eul,t] = estimate_pose(sensor, varargin)
%ESTIMATE_POSE 6DOF pose estimator based on apriltags
%   sensor - struct stored in provided dataset, fields include
%          - is_ready: logical, indicates whether sensor data is valid
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
%              estimate_pose_handle = ...
%                  @(sensor) estimate_pose(sensor, your personal input arguments);
%   pos - 3x1 position of the quadrotor in world frame
%   eul - 3x1 euler angles of the quadrotor

if isempty(sensor.id)
    pos = [];
    eul = [];
    t=[];
    return
end
%initialize some stuff


% parameters

% Camera Matrix (zero-indexed):
K=  [314.1779 0         199.4848; ...
    0         314.2218  113.7838; ...
    0         0         1];

% Camera-IMU Calibration (see attached images for details):
XYZ = [-0.04, 0.0, -0.03];
Yaw = pi/4;

% Tag ids:
ids=[0, 12, 24, 36, 48, 60, 72, 84,  96;
 1, 13, 25, 37, 49, 61, 73, 85,  97;
 2, 14, 26, 38, 50, 62, 74, 86,  98;
 3, 15, 27, 39, 51, 63, 75, 87,  99;
 4, 16, 28, 40, 52, 64, 76, 88, 100;
 5, 17, 29, 41, 53, 65, 77, 89, 101;
 6, 18, 30, 42, 54, 66, 78, 90, 102;
 7, 19, 31, 43, 55, 67, 79, 91, 103;
 8, 20, 32, 44, 56, 68, 80, 92, 104;
 9, 21, 33, 45, 57, 69, 81, 93, 105;
10, 22, 34, 46, 58, 70, 82, 94, 106;
11, 23, 35, 47, 59, 71, 83, 95, 107];

% correct the image coordinates
[I,J] = ind2sub([12,9],sensor.id+1); % returns indices of ids
target_p0_x=(0.152/2)+(2*0.152).*(I-1);
target_p0_y=(0.152/2)+(2*0.152).*(J-1)+(floor((J-1)./3))*(.178-.152);
%bottom left
target_p1_x=target_p0_x+(0.152/2);
target_p1_y=target_p0_y-(0.152/2);
%bottom right
target_p2_x=target_p0_x+(0.152/2);
target_p2_y=target_p0_y+(0.152/2);
%top right
target_p3_x=target_p0_x-(0.152/2);
target_p3_y=target_p0_y+(0.152/2);
%top left
target_p4_x=target_p0_x-(0.152/2);
target_p4_y=target_p0_y-(0.152/2);

% all the target coordinates
target_X=[target_p0_x,target_p1_x,target_p2_x,target_p3_x,target_p4_x]';
target_Y=[target_p0_y,target_p1_y,target_p2_y,target_p3_y,target_p4_y]';

% estimate homography
camera_x=[sensor.p0(1,:),sensor.p1(1,:),sensor.p2(1,:),sensor.p3(1,:),sensor.p4(1,:)]';
camera_y=[sensor.p0(2,:),sensor.p1(2,:),sensor.p2(2,:),sensor.p3(2,:),sensor.p4(2,:)]';

H = est_homography(camera_x,camera_y,target_X,target_Y);

this_matrix=K\H;
h1=this_matrix(:,1);
h2=this_matrix(:,2);
h3=this_matrix(:,3);

this_other_matrix=[h1,h2,cross(h1,h2)];

[U,S,V]=svd(this_other_matrix);

R=U*[1,0,0;0,1,0;0,0,det(U*V')]*V';
T=R'*(-h3/norm(h1));
t=h3/norm(h1);

pos=T+R'*XYZ';
R2=[sind(45), -sind(45), 0; -sind(45),-sind(45),0;0,0,-1]*R;

[phi,theta,psi] = RotToRPY_ZXY(R2);
eul=[phi;theta;psi];




% need to double check coordinates of target is spot on
% rotation angles, convert to euler angles 

%T=R'*(-h3/norm(h1));
%pos=T+R'*XYZ';


end
