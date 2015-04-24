function [VOmega] = estimate_vel_3_pts(point_indices, valid_new_points, valid_old_points, delta_T, R_c, pos,t) 

K=[314.1779 0         199.4848; ...
0         314.2218  113.7838; ...
0         0         1];

[m,n]=size(point_indices);

% locations
p=valid_new_points;%(point_indices,:)
p_old=valid_old_points;%(point_indices,:)
%p1=(valid_new_points(p1_index,:));%+valid_old_points(p1_index))/2;
%p2=(valid_new_points(p2_index,:));%+valid_old_points(p2_index))/2;
%p3=(valid_new_points(p3_index,:));%+valid_old_points(p3_index))/2;

vp=[K\[p,ones(length(p),1)]']';
vp_old=[K\[p_old,ones(length(p),1)]']';


% velocities
p_dot=(vp(:,1:2)-vp_old(:,1:2))/delta_T;

% p1_dot=(valid_new_points(p1_index,:)-valid_old_points(p1_index,:))./delta_T;
% p2_dot=(valid_new_points(p2_index,:)-valid_old_points(p2_index,:))./delta_T;
% p3_dot=(valid_new_points(p3_index,:)-valid_old_points(p3_index,:))./delta_T;
% p_dot=[p1_dot';p2_dot';p3_dot'];

% calculate the corresponding depths
Z=(R_c(:,3)'*t)./(R_c(:,3)'*inv(K)*[p(:,1),p(:,2),ones(n,1)]'); %HERE
%Z2=(R_c(:,3)'*t)/(R_c(:,3)'*inv(K)*[p2(1);p2(2);1]);
%Z3=(R_c(:,3)'*t)/(R_c(:,3)'*inv(K)*[p3(1);p3(2);1]);

% x1=p1(1);
% y1=p1(2);
% x2=p2(1);
% y2=p2(2);
% x3=p3(1);
% y3=p3(2);

f_matrix=zeros(2*n,6);
p_dot2=zeros(2*n,1);
for w=1:n;
    f_matrix(2*w-1,:)=[-1/Z(w), 0, vp(w,1)/Z(w), vp(w,1)*vp(w,2),-(1+vp(w,1)^2), vp(w,2)];
    f_matrix(2*w,:)=[0, -1/Z(w), vp(w,2)/Z(w), (1+vp(w,2)^2), -vp(w,1)*vp(w,2), -vp(w,1)];
    p_dot2(2*w-1)=p_dot(w,1);
    p_dot2(2*w)=p_dot(w,2);

end

% % f_matrix used for solving for VOmega
% f_matrix=[-1/Z1,      0,  x1/Z1,     x1*y1,  -(1+x1^2),   y1;
%               0,  -1/Z1,  y1/Z1,  (1+y1^2),     -x1*y1,  -x1;
%           -1/Z2,      0,  x2/Z2,     x2*y2,  -(1+x2^2),   y2;
%               0,  -1/Z2,  y2/Z2,  (1+y2^2),     -x2*y2,  -x2;
%           -1/Z3,      0,  x3/Z3,     x3*y3,  -(1+x3^2),   y3;
%               0,  -1/Z3,  y3/Z3,  (1+y3^2),     -x3*y3,  -x3];

    
VOmega = f_matrix\p_dot2;