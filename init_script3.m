% Add additional inputs after the given ones if you want to
% Example:
% your_input = 1;
% ekf_handle1 = @(sensor, vic) ekf1(sensor, vic, your_input);
% ekf_handle2 = @(sensor) ekf2(sensor, your_input);
%
% We will only call ekf_handle in the test function.
% Note that this will only create a function handle, but not run the function



syms theta1 phi1 psi1 b_g omega v_m1 v_m2 v_m3 n_v x_3 n_g p1 p2 p3 
syms nb1 nb2 nb3 nb4 nb5 nb6 nb7

G=[cos(theta1) 0 -cos(phi1)*sin(theta1);
     0 1 sin(phi1);
     sin(theta1) 0 cos(phi1)*cos(theta1)];

 
R = [cos(psi1)*cos(theta1) - sin(phi1)*sin(psi1)*sin(theta1), ...
     cos(theta1)*sin(psi1) + cos(psi1)*sin(phi1)*sin(theta1), ...
     -cos(phi1)*sin(theta1); ...
     -cos(phi1)*sin(psi1),...
     cos(phi1)*cos(psi1), ...
     sin(phi1);...
     cos(psi1)*sin(theta1) + cos(theta1)*sin(phi1)*sin(psi1),...
     sin(psi1)*sin(theta1) - cos(psi1)*cos(theta1)*sin(phi1),...
     cos(phi1)*cos(theta1)]; 
 
omega=sym('omega',[3,1]);
x3=sym('x3',[3,1]);
x4=sym('x4',[3,1]);
n_g=sym('n_g',[3,1]);

g=sym('g',[3,1]);
a_m=sym('a_m',[3,1]);
x5=sym('x5',[3,1]);
n_a=sym('n_a',[3,1]);

n_bg=sym('n_bg',[3,1]);

n_ba=sym('n_ba',[3,1]);


x_dot1=[x3; inv(G)*(omega-n_g); g+R*(a_m-n_a)]

%% jacobian with rescpect to X

jacobian_x=jacobian(x_dot1, [p1;p2;p3;v_m1;v_m2;v_m3;phi1;theta1;psi1]);
jacobian_x_simp=simplify(jacobian_x);
jacobian_x_simp_func=matlabFunction(jacobian_x_simp);

%% jacobian with respect to n

jacobian_n=jacobian(x_dot1, [n_g; n_a; nb1;nb2;nb3]);
jacobian_n_simp=simplify(jacobian_n);
jacobian_n_simp_func=matlabFunction(jacobian_n_simp);

%% xdot funct
xdot_func=matlabFunction(x_dot1);

%% h jacobian
% g=[p1;p2;p3;phi1;theta1;psi1];
% h_jacobian_x=simplify(jacobian(g,[p1;p2;p3;phi1;theta1;psi1]))
% h_jacobian_x_func=matlabFunction(h_jacobian_x)





