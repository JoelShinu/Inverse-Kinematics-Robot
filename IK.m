function Q = IK(DH_params, jtype, q, pdes, odes)
% IK calculates the inverse kinematics of a manipulator
% Q = IK(DH_params, jtype, q, pdes, odes) calculates the joint
% values of a manipulator robot given the desired end-effector's position
% and orientation by solving the inverse kinematics, iteratively.  The
% inputs are DH_params, jtype and q, pdes, and odes.  DH_params is
% nx4 matrix including Denavit-Hartenberg parameters of the robot.  jtype
% and q are n-dimensional vectors.  jtype describes the joint types of the
% robot and its values are either 0 for revolute or 1 for prismatic joints.
% q is the joint values.  pdes and odes are desired position and
% orientation of the end-effector.  pdes is 2x1 for a planar robot and 3x1
% for a spatial one.  odes is a scalar for a planar robot and a 3x1 vector
% for a 3D robot.  Orientation is defined as roll-pitch-yaw angles.  Use
% Jacobian Transpose to compute the inverse of the Jacobian.

n = size(q,1);  % robot's DoF
% consistency check
if (n~=size(DH_params,1)) || (n~=size(jtype,1))
    error('inconsistent in dimensions');
end

% robot's dimension (planar or spatial)
dim = size(pdes,1);
if ~((dim==2) || (dim==3))
    error('size of pdes must be either 2x1 for a planar robot or 3x1 for a spatial robot');
end

% check the size of orientation parameter
if dim==2
    if size(odes,1)~=1
        error('desired orientation must be scalar for a planar robot');
    end
else
    if size(odes,1)~=3
        error('desired orientation must be 3x1 for a spatial robot');
    end
end

%% iterative inverse kinematics
%%% Complete here %%%

% Parameters %
max_iteration=100;
K=0.2;

x_dot = [0; 0; 0];
q_dot = [0; 0; 0];  

iteration=0;

% Jacobian Calculation %
[trans, jacob] = FK(DH_params,jtype,q)

x_des = [pdes;odes];
Q_log = q;

while iteration < max_iteration
        [trans, jacob] = FK(DH_params,jtype,q);
        jacob(3,:) = [];
        jacob(3,:) = [];
        jacob(3,:) = [];
        trans;

        J_plus = transpose(jacob)*inv(jacob*transpose(jacob));
        p = trans(1:2,4);
        r = trans(1:3,1:3);
        phi = atan2(r(2,1),r(1,1));

        x = [p;phi];

        x_delta=x_des-x;

        q_dot=J_plus*(K*x_delta + x_dot);

        q=q+q_dot;
        
        Q_log=cat(2,Q_log,q);

end

