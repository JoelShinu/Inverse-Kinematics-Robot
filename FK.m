function [T, J] = FK(DH_params, jtype, q)
% FK calculates the forward kinematics of a manipulator
% [T, J] = FK(DH_params, jtype, q) calculates the Homogeneous
% transformation matrix from base to the end-effector (T) and also the
% end-effector Jacobian with respect to the base frame of a manipulator
% robot.  The inputs are DH_params, jtype and q.  DH_params is nx4 matrix
% including Denavit-Hartenberg parameters of the robot.  jtype and q are
% n-dimensional vectors.  jtype describes the joint types of the robot.
% Its values are either 0 for revolute or 1 for prismatic joints.  q is the
% vector of joint values.

n = size(q,1);  % robot's DoF
% consistency check
if (n~=size(DH_params,1)) || (n~=size(jtype,1))
    error('inconsistent in dimensions');
end

% initialisation
T = eye(4,4);
J = zeros(6,n);
%% T homogeneous matrix calculations
out_Homogenous = eye(4,4);
for p = 1:n
    %setting the variables
    curr_ai = DH_params(p,1);
    curr_vi = q(p);
    curr_alphai = DH_params(p,2);
    curr_di = DH_params(p,3);

    %intialisation of current homogenous equation
    curr_hom = zeros(4,4);
    
    %column 1
    curr_hom(1,1)  = cos(curr_vi);
    curr_hom(2,1) = sin(curr_vi);
    %column 2
    curr_hom(1,2) = -sin(curr_vi)*cos(curr_alphai);
    curr_hom(2,2) = cos(curr_vi)*cos(curr_alphai);
    curr_hom(3,2) = sin(curr_alphai);
    %column 3
    curr_hom(1,3) = sin(curr_vi)*sin(curr_alphai);
    curr_hom(2,3) = -cos(curr_vi)*sin(curr_alphai);
    curr_hom(3,3) = cos(curr_alphai);
    %column 4
    curr_hom(1,4) = curr_ai*cos(curr_vi);
    curr_hom(2,4) = curr_ai*sin(curr_vi);
    curr_hom(3,4) = curr_di;
    curr_hom(4,4) = 1;

    %multiply the previous homogenous with current
    out_Homogenous = out_Homogenous*curr_hom;
    %each step is saved to be used ofr the jacobian
    homogenous{p+1} = out_Homogenous;
   
end
T=out_Homogenous;

%% J Jacobian calculations
%inital values
R{1} = eye(3,3);
d{1} = [0;0;0];

%creates rotational matrices
for k = 1 : n
    ch = homogenous{k+1};
    R{k+1} = ch(1:3,1:3);
    ch(1:3,4);
    d{k+1} = ch(1:3,4);
end
jacob = 0;
for j =1 : n
    if jtype(j)==1  %prismatic
        %%%% complete here %%%%
        linear = R{j}*[0;0;1];
        rotational = [0;0;0];
    else             % revolute  
        %%%% complete here %%%%
        linear = cross((R{j}*[0;0;1]),(d{n+1}-d{j}));
        rotational = R{j}*[0;0;1]; 
    end

    %concatenating the matrices to create the jacobian
    curr_colm = cat(1,linear,rotational);
    if j ==1
        jacob = curr_colm;
    else
        jacob = cat(2,jacob,curr_colm);

    end   

end
J = jacob;
end
