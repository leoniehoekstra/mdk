% Computes the configuration all the links of a serial chain mechanism relative to the base frame Psi_0
% Input: unit_twists - cell array of n 6x1 matrices - unit twists of the joints T_i^{(i-1), (i-1)}
%        H0s - cell array of n 4x4 matrices - reference (q = 0) configurations of the links H_i^{i-1}(0)
%        q - nx1 matrix - joint variables q^i
% Output: Hs - cell array of n+1 4x4 matrices - current configurations of the links 0 to n of the serial chain mechanism relative to the base frame Psi_0

function Hs = direct_kinematics(unit_twists, H0s, q)
    n = length(unit_twists);
    Hs = cell(1, n);
    T_base = cell(1, n+1);
    eT = cell(1, n+1);

    % >>>>>> INSERT YOUR CODE HERE <<<<<<

    Hs{1} = eye(4); % Initial base frame
    H0c = cumulativeH(H0s);

    for i = 1:n
        % Get unit twist in base frame

        if i == 1
            T_base{i} = express_in_base(unit_twists{i}, eye(4));
        else
            T_base{i} = express_in_base(unit_twists{i}, H0c{i-1});
        end
       

        % Chain multiplication
        if i == 1
            eT{i} = twist_exp(T_base{i}, q(i));
        else
            eT{i} = twist_exp(T_base{i}, q(i)) * eT{i-1};
        end

        Hs{i} = eT{i} * H0c{i};
    end

end

function T_base = express_in_base(unit_twist, H0c)
    % This function takes a unit twist and H0 and returns the
    % unit twist expressed in the base frame.
    
    Ad_H = get_adjoint(H0c);
    T_base = Ad_H * unit_twist;
end


function Ad_H = get_adjoint(H0c)
    % This function takes a H0 and returns its adjoint matrix.

    R = H0c(1:3, 1:3);       % Extract rotation matrix
    o = H0c(1:3, 4);         % Extract position vector

    o_skew = skew(o);       % Get skew matrix from o

    Ad_H = [R, zeros(3,3); o_skew*R, R];    % return Ad_H
end


function S = skew(o)
    % This function takes a vector o and returns its skew matrix
    n = length(o);

    if n == 3
        S = [0, -o(3),  o(2);
            o(3),   0, -o(1);
            -o(2), o(1),   0];

    elseif n == 6
        omega = o(1:3);
        v = o(4:6);
        S = [skew(omega), v; 0 0 0 0];

    else
        error('length of o should be 3 or 6');

    end 
end


function H_c = cumulativeH(H0s)
    % Comuputes the cumulative H matrix from the H0's

    n = numel(H0s);
    H_c = cell(1, n);
    H_c{1} = H0s{1};

    for i = 2:n
        H_c{i} = H_c{i-1} * H0s{i};
    end
end


function eT = twist_exp(T, q)
    omega = T(1:3);             % angular part
    omega_skew = skew(omega);   % skew of angular part
    v = T(4:6);                 % linear part
    
    % Rodrigues' formula for rotation matrix e^omega_skew:
    e_omega_skew = eye(3) + sin(q)*omega_skew + (1-cos(q))*(omega_skew*omega_skew);

    % Calculate translation vector e^v:
    t = (1/(norm(omega) * norm(omega)) * ((eye(3)-e_omega_skew)*cross(omega,v)+transpose(omega)*v*omega));

    % Return matrix eT:
    eT = [e_omega_skew, t; 0 0 0 1];
end

