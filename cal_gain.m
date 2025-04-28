function [K1, K2, K3, K4] = cal_gain()
    % 1. Define system parameters
    r = 4;   % Number of fuzzy rules
    n_x = 2; % State dimension
    n_w = 2; % Disturbance dimension (B_w is a 2x2 matrix)
    n_u = 1; % Input dimension
    n_z = 2; % Output dimension (dimension of z)
    
    % System matrices (based on provided values)
    A1 = [-9 12; 0 -10];       % State matrix for rule 1 of y2
    A2 = [-7 8; 0 -6];         % State matrix for rule 2 of y2
    Ad1 = [-0.4 0; -0.4 -0.4]; % Delay matrix for y2(x,t-tau1) rule 1
    Ad2 = [0.4 0; -0.4 -0.4];  % Delay matrix for y2(x,t-tau1) rule 2
    Bi = [1; 1];               % Input matrix
    Bw = [-0.5 1; -1 0.2];     % Disturbance input matrix
    C = [0.2 -1; 1 -0.2];      % Output matrix
    
    % Given parameters
    tau1_bar = 0.1*pi;     % tau1_bar
    tau2_bar = 0.15*pi;    % tau2_bar
    gamma = 1;             % Disturbance attenuation level gamma
    epsilon = 1;           % epsilon (ε)
    
    % 2. Initialize LMI system
    setlmis([]); % Clear LMI system
    
    % 3. Define LMI variables
    struct = [ones(n_x, 1), zeros(n_x, 1)];
    P_bar = lmivar(1, struct);   % Diagonal positive definite matrix P_bar (type 1: symmetric matrix, n_x x n_x)
    Q1_bar = lmivar(1, [n_x 1]); % Symmetric positive definite matrix Q1_bar
    Q2_bar = lmivar(1, [n_x 1]); % Symmetric positive definite matrix Q2_bar
    X1 = lmivar(2, [n_u n_x]);   % Matrix X1 (type 2: rectangular matrix, n_u x n_x)
    X2 = lmivar(2, [n_u n_x]);   % Matrix X2
    X3 = lmivar(2, [n_u n_x]);   % Matrix X3
    X4 = lmivar(2, [n_u n_x]);   % Matrix X4
    varepsilon = lmivar(1, [1 1]); % Scalar varepsilon (type 1: scalar)
    
    % 4. Define LMI constraints (for all i, j combinations)
    for i = 1:r
        for j = 1:r
            % Select corresponding matrices
            if i == 1
                Ai = A1;
                Adi = Ad1;
            elseif i == 2
                Ai = A1;
                Adi = Ad2;
            elseif i == 3
                Ai = A2;
                Adi = Ad1;
            else
                Ai = A2;
                Adi = Ad2;
            end
    
            if j == 1
                Xj = X1;
            elseif j == 2
                Xj = X2;
            elseif j == 3
                Xj = X3;
            else
                Xj = X4;
            end
    
            % LMI matrix (9x9 matrix)
            % The k-th LMI constraint
            k = (i-1)*r + j; % Unique identifier for the (i,j)-th constraint
            % 4.1 Sigma_11 (decomposed definition)
            lmiterm([k 1 1 P_bar], Ai, 1, 's'); % Ai*P_bar + P_bar*Ai'
            lmiterm([k 1 1 Q1_bar], 1, 1); % Q1_bar
            lmiterm([k 1 1 Q2_bar], 1, 1); % Q2_bar
            % 4.2 Adi*P_bar
            lmiterm([k 1 2 P_bar], Adi, 1); % Adi*P_bar
            % 4.3 Bi*Xj
            lmiterm([k 1 3 Xj], Bi, 1); % Bi*Xj
            % 4.4 B_w
            lmiterm([k 1 4 0], Bw); % B_w
            % 4.5 varepsilon*epsilon*I
            lmiterm([k 1 5 varepsilon], epsilon, eye(n_x)); % varepsilon*epsilon*I
            % 4.6 epsilon*P_bar
            lmiterm([k 1 6 P_bar], 1, 1); % epsilon*P_bar
            % 4.7 P_bar*C'
            lmiterm([k 1 9 P_bar], 1, C'); % P_bar*C'
            % 4.8 Sigma_22 (decomposed definition)
            lmiterm([k 2 2 Q1_bar], -(1 - tau1_bar), 1); % -(1 - tau1_bar)*Q1_bar
            % 4.9 epsilon*P_bar
            lmiterm([k 2 7 P_bar], 1, 1); % epsilon*P_bar
            % 4.10 Sigma_33 (decomposed definition)
            lmiterm([k 3 3 Q2_bar], -(1 - tau2_bar), 1); % -(1 - tau2_bar)*Q2_bar
            % 4.11 epsilon*Xj'
            lmiterm([k 3 8 -Xj], 1, 1); % epsilon*Xj' (using -Xj to represent transpose)
            % 4.12 -gamma^2*I
            lmiterm([k 4 4 0], -gamma^2*eye(n_w)); % -gamma^2*I
            % 4.13 -varepsilon*I
            lmiterm([k 5 5 varepsilon], -1, eye(n_x)); % -varepsilon*I
            % 4.14 -varepsilon*I
            lmiterm([k 6 6 varepsilon], -1, eye(n_x)); % -varepsilon*I
            % 4.15 -varepsilon*I
            lmiterm([k 7 7 varepsilon], -1, eye(n_u)); % -varepsilon*I
            % 4.16 -varepsilon*I
            lmiterm([k 8 8 varepsilon], -1, eye(n_u)); % -varepsilon*I (consistent with the number of columns of Xj')
            % 4.17 -I
            lmiterm([k 9 9 0], -eye(n_z)); % -I (dimension is n_z)
        end
    end
    
    % 5. Additional constraints: P_bar > 0, Q1_bar > 0, Q2_bar > 0, varepsilon > 0
    % P_bar > 0
    lmiterm([r*r+1 1 1 P_bar], -1, 1); 
    % Q1_bar > 0
    lmiterm([r*r+2 1 1 Q1_bar], -1, 1); 
    % Q2_bar > 0
    lmiterm([r*r+3 1 1 Q2_bar], -1, 1); 
    % varepsilon > 0
    lmiterm([r*r+4 1 1 varepsilon], -1, 1); 
    
    % 6. Get LMI system
    lmis = getlmis;
    
    % 7. Solve LMI (using feasp to find a feasible solution)
    [tmin, xfeas] = feasp(lmis);
    
    % Define the objective function 
    c = mat2dec(lmis, zeros(n_x), zeros(n_x), zeros(n_x), [1 2], [-2 -3], [-1 2], [3 4], 0);
    
    % Set optimization parameters
    options = [1e-5, 1000, 0, 10, 0]; 
    
    % Solve the optimization problem
    [copt, xopt] = mincx(lmis, c, options);
    
    % 8. Check the solution result
    if tmin < 0
        disp('LMI solution successful!');
        
        % Extract variable values
        P_bar_value = dec2mat(lmis, xopt, P_bar);
        X1_value = dec2mat(lmis, xopt, X1);
        X2_value = dec2mat(lmis, xopt, X2);
        X3_value = dec2mat(lmis, xopt, X3);
        X4_value = dec2mat(lmis, xopt, X4);
        varepsilon_value = dec2mat(lmis, xopt, varepsilon);
        
        % Display the value of varepsilon
        disp('Solved varepsilon:');
        disp(varepsilon_value);
        
        % Calculate feedback gains K_i = X_i * P_bar^(-1)
        K1 = X1_value / P_bar_value;
        K2 = X2_value / P_bar_value;
        K3 = X3_value / P_bar_value;
        K4 = X4_value / P_bar_value;   

        disp('Feedback gain K1:');
        disp(K1);
        disp('Feedback gain K2:');
        disp(K2);
        disp('Feedback gain K3:');
        disp(K3);
        disp('Feedback gain K4:');
        disp(K4);
    else
        disp('LMI solution failed!');
        if tmin == 0
            disp('LMI problem may be on the boundary, try adjusting numerical parameters.');
        else
            disp('LMI has no solution, please check the system matrices or LMI conditions!');
        end
    end
end