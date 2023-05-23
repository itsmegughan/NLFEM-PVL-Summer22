%% Author: Gughan Dhanasekaran
%% Title: Program for PVL-Summer term 2022

clc
clear

%% Values from the table (constants) and user input
psi0 = 0.01;                                                                % mJ/mm^2
chi = 1.4;                                                                  % no unit
lambda = 0.0064;                                                            % mJ/mm
M0 = 0.4;                                                                   % mm^5/mJs
L = 20;                                                                     % mm
cl = 0.18;                                                                  % no unit, initial condition for <L/2
cu = 0.82;                                                                  % no unit, initial condition for >=L/2
si = [0.339981043584856,-0.339981043584856,0.861136311594053,-0.861136311594053];        % xi values for Galerkin method for n=4
wi = [0.652145454862546,0.652145454862546,0.347854845137454,0.347854845137454];          % wi values for Galerkin method for n=4
c_tilde = 0;                                                                % no unit
c_u = 1;                                                                    % no unit
c_dash_half_length = 0.522913;                                              % mm^-1
s = (c_u - c_tilde) / c_dash_half_length;                                   % thickness of the interface
fprintf("READ THE DOCUMENTATION BEFORE RUNNING THE PROGRAM!!!")
n = input("Enter the number of elements (to satisfy heuristic condition, n >= 105): ");
del_t = input("Enter the time step value: ");
t_tot = input("Enter the total time in seconds: ");
syms xi                                                                     % Declaring xi as a symbol to be easily parsed
he = L/n;                                                                   % Compute the heuristic value according to the number of elements
c = zeros(2*(n+1),1);                                                       % global c

for d = 0:n-1                                                               % Building the c global as per the condition given
    if (d*he) < L/2
        c((2*d)+1,1) = cl;                                                  
        c((2*d)+3,1) = cl;
    elseif (n*he) >= L/2
        c((2*d)+1,1) = cu;
        c((2*d)+3,1) = cu;
    end     % end of conditional statement
end         % end of element-wise loop
cm = c;                                                                     % initial condition for first time step
cmplus1 = cm;                                                               % initial condition for first time step
counter = 0;                                                                % counter initialization
%% Main function
% Looping for each time step 
for z = 0:del_t:t_tot                                                       % c = index of time step
    counter =  counter + 1;                                                 % counter is started to identify the Newton-Raphson iteration
    dummy = 0;                                                              % initialising a temporary variable to be used later
    NR_count = 0;                                                           % initialising a variable to store the total number of Newton-Raphson iterations run

    while true
        NR_count = NR_count + 1;
        R_global = zeros(2*(n+1),1);
        K_global = zeros(2*(n+1),2*(n+1));
        F_global = zeros(2*(n+1),1);
        
        % Looping through each element
        for b = 0:n-1                                                       % b = index of element
            xe = [b*he,(b+1)*he]';                                          % compute xe as it varies only per element
            J = [-0.5, 0.5] * xe;                                           % compute J as it varies only per element
            A = zeros(4,2*(n+1));                                           % Initialize the global assembly matrix
            A(1,(2*b)+1) = 1;                                               % Assigning 1 to the corresponding index according to elemental c
            A(2,(2*b)+2) = 1;                                               % Assigning 1 to the corresponding index according to elemental c
            A(3,(2*b)+3) = 1;                                               % Assigning 1 to the corresponding index according to elemental c
            A(4,(2*b)+4) = 1;                                               % Assigning 1 to the corresponding index according to elemental c
            % Initialising the elemental matrices
            R_element = zeros(4,1);
            K_element = zeros(4,1);
            F_element = zeros(4,1);
            
            % Galerkin method - element routine
            for a = 1:4                                                     % a = index of Galerkin's method
               cmplus1_e = A * cmplus1;                                     % Compute the elemental cmplus1
               cm_e = A * cm;                                               % Compute the elemental cm
               [rtemp,ktemp,ftemp] = residue_stiffness_fint(lambda,si(a),J,he,M0,psi0,chi,del_t,cmplus1_e,cm_e);        % Function call to compute the residual, tangent stiffness matrices and internal force vector
               r = rtemp*wi(a);                                             % multiplying by the weight and storing it in a temporary variable
               k = ktemp*wi(a);                                             % multiplying by the weight and storing it in a temporary variable
               f = ftemp*wi(a);                                             % multiplying by the weight and storing it in a temporary variable
               R_element = R_element + r;                                   % incrementing for each iteration of Galerkin method
               K_element = K_element + k;                                   % incrementing for each iteration of Galerkin method
               F_element = F_element + f;                                   % incrementing for each iteration of Galerkin method
            end         % end of Galerkin loop

            % Material routine - formation of Global terms
            rgtemp = A' * R_element;                                        % Compute the R for b element
            R_global = R_global + rgtemp;                                   % Increment the temporary variable for each element
            kgtemp = A' * K_element * A;                                    % Compute the K for b element
            K_global = K_global + kgtemp;                                   % Increment the temporary variable for each element
            fgtemp = A' * F_element;                                        % Compute the Fint for b element
            F_global = F_global + fgtemp;                                   % Increment the temporary variable for each element
        end             % end of element loop
        
        del_c = K_global\(-R_global);                                       % inverse computation using backslash operator        
        cmplus1 = cmplus1 + del_c;                                          % Computing the global cmplus1 to be used in the next time step
        
        if (NR_count == 1) && (z == 0)                                      % Store the first value of the internal force vector to be used in the convergence condition 
            prev_F_global = F_global;
        end

        % Displaying the convergence terms to visually check if they are indeed converging
        u = norm(del_c,inf);
        v = (1e-5)*(norm(cmplus1,inf));
        disp(u);        disp(v);
        w = norm(R_global,inf);
        x = 0.005*max(norm(prev_F_global,inf),1e-8);
        disp(w);        disp(x)
       
        if  ( (u < v) && (w < x) ) == true                                  % convergence criteria
            break
        end             % end of convergence condition handle
    end                 % end of Newton-Raphson loop 
    cm = cmplus1;                                                           % reassign the value of the variable to be used in the next time step
    
    for h = 0:n-1                                                           % Element routine for total energy
       [pelement] = stat_solution(psi0,chi,lambda,si,wi,he,n,cmplus1,h);
       dummy = dummy + pelement;                                            % Incrementing the variable w.r.t. previous value
    end                 % end of elementwise loop for total energy
    
    pglobal(counter) = dummy;                                               % store the value of the total energy with time step counter's index
    NR(counter) = NR_count;                                                 % store the number of the Newton-Raphson interations with time step counter's index
end                     % end of time step loop
%% Conclusion
% Plotting the desired values
subplot(3,1,1);
plot(cmplus1(1:2:end)); title("Concentration variation"); xlabel("Number of elements"); ylabel("Concentration")     % Plotting the concentration
subplot(3,1,2);                                                                     
plot(0:del_t:t_tot,pglobal); title("Total energy variation"); xlabel("Time (seconds)"); ylabel("Total energy")                % Plotting the total energy
subplot(3,1,3);
plot(0:del_t:t_tot,NR); title("Newton-Raphson visualisation"); xlabel("Time (seconds)"); ylabel("Iterations")   % Plotting the Newton-Raphson iterations
q = gcf; q.Position = [100 100 1280 720];                                            % Changing the size of the plot
%% Functions
% Function to compute R and K for each element
function [R_element,K_element,F_int] = residue_stiffness_fint(lambda,xi,J,he,M0,psi0,chi,del_t,cmplus1_e,cm_e)
    
    Nherm = [(1/4)*((xi^3)-(3*xi)+2), ...                                   
                (he/8)*((xi^3)-(xi^2)-xi+1), ...
                    (1/4)*(-(xi^3)+(3*xi)+2), ...
                        (he/8)*((xi^3)+(xi^2)-xi+1)];                       % Compute the Hermite shape function
    
    c_dot_e = cmplus1_e - cm_e;                                             % Computing the difference between cmplus1 and cm to be used in the residual equation
    c = Nherm*cmplus1_e;                                                    % Computing the global c matrix
    M = 4 * M0 * ( (c^2) * ((1-c)^2) );                                     % Compute the mobility function
    M1d = 8 * M0 * ( (2*(c^3)) - (3*(c^2)) + c );                           % Compute the first derivative of mobility function
    M2d = 8 * M0 * ( (6*(c^2))- (6*c) + 1 );                                % Compute the second derivative of mobility function
    psi2d = 2*psi0*chi*( (6*(c^2)) - (6*c) + 1 );                           % Compute the second derivative of the free energy function
    psi3d = 2*psi0*chi*(12*c - 6);                                          % Compute the third derivative of the free energy function

    B = [(3/4)*( (xi^2) - 1 ), ... 
            (he/8)*( (3*(xi^2))- (2*xi) + 1 ), ... 
                (-3/4)*( (xi^2)-1 ), ...
                    (he/8)*( (3*(xi^2)) + (2*xi) - 1 )] * inv(J);           % Compute B matrix  
    G = [(3/2)*xi, ...
            (he/4)*( (3*xi) - 1 ), ...
                (-3/2)*xi, ...
                    (he/4)*( (3*xi) + 1 )] * (inv(J)^2);                    % Compute the G matrix
  
    c1d = B*cmplus1_e;                                                      % c dash
    c2d = G*cmplus1_e;                                                      % c double dash

    R_element = (he/2) * (((Nherm'*(Nherm*c_dot_e)) / del_t) + (M*psi2d*c1d*B') + (M1d*lambda*c2d*c1d*B') + (M*lambda*c2d*G')); % Compute the residual matrix for the element
    
    % Computing each element of the tangent stiffness matrix individually and then assembling them into a single equation
    a1 = (M1d * psi2d * c1d) * (B' * Nherm);
    a2 = (M * psi3d * c1d) * (B' * Nherm);
    a3 = (M * psi2d) * (B' * B);
    b1 = (M2d * c2d * c1d * lambda) * (B' * Nherm);
    b2 = (M1d * c1d * lambda) * (B' * G);
    b3 = (M1d * c2d * lambda) * (B' * B);
    c1 = (M1d * lambda * c2d) * (G' * Nherm);
    c2 = (M * lambda) * (G' * G);
    z1 = (Nherm' * Nherm) / del_t; 

    K_element = (he/2) * (z1 + a1 + a2 + a3 + b1 + b2 + b3 + c1 + c2);      % Compute the tangent stiffness matrix
    F_int = (he/2) * ((M*psi2d*c1d*B') + (M1d*lambda*c2d*c1d*B') + (M*lambda*c2d*G'));  % Compute the internal force vector
end

% Function to compute the stationary solution
function [pelement] = stat_solution(psi0,chi,lambda,si,wi,he,n,cmplus1,h)
    xe = [h*he,(h+1)*he]';                                                  % compute xe as it varies only per element
    J = [-0.5, 0.5] * xe;                                                   % compute the Jacobian
     % Initialize the global assembly matrix
        A1 = zeros(4,2*(n+1));
        A1(1,(2*h)+1) = 1;                                                  % Assigning 1 to the corresponding index according to elemental c
        A1(2,(2*h)+2) = 1;                                                  % Assigning 1 to the corresponding index according to elemental c
        A1(3,(2*h)+3) = 1;                                                  % Assigning 1 to the corresponding index according to elemental c
        A1(4,(2*h)+4) = 1;                                                  % Assigning 1 to the corresponding index according to elemental c
        ce = A1 * cmplus1;                                                  % Compute the elemental c w.r.t. cmplus1
        pelement = 0;                                                       % Initialize a counter
        
        for g = 1:4
            Nherm = [(1/4)*((si(g)^3)-(3*si(g))+2), ...                     
                        (he/8)*((si(g)^3)-(si(g)^2)-si(g)+1), ...
                            (1/4)*(-(si(g)^3)+(3*si(g))+2), ...
                                (he/8)*((si(g)^3)+(si(g)^2)-si(g)+1)];      % Compute the Hermite shape function
            c = Nherm * ce;
            psi = psi0 * chi * (c^2)*((1-c)^2);                             % Compute the free energy function
            B = [(3/4)*( (si(g)^2) - 1 ), ... 
                    (he/8)*( (3*(si(g)^2))- (2*si(g)) + 1 ), ... 
                        (-3/4)*( (si(g)^2)-1 ), ...
                            (he/8)*( (3*(si(g)^2)) + (2*si(g)) - 1 )] * inv(J);           % Compute B matrix  
            psi_temp = det(J) * (psi + (0.5*lambda*abs((B*ce)^2)));
            pgal = wi(g) * psi_temp;                                        % Multiplying by the corresponding weight 
            pelement = pelement + pgal;                                     % Incrementing the result of each Gauss point
        end
end