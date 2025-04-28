clear all;

% system parameter
Bw = [-0.5 1;
      -1 0.2];
beta = 0.4;

% Upwind method
ds = 0.01;
s = 0:ds:1;
dt = 0.001;
tf = 50;
t = 0:dt:tf;

% Given the maximum velocity of production
V = [0.5 0;
    0 0.3];

% Initialize of system state
y = zeros(length(s), length(t),2);

% Define new state type matrix
A = eye(length(s))-diag(ones([1,length(s)-1]),-1);

% Define the disturbance channel
for i=1:length(s)
    for j=1:length(t)
        w(i,j,1) = (i-1)*ds*exp(-0.1*(j-1)*dt)*sin(2*(j-1)*dt);
        w(i,j,2) = cos((j-1)*dt)/((j-1)*dt+1);
    end
end

% Given state delay
tau1 = round((0.2+0.1.*sin(pi.*t))./dt);

% Given input delay
tau2 = round((0.3+0.15.*cos(pi.*t))./dt);

% Given the constant initial condition
y(:,1,1) = 2*sin(2*s);
y(:,1,2) = 3-cos((s+1)/2);

% Define the boundary condition
bl = zeros([1, length(t)]);

% Calculate feedback gain
[K1, K2, K3, K4] = cal_gain;

for n=1:length(t)-1
    if n-tau1(n)-tau2(n) >= 1
        theta1 = n-tau1(n)-tau2(n); % as the premise variable for control rule will contain the input delay
    else
        theta1 = 1;
    end

    if n-tau2(n) >= 1
        theta2 = n-tau2(n);
    else
        theta2 = 1;
    end

    % Fuzzy sets of y2 without delay but as the premise variable for control rule will contain the input delay
    M11 = (4-y( :,theta2,2))./8;
    M12 = (y( :,theta2,2)+4)./8;

    % Fuzzy sets of y2 with state delays and as the premise variable for control rule will contain the input delay
    M21 = (4-y( :,theta1,2))./8;
    M22 = (y( :,theta1,2)+4)./8;

    % Define membership functions
    Mu1 = M11.*M21;
    Mu2 = M11.*M22;
    Mu3 = M12.*M21;
    Mu4 = M12.*M22;
    
    h1 = Mu1./(Mu1+Mu2+Mu3+Mu4);
    h2 = Mu2./(Mu1+Mu2+Mu3+Mu4);
    h3 = Mu3./(Mu1+Mu2+Mu3+Mu4);
    h4 = Mu4./(Mu1+Mu2+Mu3+Mu4);

    % State equations of RMS
    f1 = -8.*y(:,n,1) + y(:,n,1).*y(:,n,2)./4 + sin(y(:,n,2)).*y(:,n,1) + 10.*y(:,n,2) - y(:,n,2).^2./2 ..., 
        + Bw(1,1).*w(:,n,1) + Bw(1,2).*w(:,n,2) ...,
        - beta.*y(:,theta1,1).*y(:,theta1,2)./4 ...,
        + (K1(1).*h1+K2(1).*h2+K3(1).*h3+K4(1).*h4).*y(:,theta2,1) + (K1(2).*h1+K2(2).*h2+K3(2).*h3+K4(2).*h4).*y(:,theta2,2);
        
    f2 = -8.*y(:,n,2) + y(:,n,2).^2./2 + cos(y(:,n,2)).*y(:,n,2)..., 
        + Bw(2,1).*w(:,n,1) +Bw(2,2).*w(:,n,2) ...,
        - beta.*y(:,theta1,1) - beta.*y(:,theta1,2) ...,
        + (K1(1).*h1+K2(1).*h2+K3(1).*h3+K4(1).*h4).*y(:,theta2,1) + (K1(2).*h1+K2(2).*h2+K3(2).*h3+K4(2).*h4).*y(:,theta2,2);
    
    % Calculate the velocity varies with the WIP in the factory production line
    W1 = ds/48*(17*y(1,n,1)+59*y(2,n,1)+43*y(3,n,1)+49*y(4,n,1)+48*sum(y(5:length(s)-4,n,1))+49*y(length(s)-3,n,1)+43*y(length(s)-2,n,1)+59*y(length(s)-1,n,1)+17*y(length(s),n,1));
    W2 = ds/48*(17*y(1,n,2)+59*y(2,n,2)+43*y(3,n,2)+49*y(4,n,2)+48*sum(y(5:length(s)-4,n,2))+49*y(length(s)-3,n,2)+43*y(length(s)-2,n,2)+59*y(length(s)-1,n,2)+17*y(length(s),n,2));
    v1 = V(1,1)/(1+W1+4);
    v2 = V(2,2)/(1+W2+4);

    % Calculate the next step of the error system
    y( :,n+1,1) = y( :,n,1) - v1*dt/ds*A*y( : ,n,1) + f1*dt;
    y( :,n+1,2) = y( :,n,2) - v2*dt/ds*A*y( : ,n,2) + f2*dt;

    y( 1,n+1,1)=bl(n+1);
    y( 1,n+1,2)=bl(n+1);
end

% Plot the state trajectories at s=1
figure;
plot(t,y(end,:,1),'r--')
hold on
plot(t,y(end,:,2),'b')
xlabel('Time $t$','interpreter','latex', 'Fontsize',15);
ylabel('State responses at $s=1$','interpreter','latex', 'Fontsize',15);
legend('y_1(1,t)','y_2(1,t)');
hold off;

% Plot the state trajectories of y1
figure;
surf(t,s,y(:,:,1));
shading interp;
ylabel('Stage $s$','interpreter','latex', 'Fontsize',15),
xlabel('Time $t$','interpreter','latex', 'Fontsize',15),
xlim([0 tf]);
zlabel('$y_1(s,t)$','interpreter','latex', 'Fontsize',15);
view(45,25);
colorbar;

% Plot the state trajectories of y2
figure;
surf(t,s,y(:,:,2));
shading interp;
ylabel('Stage $s$','interpreter','latex', 'Fontsize',15),
xlabel('Time $t$','interpreter','latex', 'Fontsize',15),
xlim([0 tf]);
zlabel('$y_2(s,t)$','interpreter','latex', 'Fontsize',15);
view(45,25);
colorbar;
