

% EKF GPS-INS Fusion for a 6 DOF drone



clear all
close all

% load a check file with the data
load check
% This loads a matrix called A (arbitrary name, nothing to do with the real
% A) which contains the data that you need. This data has sensor data and some
% gps data to correct our estimators

% initial position in x y and z
x=[0 0 0];

% bias values, these are accelerometer and gyroscope biases
bp= 0;%.54*pi/180;
bq=-12*pi/180;
br=-.1*pi/180;
bfx = 0;
bfy = 0;
bfz = 0;


% IMU location specifier
r_imu=[-.5/12 -3/12 1/12]'*0; %% I have set this to zero, for Bonus, you can include the effect of this
r_GPS=[1.5, 0 ,0 ]; % This is the location of the GPS wrt CG, this is very important
%rotation matrix ------------------------------------------------------
phi= x(1);
theta= x(2);
psi = x(3);

%roation matrix body to inertial
L_bi = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta); ...
    sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi)  sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta); ...
    cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi)  cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];


Rt2b=L_bi';

[U,S,V]=svd(Rt2b);
R = U*V';
if 1+R(1,1)+R(2,2)+R(3,3) > 0
    b(1,1)    = 0.5*sqrt(1+R(1,1)+R(2,2)+R(3,3));
    b(2,1)    = ((R(3,2)-R(2,3))/4)/b(1);
    b(3,1)    = (R(1,3)-R(3,1))/4/b(1);
    b(4,1)    = (R(2,1)-R(1,2))/4/b(1);
    b       = b/norm(b);    
else
    R
    error('R diagonal too negative.')
    b = zeros(4,1);
end

%b =[0 0 0 0]';

% set quats
%-----------------------------------------------------------------
q1=b(1);%the quaternions are called b1-b4 in the data file that you loaded
q2=b(2);
q3=b(3);
q4=b(4);

%initialize velocity
vx = 0;
vy = 0;
vz = 0;

%set sample time
dt = .02;
tf=size(A,2);


% initialize x hat
% Note carefull the order the states appear in, this can be arbitrary, but
% we must stick to it along the entire code
%      [x y z vx vy vz     quat          gyro-bias accl-bias]
xhat = [0 0 0 0 0 0 b(1) b(2) b(3) b(4) bp bq br bfx bfy bfz]';

% noise params process noise 
Q = diag([.1 .1 .1 .1 .1 .1 .8 .8 .8 .8 .0001 .0001 .0001 .0001 .0001 .0001]);
%Q=100*Q;
% Really makes it converge very well Q=100*Q;


% noise params, measurement noise
% measurements are GPS position and velocity and mag
R = diag([9 9 9 3 3 3]);
%R=diag([10 10 10 8 8 8]);
%R=diag([0.01 0.01 0.01 0.01 0.01 0.01]);
%R = diag([10 10 10 8 8 8]);
%R = 0.01*R;

% Initialize P, the covariance matrix
P = diag([30 30 30 3 3 3 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1]);
Pdot=P*0;
tic
for k = 1:tf
    time= (k-1)*dt;
    
    %  Streaming sensor measurements and adjust for bias
    % these come from the file that is loaded in the begining
    %%%% These are still measurements but now bias has been removed
    xhatold=xhat;
    p = (A(1,k)*pi/180-xhat(11));
    q = (A(2,k)*pi/180-xhat(12));
    r = A(3,k)*pi/180-xhat(13);
    fx = (A(4,k)-xhat(14));
    fy = (A(5,k)-xhat(15));
    fz = -A(6,k)-xhat(16);  
    
    
    
    % Raw sensor measurments for plotting
    p_raw = A(1,k)*pi/180;
    q_raw = A(2,k)*pi/180;
    r_raw = A(3,k)*pi/180;
    fx_raw = A(4,k);
    fy_raw = A(5,k);
    fz_raw = A(6,k);
    
    quat = [xhat(7) xhat(8) xhat(9) xhat(10)]';
    quatmag= sqrt(quat(1)^2+quat(2)^2+quat(3)^2+quat(4)^2);
    
    quat = quat/quatmag;

    q1 = quat(1);
    q2 = quat(2);
    q3 = quat(3);
    q4 = quat(4);

    L_bl = [q1^2+q2^2-q3^2-q4^2    2*(q2*q3+q1*q4)      2*(q2*q4-q1*q3)
        2*(q2*q3-q1*q4)    q1^2-q2^2+q3^2-q4^2    2*(q3*q4+q1*q2)
        2*(q2*q4+q1*q3)      2*(q3*q4-q1*q2)    q1^2-q2^2-q3^2+q4^2];
    L_lb = L_bl';
    
    
    %% Implement your code here: 
    
    %% Prediction step
    %First write out all the dots, e.g. pxdot, pydot, q1dot etc
     
    %Now integrate Euler Integration for Process Updates and Covariance Updates
    % Euler works fine

    %Remember again the state vector [ px py pz vx vy vz q1 q2 q3 q4 bp bq br bx by bz]
    
    % Now write out all the partials to compute the transition matrix Phi
    %delV/delQ
    vq_term_11=2*(q1*fx-q4*fy+q3*fz);
    vq_term_12=2*(q2*fx+q3*fy+q4*fz);
    vq_term_13=2*(-q3*fx+q2*fy+q1*fz);
    vq_term_14=2*(-q4*fx-q1*fy+q2*fz);
    
    vq_term_21=2*(q4*fx+q1*fy-q2*fz);
    vq_term_22=2*(q3*fx-q2*fy-q1*fz);
    vq_term_23=2*(q2*fx+q3*fy+q4*fz);
    vq_term_24=2*(q1*fx-q4*fy+q3*fz);
    
    vq_term_31=2*(-q3*fx+q2*fy+q1*fz);
    vq_term_32=2*(q4*fx+q1*fy-q2*fz);
    vq_term_33=2*(-q1*fx+q4*fy-q3*fz);
    vq_term_34=2*(q2*fx+q3*fy+q4*fz);
    
    delv_delq=[[vq_term_11 vq_term_12 vq_term_13 vq_term_14];[vq_term_21 vq_term_22 vq_term_23 vq_term_24];[vq_term_31 vq_term_32 vq_term_33 vq_term_34]];
    
    %delV/del_abias
    delv_delba=(-1)*(L_lb);
    
    %delQ/delQ
    delq_delq=(-0.5)*[[0 p q r];[-p 0 -r q];[-q r 0 -p];[-r -q p 0]];
    
 
    %delQ/del_gyrobias
    delq_delbw=0.5*([[q2 q3 q4];[-q1 q4 -q3];[-q4 -q1 q2];[q3 -q2 -q1]]);
     
      
    % Now assemble the Transition matrix
    I=eye(3,3);
    trans_row1=[zeros(3,3) I zeros(3,4) zeros(3,3) zeros(3,3)];
    trans_row2=[zeros(3,3) zeros(3,3) delv_delq zeros(3,3) delv_delba];
    trans_row3=[zeros(4,3) zeros(4,3) delq_delq delq_delbw zeros(4,3)];
    trans_row4=[zeros(3,3) zeros(3,3) zeros(3,4) zeros(3,3) zeros(3,3)];
    trans_row5=[zeros(3,3) zeros(3,3) zeros(3,4) zeros(3,3) zeros(3,3)];
    
    phi_matrix=[trans_row1;trans_row2;trans_row3;trans_row4;trans_row5];
    
    phi_matrix_for_predict=phi_matrix(1:10,1:10);
   
    % Prediction of next state
    xhat(1:10)=(expm(phi_matrix_for_predict*dt))*xhat(1:10)+[zeros(1,5) 32.2*dt zeros(1,4)]';
    quat_update=[xhat(7) xhat(8) xhat(9) xhat(10)]';
    quat_update_norm=norm(quat_update);
    xhat(7)=xhat(7)/quat_update_norm;
    xhat(8)=xhat(8)/quat_update_norm;
    xhat(9)=xhat(9)/quat_update_norm;
    xhat(10)=xhat(10)/quat_update_norm;
    
   % Propagate the error covariance matrix, I suggest using the continuous integration since Q, R are not discretized 
%   Pdot = Phi*P+P*Phi' + Q;
    %P1 = Pdot*dt;
    %P = P +Pdot*dt;
    Pdot=phi_matrix*P+P*phi_matrix'+Q;
    P1=Pdot*dt;
    P=P+P1;
    
    %% Correction step
    % Get your measurements, 3 positions and 3 velocities from GPS
    z =[ A(7,k) A(8,k) A(9,k) A(10,k) A(11,k) A(12,k)]'; % x y z vx vy vz
       
    % Write out the measurement matrix linearization to get H
    
    %del P/del q
    Hxq_row_1=[-r_GPS(1)*2*q1 -r_GPS(1)*2*q2 r_GPS(1)*2*q3 r_GPS(1)*2*q4];
    Hxq_row_2=[-r_GPS(1)*2*q4 -r_GPS(1)*2*q3 -r_GPS(1)*2*q2 -r_GPS(1)*2*q1];
    Hxq_row_3=[r_GPS(1)*2*q3 -r_GPS(1)*2*q4 r_GPS(1)*2*q1 -r_GPS(1)*2*q2];
    H_delp_delq=[Hxq_row_1;Hxq_row_2;Hxq_row_3];
    
    % del v/del q
    
    Hvq_row_1=[(r_GPS(1)*2*q3*q+r_GPS(1)*2*q4*r) (r_GPS(1)*2*q4*q-r_GPS(1)*2*q3*r) (r_GPS(1)*2*q1*q-r_GPS(1)*2*q2*r) (r_GPS(1)*2*q2*q+r_GPS(1)*2*q1*r)];
    Hvq_row_2=[(-r_GPS(1)*2*q2*q-r_GPS(1)*2*q1*r) (r_GPS(1)*2*q2*r-r_GPS(1)*2*q1*q) (r_GPS(1)*2*q4*q-r_GPS(1)*2*q3*r) (r_GPS(1)*2*q3*q+r_GPS(1)*2*q4*r)];
    Hvq_row_3=[(r_GPS(1)*2*q1*q-r_GPS(1)*2*q2*r) (-r_GPS(1)*2*q2*q-r_GPS(1)*2*q1*r) (-r_GPS(1)*2*q3*q-r_GPS(1)*2*q4*r) (r_GPS(1)*2*q4*q-r_GPS(1)*2*q3*r)];
    H_delv_delq=[Hvq_row_1;Hvq_row_2;Hvq_row_3];
    
    % Assemble H
    H=[[I zeros(3,3) H_delp_delq zeros(3,6)];[zeros(3,3) I H_delv_delq zeros(3,6)]];
    rank(obsv(phi_matrix,H))
    %Compute Kalman gain
    S=H*P*H'+R;
    K=P*H'*inv(S);
    % Perform xhat correction    xhat = xhat + K*(z - H*xhat);

    xhat=xhat+K*(z-H*xhat);
    
    P=(eye(16,16)-K*H)*P;
    
    % We are predicting state for next time instant. But in data we are
    % plotting at current time instant. Hence the previous state is stored
    % (which would be current) now
    quatprev=[xhatold(7) xhatold(8) xhatold(9) xhatold(10)]';
    quatprev=quatprev/norm(quatprev);
     
    
    %% Now let us do some book-keeping 
    % Get some Euler angles
    [phi(k),theta(k),psi(k)]=quat2euler(quatprev');
    phi(k)=phi(k)*180/pi;
    theta(k)=theta(k)*180/pi;
    psi(k)=psi(k)*180/pi;
    
    quat1 = A(13:16,k);
    [phi_raw(k),theta_raw(k),psi_raw(k)]=quat2euler(quat1');
    phi_raw(k)=phi_raw(k)*180/pi;
    theta_raw(k)=theta_raw(k)*180/pi;
    psi_raw(k)=psi_raw(k)*180/pi;
    xhatR(k,:)= [xhatold(1:6)' quatprev(1) quatprev(2) quatprev(3) quatprev(4) xhatold(11:16)'];
    P_R(k,:) = diag(P);
    z_R(k,:) = z;
    OMEGA_raw(k,:)=[p_raw,q_raw,r_raw]';
    OMEGA(k,:)=[p,q,r]';
    FX(k,:)=[fx_raw,fy_raw,fz_raw]';
  
end
toc
t = 1:(k);
figure(1)
plot(t,P_R(:,1:3),'LineWidth',2)
title('Covariance of Position')
legend('px','py','pz')
grid on

figure(2)
plot(t,P_R(:,4:6),'LineWidth',2)
legend('pxdot','pydot','pzdot')
title('Covariance of Velocities')
grid on

figure(3)
plot(t,P_R(:,7:10),'LineWidth',2)
title('Covariance of Quaternions')
legend('q0cov','q1cov','q2cov','q3cov')
grid on

figure(8)
plot(t,phi,t,theta,t,psi,t,phi_raw,'b:',t,theta_raw,'g:',t,psi_raw,'r:','LineWidth',2)
title('Phi, Theta, Psi')
legend('phi','theta','psi','phiraw', 'thetaraw', 'psiraw')
grid on

figure(4)
hold on
plot(t,xhatR(:,1:3),t,A(7:9,:),':','LineWidth',2)
%plot(t,z_R(:,1),'r')
title('Position')
legend('px','py','pz','pxgps','pygps','pzgps')
grid on

figure(5)
plot(t,xhatR(:,4:6),t,A(10:12,:),':','LineWidth',2)
title('vel x y z')
legend('vx','vy','vz','vxgps','vygps','vzgps')
grid on

figure(6)
plot(t,xhatR(:,7:10),t,A(13:16,:),':','LineWidth',2)
legend('qo','q1','q2','q3','qodata','q1data','q2data','q3data')
title('Quaternions')
grid on


figure(9)
plot(t,xhatR(:,11:16),'LineWidth',2)
title('Bias')
legend('bp','bq','br','bfx','bfy','bfz')
grid on

figure(7)
plot(t,OMEGA(:,1:3),'LineWidth',2)
title('OMEGA without bias')
legend('p','q','r')
grid on

figure(10)
plot(t,OMEGA_raw(:,1),t,OMEGA_raw(:,2),t,OMEGA_raw(:,3),'LineWidth',2)
title('OMEGA raw without Bias')
legend('p','q','r')
grid on

figure(11)
plot(t,FX(:,1),t,FX(:,2),t,FX(:,3),'LineWidth',2)
title('accelerometer')
legend('ax','ay','az')
grid on
