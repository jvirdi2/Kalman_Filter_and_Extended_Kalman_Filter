function [err_state_matrix,for_plot_covariance]=ass3_q2
% Part (a), calculation of difference between state estimator and true
% state. The state vector is [p,v,b] where p is postion, v is velocity and b is bias
% pdot=v
% vdot=a
% acc_measurement=actual acc+bias+noise
% bias_dot=wb (zero mean Gaussian random walk)
% Plots difference between actual state and estimated state. It will diverge due to noise
delt=0.01;
initial_state_err=[0;0;0];
err_state_matrix=zeros(3,(500/delt)+1);
err_state_matrix(:,1)=initial_state_err; % at time 0
phi=[0 1 0;0 0 -1;0 0 0];
B=[0;0;1];
C=[0;-1;0];
noise_b=normrnd(0,10^(-3),[1,500/delt]);
noise_gamma=normrnd(0,0.05,[1,500/delt]);

% Part 2_covariance matrix
bias=[0 0 0;0 0 0;0 0 1]; % B*wb*wb'*B'
gamma=[0 0 0;0 1 0;0 0 0]; % C*gamma*gamma'*C'

covariance_matrix=zeros(3,3);
for_plot_covariance=zeros(3,(500/delt)+1);
iden=[1 0 0;0 1 0;0 0 1];
phi_discrete=(phi*delt+iden);
for i=2:50001
    col_no=i;
    err_state_matrix(:,col_no)=phi_discrete*err_state_matrix(:,col_no-1)+B*delt*noise_b(col_no-1)+C*delt*noise_gamma(col_no-1);
    
    covariance_matrix=(phi_discrete*covariance_matrix*phi_discrete')+(bias)*(noise_b(col_no-1)*delt)^2+(gamma)*(noise_gamma(col_no-1)*delt)^2;
    for_plot_covariance(:,col_no)=[covariance_matrix(1,1);covariance_matrix(2,2);covariance_matrix(3,3)];
end
plot(0:0.01:500,err_state_matrix(1,1:50001),'r')
xlabel('Time in seconds')
ylabel('Error')
title('Error in position with only INS')
%plot(0:0.01:500,err_state_matrix(2,1:50001),'k')
%plot(0:0.01:500,err_state_matrix(3,1:50001))
hold on
grid on
