function [err_state_matrix,for_plot_covariance]=ass3_q3_kf
% This is a KF implementation
% Have invoked a GPS sensor in this. In continuation of ass3_q3_q2
% measured position=actual_position+measurement_noise
% Here difference between actual and estimated state is plotted & it is seen clearly that it eventually becomes 0
delt=0.01;
final_time=500;
initial_state=[0;0;0];
err_state_matrix=zeros(3,(final_time/delt)+1);
err_state_matrix(:,1)=initial_state; % at time 0
phi=[0 1 0;0 0 -1;0 0 0];
B=[0;0;1];
C=[0;-1;0];

noise_b=normrnd(0,10^(-3),[1,final_time/delt+1]);
noise_gamma=normrnd(0,0.05,[1,final_time/delt+1]);
noise_position=normrnd(0,sqrt(3),[1,final_time/delt+1]);

bias=[0 0 0;0 0 0;0 0 1]; % B*wb*wb'*B'
gamma=[0 0 0;0 1 0;0 0 0]; % C*gamma*gamma'*C'
H=[1 0 0];

covariance_matrix=zeros(3,3);
for_plot_covariance=zeros(3,(final_time/delt)+1);
counter=1;
iden=[1 0 0;0 1 0;0 0 1];
phi_discrete=(phi*delt+iden);

% Prediction Step
for i=2:(final_time/delt+1)
  
    col_no=i;
    
    err_state_matrix(:,col_no)=(phi_discrete*err_state_matrix(:,col_no-1))+(B*delt*noise_b(col_no-1))+(C*delt*noise_gamma(col_no-1));
    
    covariance_matrix=(phi_discrete*covariance_matrix*phi_discrete')+(bias)*(noise_b(col_no-1)*delt)^2+(gamma)*(noise_gamma(col_no-1)*delt)^2;
    
    
    % GPS measurement at 1 second interval. Correction Step
    if col_no==(100*counter+1)
        counter=counter+1;
        Sk=H*covariance_matrix*H'+((noise_position(i))^2);
        Lk=(covariance_matrix*H')*inv(Sk);
        covariance_matrix=covariance_matrix-Lk*Sk*Lk';
        err_state_matrix(:,col_no)=(iden-Lk*H)*err_state_matrix(:,col_no)-Lk*(noise_position(i));
    end
    
    for_plot_covariance(:,col_no)=[covariance_matrix(1,1);covariance_matrix(2,2);covariance_matrix(3,3)];
        
end
plot(0:0.01:final_time,err_state_matrix(1,1:(final_time/delt+1)),'r')
xlabel('Time in seconds')
ylabel('Error')
title('Position Error (INS+GPS)')
grid on
