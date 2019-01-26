function [phi, theta, psi] = quat2euler(q)
% Takes quaternion input and outputs euler angles
% The sign convention is the aero space frame
[n, junk] = size(q);
m=eye(3);
phi=zeros(n,1);
theta=phi;
psi = theta;

for i=1:n,...
      q_1 = q(i,1);
      q_2 = q(i,2);
      q_3 = q(i,3);
      q_4 = q(i,4);
  m(1,1) = 1.0 - 2.0*( q_3*q_3 + q_4*q_4 );
  m(1,2) = 2.0*( q_2*q_3 - q_1*q_4 );
  m(1,3) = 2.0*( q_2*q_4 + q_1*q_3 );
  m(2,1) = 2.0*( q_2*q_3 + q_1*q_4 );
  m(2,2) = 1.0 - 2.0*( q_2*q_2 + q_4*q_4 );
  m(2,3) = 2.0*( q_3*q_4 - q_1*q_2 );
  m(3,1) = 2.0*( q_2*q_4 - q_1*q_3 );
  m(3,2) = 2.0*( q_3*q_4 + q_1*q_2 );
  m(3,3) = 1.0 - 2.0*( q_2*q_2 + q_3*q_3 );
    phi(i,1)   = atan2( m(3,2), m(3,3) );
    theta(i,1) = -asin( m(3,1) );
    psi(i,1)   = atan2( m(2,1), m(1,1) );
end
 
