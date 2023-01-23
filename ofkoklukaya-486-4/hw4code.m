%% init

clc
clear
close all
load hw2_data.mat;

q_bi=zeros(4,2500);
q_bi(:,1)=transpose([1 0 0 1])/sqrt(2);
w_bi=zeros(3,25000);
w_bi(:,1)=transpose([0.1 0.1 0.1]);
mag_body=zeros(3,25000);
M_b=zeros(3,25000);
J=[6.9,0,0;...
    0,7.5,0;...
    0,0,8.4];
J_inv=inv(J);
t=zeros(25000);
A_bi=zeros(3,3,25000);
dt=1;
k=(4*pi*(1+sin(deg2rad(98)))*6.9)/5800;

for i=1:25000
    t(i)=i;
end

%% magnetic bdot

limit=6000;

for i=1:limit


    

    q_bi(:,i)=[q_bi(4,i),q_bi(1,i),q_bi(2,i),q_bi(3,i)];
    A_bi(:,:,i)=quat2dcm(transpose(q_bi(:,i)));
    q_bi(:,i)=[q_bi(2,i),q_bi(3,i),q_bi(4,i),q_bi(1,i)];
    
    mag_body(:,i)=A_bi(:,:,i)*(mag_eci(:,i));

    M_b(:,i) = ( k / norm(mag_body(:,i))) *...
        cross(w_bi(:,i),mag_body(:,i)/norm(mag_body(:,i)));
    
    N_md=cross(M_b(:,i),mag_body(:,i));

    qk1=0.5*omgof(w_bi(:,i))*q_bi(:,i); %rk4
    wk1=J_inv*(N_md-(cross(w_bi(:,i),(J*w_bi(:,i)))));
    qrk=q_bi(:,i)+qk1*(dt/2);
    wrk=w_bi(:,i)+wk1*(dt/2);
    qk2=0.5*omgof(wrk)*qrk;
    wk2=J_inv*(N_md-(cross(wrk,(J*wrk))));
    qrk=q_bi(:,i)+qk2*(dt/2);
    wrk=w_bi(:,i)+wk2*(dt/2);
    qk3=0.5*omgof(wrk)*qrk;
    wk3=J_inv*(N_md-(cross(wrk,(J*wrk))));
    qrk=q_bi(:,i)+qk3*(dt);
    wrk=w_bi(:,i)+wk3*(dt);
    qk4=0.5*omgof(wrk)*qrk;
    wk4=J_inv*(N_md-(cross(wrk,(J*wrk))));
    q_bi(:,i+1)=q_bi(:,i)+((qk1+2*qk2+2*qk3+qk4)/6)*dt;
    w_bi(:,i+1)=w_bi(:,i)+((wk1+2*wk2+2*wk3+wk4)/6)*dt;

    q_bi(:,i+1)=q_bi(:,i+1)/norm(q_bi(:,i+1));
end

figure(1);
subplot(3,2,1)
plot(t(1:limit),w_bi(1,1:limit));
title("x-components of the angular rate in the body frame");
ylabel("w_{b_x}");xlabel("t(s)");
subplot(3,2,3)
plot(t(1:limit),w_bi(2,1:limit));
title("y-components of the angular rate in the body frame");
ylabel("w_{b_y}");xlabel("t(s)");
subplot(3,2,5)
plot(t(1:limit),w_bi(3,1:limit));
title("z-components of the angular rate in the body frame");
ylabel("w_{b_z}");xlabel("t(s)");

subplot(3,2,2)
plot(t(1:limit),M_b(1,1:limit));
title("x-components of the dipole moment in the body frame");
ylabel("M_{d_x}");xlabel("t(s)");
subplot(3,2,4)
plot(t(1:limit),M_b(2,1:limit));
title("y-components of the dipole moment in the body frame");
ylabel("M_{d_y}");xlabel("t(s)");
subplot(3,2,6)
plot(t(1:limit),M_b(3,1:limit));
title("z-components of the dipole moment in the body frame");
ylabel("M_{d_z}");xlabel("t(s)");

%% Q1-2

limit=18000;
M_bmax=3*1e-9;

for i=1:limit

    q_bi(:,i)=[q_bi(4,i),q_bi(1,i),q_bi(2,i),q_bi(3,i)];
    A_bi(:,:,i)=quat2dcm(transpose(q_bi(:,i)));
    q_bi(:,i)=[q_bi(2,i),q_bi(3,i),q_bi(4,i),q_bi(1,i)];
    
    mag_body(:,i)=A_bi(:,:,i)*(mag_eci(:,i));

    M_b(:,i) = ( k / norm(mag_body(:,i))) *...
        cross(w_bi(:,i),mag_body(:,i)/norm(mag_body(:,i)));
    
    for j=1:3
        if abs(M_b(j,i))>M_bmax
            M_b(j,i) =  sign(M_b(j,i))*M_bmax;
        end
    end

    N_md=cross(M_b(:,i),mag_body(:,i));

    qk1=0.5*omgof(w_bi(:,i))*q_bi(:,i); %rk4
    wk1=J_inv*(N_md-(cross(w_bi(:,i),(J*w_bi(:,i)))));
    qrk=q_bi(:,i)+qk1*(dt/2);
    wrk=w_bi(:,i)+wk1*(dt/2);
    qk2=0.5*omgof(wrk)*qrk;
    wk2=J_inv*(N_md-(cross(wrk,(J*wrk))));
    qrk=q_bi(:,i)+qk2*(dt/2);
    wrk=w_bi(:,i)+wk2*(dt/2);
    qk3=0.5*omgof(wrk)*qrk;
    wk3=J_inv*(N_md-(cross(wrk,(J*wrk))));
    qrk=q_bi(:,i)+qk3*(dt);
    wrk=w_bi(:,i)+wk3*(dt);
    qk4=0.5*omgof(wrk)*qrk;
    wk4=J_inv*(N_md-(cross(wrk,(J*wrk))));
    q_bi(:,i+1)=q_bi(:,i)+((qk1+2*qk2+2*qk3+qk4)/6)*dt;
    w_bi(:,i+1)=w_bi(:,i)+((wk1+2*wk2+2*wk3+wk4)/6)*dt;

    q_bi(:,i+1)=q_bi(:,i+1)/norm(q_bi(:,i+1));
end

figure(2);
subplot(3,2,1)
plot(t(1:limit),w_bi(1,1:limit));
title("x-components of the angular rate in the body frame");
ylabel("w_{b_x}");xlabel("t(s)");
subplot(3,2,3)
plot(t(1:limit),w_bi(2,1:limit));
title("y-components of the angular rate in the body frame");
ylabel("w_{b_y}");xlabel("t(s)");
subplot(3,2,5)
plot(t(1:limit),w_bi(3,1:limit));
title("z-components of the angular rate in the body frame");
ylabel("w_{b_z}");xlabel("t(s)");

subplot(3,2,2)
plot(t(1:limit),M_b(1,1:limit));
title("x-components of the dipole moment in the body frame");
ylabel("M_{d_x}");xlabel("t(s)");
ylim([-3.5 3.5]*1e-9)
subplot(3,2,4)
plot(t(1:limit),M_b(2,1:limit));
title("y-components of the dipole moment in the body frame");
ylabel("M_{d_y}");xlabel("t(s)");
ylim([-3.5 3.5]*1e-9)
subplot(3,2,6)
plot(t(1:limit),M_b(3,1:limit));
title("z-components of the dipole moment in the body frame");
ylabel("M_{d_z}");xlabel("t(s)");
ylim([-3.5 3.5]*1e-9)


%% Q21

limit=300;
h_w=zeros(3,25000);
h_w(:,1)=0;
q_bi(:,1)=transpose([0.6853 0.6953 0.1531 0.1531]);
w_bi(:,1)=transpose([deg2rad(0.53) deg2rad(0.53) deg2rad(0.053)]);
q_cio=quatinv([1 0 0 0]);%inverse of q_c in other notation
dq=zeros(4,25000);
N_w=zeros(3,25000);

kp=50/1000;%approximated values based on text book example
kd=500/1000;

for i=1:limit
    
    q_bi(:,i)=[q_bi(4,i),q_bi(1,i),q_bi(2,i),q_bi(3,i)];
    dq(:,i)=transpose(quatmultiply(transpose(q_bi(:,i)),q_cio));
    q_bi(:,i)=[q_bi(2,i),q_bi(3,i),q_bi(4,i),q_bi(1,i)];
    dq(:,i)=[dq(2,i),dq(3,i),dq(4,i),dq(1,i)];

    N_w(:,i)= -kp * dq(1:3,i) - kd * w_bi(:,i);%eqn(7.7)
    
    qk1=0.5*omgof(w_bi(:,i))*q_bi(:,i); %rk4
    wk1=J_inv*(N_w(:,i)-(cross(w_bi(:,i),(J*w_bi(:,i)))));
    qrk=q_bi(:,i)+qk1*(dt/2);
    wrk=w_bi(:,i)+wk1*(dt/2);
    qk2=0.5*omgof(wrk)*qrk;
    wk2=J_inv*(N_w(:,i)-(cross(wrk,(J*wrk))));
    qrk=q_bi(:,i)+qk2*(dt/2);
    wrk=w_bi(:,i)+wk2*(dt/2);
    qk3=0.5*omgof(wrk)*qrk;
    wk3=J_inv*(N_w(:,i)-(cross(wrk,(J*wrk))));
    qrk=q_bi(:,i)+qk3*(dt);
    wrk=w_bi(:,i)+wk3*(dt);
    qk4=0.5*omgof(wrk)*qrk;
    wk4=J_inv*(N_w(:,i)-(cross(wrk,(J*wrk))));
    q_bi(:,i+1)=q_bi(:,i)+((qk1+2*qk2+2*qk3+qk4)/6)*dt;
    w_bi(:,i+1)=w_bi(:,i)+((wk1+2*wk2+2*wk3+wk4)/6)*dt;
    
    q_bi(:,i+1)=q_bi(:,i+1)/norm(q_bi(:,i+1));

end

figure(3);
subplot(4,1,1)
plot(t(1:limit),dq(1,1:limit));
ylabel("dq_1");xlabel("t(s)");
title("x-component of the quaternion vectors");
subplot(4,1,2)
plot(t(1:limit),dq(2,1:limit));
ylabel("dq_2");xlabel("t(s)");
title("y-component of the quaternion vectors");
subplot(4,1,3)
plot(t(1:limit),dq(3,1:limit));
ylabel("dq_3");xlabel("t(s)");
title("z-component of the quaternion vectors");
subplot(4,1,4)
plot(t(1:limit),dq(4,1:limit));
ylabel("dq_4");xlabel("t(s)");
title("scalar component of the quaternion vectors");
ylim([0 1.2]);


figure(4);
subplot(3,1,1)
plot(t(1:limit),N_w(1,1:limit));
ylabel("N_1");xlabel("t(s)");
title("x-component of the wheel torques");
subplot(3,1,2)
plot(t(1:limit),N_w(2,1:limit));
ylabel("N_2");xlabel("t(s)");
title("y-component of the wheel torques");
subplot(3,1,3)
plot(t(1:limit),N_w(3,1:limit));
ylabel("N_3");xlabel("t(s)");
title("z-component of the wheel torques");


%% Q22

limit=300;
h_w=zeros(3,25000);
h_w(:,1)=0;
q_bi(:,1)=transpose([0.6853 0.6953 0.1531 0.1531]);
w_bi(:,1)=transpose([deg2rad(0.53) deg2rad(0.53) deg2rad(0.053)]);
q_cio=quatinv([1 0 0 0]);%inverse of q_c in other notation
dq=zeros(4,25000);
h_w4=zeros(4,25000);
N_w=zeros(3,25000);

kp=50/1000;
kd=500/1000;

mapW4=pinv([1,-1,0,0;...
        1,1,1,1;...
        0,0,1,-1]/sqrt(2));


for i=1:limit
    
    q_bi(:,i)=[q_bi(4,i),q_bi(1,i),q_bi(2,i),q_bi(3,i)];
    dq(:,i)=transpose(quatmultiply(transpose(q_bi(:,i)),q_cio));
    q_bi(:,i)=[q_bi(2,i),q_bi(3,i),q_bi(4,i),q_bi(1,i)];
    dq(:,i)=[dq(2,i),dq(3,i),dq(4,i),dq(1,i)];
    
    N_w(:,i)= -kp * sign(dq(4,i)) * dq(1:3,i) - kd * w_bi(:,i);%eqn(7.7)    
    
    qk1=0.5*omgof(w_bi(:,i))*q_bi(:,i); %rk4
    wk1=J_inv*(N_w(:,i)-(cross(w_bi(:,i),(J*w_bi(:,i)))));
    qrk=q_bi(:,i)+qk1*(dt/2);
    wrk=w_bi(:,i)+wk1*(dt/2);
    qk2=0.5*omgof(wrk)*qrk;
    wk2=J_inv*(N_w(:,i)-(cross(wrk,(J*wrk))));
    qrk=q_bi(:,i)+qk2*(dt/2);
    wrk=w_bi(:,i)+wk2*(dt/2);
    qk3=0.5*omgof(wrk)*qrk;
    wk3=J_inv*(N_w(:,i)-(cross(wrk,(J*wrk))));
    qrk=q_bi(:,i)+qk3*(dt);
    wrk=w_bi(:,i)+wk3*(dt);
    qk4=0.5*omgof(wrk)*qrk;
    wk4=J_inv*(N_w(:,i)-(cross(wrk,(J*wrk))));
    q_bi(:,i+1)=q_bi(:,i)+((qk1+2*qk2+2*qk3+qk4)/6)*dt;
    w_bi(:,i+1)=w_bi(:,i)+((wk1+2*wk2+2*wk3+wk4)/6)*dt;
    
    q_bi(:,i+1)=q_bi(:,i+1)/norm(q_bi(:,i+1));
    
    h_w4(:,i)=mapW4*h_w(:,i);
    hdot=-cross(w_bi(:,i),h_w(:,i))-N_w(:,i);
    h_w(:,i+1)=h_w(:,i)+hdot*dt;

end

figure(5);
subplot(4,1,1)
plot(t(1:limit),dq(1,1:limit));
ylabel("dq_1");xlabel("t(s)");
title("x-component of the quaternion vectors");
subplot(4,1,2)
plot(t(1:limit),dq(2,1:limit));
ylabel("dq_2");xlabel("t(s)");
title("y-component of the quaternion vectors");
subplot(4,1,3)
plot(t(1:limit),dq(3,1:limit));
ylabel("dq_3");xlabel("t(s)");
title("z-component of the quaternion vectors");
subplot(4,1,4)
plot(t(1:limit),dq(4,1:limit));
ylabel("dq_4");xlabel("t(s)");
title("scalar component of the quaternion vectors");
ylim([0 1.2]);


figure(6);
subplot(3,1,1)
plot(t(1:limit),N_w(1,1:limit));
ylabel("N_1");xlabel("t(s)");
title("x-component of the wheel torques");
subplot(3,1,2)
plot(t(1:limit),N_w(2,1:limit));
ylabel("N_2");xlabel("t(s)");
title("y-component of the wheel torques");
subplot(3,1,3)
plot(t(1:limit),N_w(3,1:limit));
ylabel("N_3");xlabel("t(s)");
title("z-component of the wheel torques");


figure(7);
subplot(3,1,1)
plot(t(1:limit),h_w(1,1:limit));
ylabel("h_1");xlabel("t(s)");
title("x-component of the wheel angular momentum");
subplot(3,1,2)
plot(t(1:limit),h_w(2,1:limit));
ylabel("h_2");xlabel("t(s)");
title("y-component of the wheel angular momentum");
subplot(3,1,3)
plot(t(1:limit),h_w(3,1:limit));
ylabel("h_3");xlabel("t(s)");
title("z-component of the wheel angular momentum");


figure(8);
subplot(4,1,1)
plot(t(1:limit),h_w4(1,1:limit));
ylabel("h_{w_1}");xlabel("t(s)");
title("1st wheel angular momentum");
subplot(4,1,2)
plot(t(1:limit),h_w4(2,1:limit));
ylabel("h_{w_2}");xlabel("t(s)");
title("2nd wheel angular momentum");
subplot(4,1,3)
plot(t(1:limit),h_w4(3,1:limit));
ylabel("h_{w_3}");xlabel("t(s)");
title("3rd wheel angular momentum");
subplot(4,1,4)
plot(t(1:limit),h_w4(4,1:limit));
ylabel("h_{w_4}");xlabel("t(s)");
title("4th wheel angular momentum");




function D = omgof(w)
    D=[-[0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0], w(:); -transpose(w(:)), 0];
end
