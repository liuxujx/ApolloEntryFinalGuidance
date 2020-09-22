function RefTrajV1()
%RefTraj 生成参考轨迹（倾侧角为常值）
%   考虑地球J2项和自转
clc; clear all; close all;
R0=6378136.3; %m
g0=9.80665; %m/s^2
GM=g0*R0^2;
% Vscale=sqrt(R0*g0); %7908.4m/s
% tscale=sqrt(R0/g0); %806.329s
W0=7.2921151e-5; % 7.2921151e-5;
J2=1.08262668e-3; % 1.08262668e-3;
rho0=1.225;
hs=7200;
d2r=pi/180;
r2d=180/pi;
mass=6800;
CL  = 0.3892;%气动升力系数
CD  = 1.3479;%气动阻力系数
Sr   = 12.88 ; %航天器参考横截面积m2   12.88
% 初始状态
r0 = R0+80e3; %初始再入地心距 100
theta0 = 0*d2r; %初始经度（弧度）%更改时下面一处也要更改
phi0 = 0*d2r; %初始纬度
V0   = 7800; %10660; %m/s初始再入速度 7450
gamma0 = -1.8*d2r; %初始再入角
psi0 = 90*d2r; %初始航向角，正北顺时针为正 60
s0 = 0; %初始航程
%仿真参数
p=1;
state(p,:)=[r0 theta0 phi0 V0 gamma0 psi0 s0];
dt=1;
time(p,1)=0;
h_thres=10e3;
while 1
    r=state(p,1); V=state(p,4);%无量纲
    L=rho0*exp(-(r-R0)/hs)*(V)^2*Sr*CL/2;% 有量纲升力 %/mass无量纲升力加速度 L/m mass*g0
    D=rho0*exp(-(r-R0)/hs)*(V)^2*Sr*CD/2;% 有量纲阻力 %/mass无量纲阻力加速度 D/m mass*g0
    ctrl(p,1)=-50*sin(0.005*p)*d2r; % 控制量sigma -80 -50*d2r
    ctrl(p,2)=L;
    ctrl(p,3)=D;
    ctrl(p,4)=L/D*cos(ctrl(p,1));
    state_temp=RK4(@EntryDyna3,time(p,1),state(p,:),ctrl(p,1),dt);
    h=state_temp(1)-R0;%先尝试积分一个周期，达到开伞高度后停止
    if h<=h_thres
        break;
        %break就是直接跳出该层循环 continue就是直接进入该层循环的下一次迭代
        %return就是直接退出程序或函数返回 大概的关系为return>break>continue
    end
    state(p+1,:)=state_temp;
    time(p+1,1)=time(p,1)+dt;
    p=p+1;
end

figure;
plot(time,(state(:,1)-R0)/1000,'k','linewidth',1.5)
xlabel('Time (s)'); ylabel('Altitude (km)')
figure
plot(time,state(:,4)/1000,'k','linewidth',1.5)
xlabel('Time (s)'); ylabel('Speed (km/s)')
figure;
plot3(state(:,2)*r2d,state(:,3)*r2d,(state(:,1)-R0)/1000,'k','linewidth',1.5)
xlabel('Longitude (deg)'); ylabel('Latitude (deg)'); zlabel('Altithde (km)');


outp=[time, state,ctrl];
save RefTrajV1.txt outp -ascii;
end

function outp=EntryDyna3(t,inp,ctrl)
R0=6378136.3; %m
g0=9.80665; %m/s^2
GM=g0*R0^2;
% Vscale=sqrt(R0*g0); %7908.4m/s
% tscale=sqrt(R0/g0); %806.329s
W0=7.2921151e-5; % 7.2921151e-5;
J2=1.08262668e-3; % 1.08262668e-3;
rho0=1.225;
hs=7200;
d2r=pi/180;
r2d=180/pi;
mass=6800;
CL  = 0.3892;%气动升力系数
CD  = 1.3479;%气动阻力系数
Sr   = 12.88 ; %航天器参考横截面积m2   12.88

%无量纲参数
r     = inp(1); %地心距
theta = inp(2); %经度
phi   = inp(3); %纬度
V     = inp(4); %速度
gamma = inp(5); %航迹倾角 飞行路径角
psi   = inp(6); %航迹偏角 航迹方位角
s     = inp(7); %航程
sigma = ctrl; %倾侧角
%相关变量
h=r-R0;
rho=rho0*exp(-h/hs);
L=rho*(V)^2*Sr*CL/2; % /mass
D=rho*(V)^2*Sr*CD/2; % /mass

%有量纲动力学模型
rdot = V*sin(gamma);
thetadot = V*cos(gamma)*sin(psi)/(r*cos(phi));
phidot = V*cos(gamma)*cos(psi)/r;
Vdot = -D/mass-GM/r^2*sin(gamma);
gammadot = L*cos(sigma)/(mass*V)-(GM/(r^2*V)-V/r)*cos(gamma);
psidot = L*sin(sigma)/(mass*V*cos(gamma))+V*cos(gamma)*sin(psi)*tan(phi)/r;
sdot = V*cos(gamma);
outp=[rdot thetadot phidot Vdot gammadot psidot sdot];
end

function outp=RK4(func,t,inp,ctrl,dt)
k1=func(t,inp,ctrl);
k2=func(t+dt/2,inp+dt/2*k1,ctrl);
k3=func(t+dt/2,inp+dt/2*k2,ctrl);
k4=func(t+dt,inp+dt*k3,ctrl);
outp=inp+ dt/6*(k1 + 2*k2 + 2*k3 + k4);
end