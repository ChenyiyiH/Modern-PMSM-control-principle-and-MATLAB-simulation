function[sys,x0,str,ts] = EKF(t,x,u,flag)
switch flag,
    % 初始化
    case 0,
        [sys,x0,str,ts] = mdlInitializeSizes;
    % 微分

    % 更新
    case 2,
        sys = mdlUpdate(t,x,u);
    % 输出
    case 3,
        sys = mdlOutputs(t,x,u);
    case {1,4,9}
        sys = [];
    otherwise
        error(['Unhandled flag = ',num2str(flag)]);
end
% end sfuntmpl
% mdllnitializeSizes
% Return the sizes, initial conditions, and sample times for the S - function
function [sys,x0,str,ts] = mdlInitializeSizes
global P0;
sizes = simsizes;

sizes.NumContStates = 0;
sizes.NumDiscStates = 4;
sizes.NumOutputs = 4;
sizes.NumInputs = 4;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0 = [0 0 0 0];
P0 = diag([0.1 0.1 0 0]);%P0的初始估计值

str = [];
ts = 1e-6;
function sys = mdlUpdate(t,x,u)
global P0;
Rs = 2.875;
Ls = 0.0085;
np = 4;
J = 0.001;
flux = 0.3;
B = 0;
Q = diag([0.1 0.1 1 0.01]);
R = diag([0.2 0.2]);
T = 1e-6;
vs_ab = [u(1) u(2)]';
is_ab = [u(3) u(4)]';
H = [1 0 0 0;0 1 0 0];
B = [1/Ls 0 0 0;0 1/Ls 0 0]';
F = [-Rs/Ls 0 flux/Ls*sin(x(4)) flux/Ls*x(3)*cos(x(4));...
     0 -Rs/Ls  -flux/Ls*cos(x(4)) flux/Ls*x(3)*sin(x(4));...
     0 0 0 0;...
     0 0 1 0];
f1 = [-Rs/Ls*x(1) + flux/Ls*x(3)*sin(x(4));...
      -Rs/Ls*x(2) - flux/Ls*x(3)*cos(x(4));...
      0;x(3)];
f2 = diag([1 1 1 1])+T*F;
x_pre = x+T*(f1+B*vs_ab);
y_pre = H*x_pre;
y = is_ab;
P_pre = f2*P0*f2'+Q;
K = P_pre*H'*inv(H*P_pre*H'+R);
sys = x_pre+K*(y-y_pre);
P0 = P_pre-K*H*P_pre;

function sys = mdlOutputs(t,x,u)
sys = x;