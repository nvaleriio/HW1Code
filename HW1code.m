clc
clear
tic   % Start timing entire MATLAB program

% Homework1

% Fixed inputs 
A=[1.4 0.485 0];
B=[1.67 0.99 0];
C=[0.255 1.035 0];
D=[0.285 0.055 0];
E=[0.195 2.54 0];
F=[-0.98 2.57 0];
G=[0.05 0.2 0];
H=[-1.71 4.26 0];

AB=norm(B-A);
BC=norm(C-B);
CD=norm(D-C);
DE=norm(E-D);
CE=norm(E-C);
EF=norm(F-E);
FG=norm(G-F);
FH=1.843;

% Input speed 
w_AB_val = 1.85753;   

% rad/s  (≈ 17.7 rpm)
% Initial link AB angle
initialAngle_AB = atan2(B(2)-A(2), B(1)-A(1));
if initialAngle_AB < 0
    angleAB_horizontal = 2*pi + initialAngle_AB;
else
    angleAB_horizontal = initialAngle_AB;
end

%First-position table
FirstPos = table( ...
    A(1),A(2), B(1),B(2), C(1),C(2), D(1),D(2), ...
    E(1),E(2), F(1),F(2), G(1),G(2), H(1),H(2), ...
    'VariableNames', {'Ax','Ay','Bx','By','Cx','Cy','Dx','Dy','Ex','Ey','Fx','Fy','Gx','Gy','Hx','Hy'});
writetable(FirstPos,'FirstPositionTable.xlsx');
disp('First-position coordinates:');
disp(FirstPos);

% Storage
theta_deg = 0:360;
% Paths
newB_x = zeros(1,361); newB_y = zeros(1,361);
newC_x = zeros(1,361); newC_y = zeros(1,361);
newE_x = zeros(1,361); newE_y = zeros(1,361);
newF_x = zeros(1,361); newF_y = zeros(1,361);
% Angular velocity & acceleration of links 
wBC_arr  = zeros(1,361); wDCE_arr = zeros(1,361);
wEF_arr  = zeros(1,361); wFGH_arr = zeros(1,361);
aBC_arr  = zeros(1,361); aDCE_arr = zeros(1,361);
aEF_arr  = zeros(1,361); aFGH_arr = zeros(1,361);
% COM acceleration magnitudes (for comparison w/ PMKS “COM accel”)
aS1 = zeros(1,361); aS2 = zeros(1,361); aS3 = zeros(1,361); aS4 = zeros(1,361); aS5 = zeros(1,361);
% Torque (dynamic)
NTorque = zeros(1,361);

%Loop 
for theta = 0:1:360

% Position analysis (circle intersections)
B_new = vpa(A + [AB*cos(angleAB_horizontal+deg2rad(theta)), AB*sin(angleAB_horizontal+deg2rad(theta)), 0]); % ;

[C_x,C_y] = circcirc(B_new(1),B_new(2),BC, D(1),D(2),CD);
if any(isnan(vpa(C_x))) || any(isnan(vpa(C_y)))
    fprintf('New position C cannot be determined at theta=%d\n',theta); break;
end
C_1=[C_x(1) C_y(1) 0]; C_2=[C_x(2) C_y(2) 0];
if norm(C_1-C) < norm(C_2-C), C_new=vpa(C_1); else, C_new=vpa(C_2); end

[E_x,E_y] = circcirc(C_new(1),C_new(2),CE, D(1),D(2),DE);
if any(isnan(vpa(E_x))) || any(isnan(vpa(E_y)))
    fprintf('New position E cannot be determined at theta=%d\n',theta); break;
end
E_1=[E_x(1) E_y(1) 0]; E_2=[E_x(2) E_y(2) 0];
if norm(E_1-E) < norm(E_2-E), E_new=vpa(E_1); else, E_new=vpa(E_2); end

[F_x,F_y] = circcirc(E_new(1),E_new(2),EF, G(1),G(2),FG);
if any(isnan(vpa(F_x))) || any(isnan(vpa(F_y)))
    fprintf('New position F cannot be determined at theta=%d\n',theta); break;
end
F_1=[F_x(1) F_y(1) 0]; F_2=[F_x(2) F_y(2) 0];
if norm(F_1-F) < norm(F_2-F), F_new=vpa(F_1); else, F_new=vpa(F_2); end

% Save paths
newB_x(theta+1)=B_new(1); newB_y(theta+1)=B_new(2);
newC_x(theta+1)=C_new(1); newC_y(theta+1)=C_new(2);
newE_x(theta+1)=E_new(1); newE_y(theta+1)=E_new(2);
newF_x(theta+1)=F_new(1); newF_y(theta+1)=F_new(2);

% Update state for this theta
B=B_new; C=C_new; E=E_new; F=F_new;
H=((F-G)/(norm(F-G)))*(norm(F-G)+1.843)+G;

%  Centers of mass 
S1=(A+B)/2; S2=(B+C)/2; S3=(C+D+E)/3; S4=(E+F)/2; S5=(H+G)/2;

%  Mass & weights (SW numbers) 
M_AB=0.44; M_BC=1.10; M_DCE=0.76; M_EF=0.12; M_FGH=2.02;
g=9.81;
W_AB=[0 -M_AB*g 0]; W_BC=[0 -M_BC*g 0]; W_DCE=[0 -M_DCE*g 0]; W_EF=[0 -M_EF*g 0]; W_FGH=[0 -M_FGH*g 0];

% Position vectors to COMs
S1_A=A-S1; S1_B=B-S1; S2_B=B-S2; S2_C=C-S2; S3_C=C-S3; S3_D=D-S3; S3_E=E-S3;
S4_E=E-S4; S4_F=F-S4; S5_F=F-S5; S5_G=G-S5; S5_H=H-S5;

%  STATICS 
syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy T
F_A=[FAx FAy 0]; F_B=[FBx FBy 0]; F_C=[FCx FCy 0]; F_D=[FDx FDy 0];
F_E=[FEx FEy 0]; F_F=[FFx FFy 0]; F_G=[FGx FGy 0]; Torque=[0 0 T];
AppliedLoad=[0 -200 0];

eqn1= F_A+F_B+W_AB ==0;
eqn2= cross(S1_A, F_A)+ cross(S1_B, F_B) + Torque==0;
eqn3= -F_B+F_C+W_BC==0;
eqn4= cross(S2_B, -F_B) +cross(S2_C, F_C) ==0;
eqn5= -F_C+F_D+F_E+W_DCE==0;
eqn6= cross(S3_C, -F_C) +cross(S3_D, F_D) +cross(S3_E, F_E)==0;
eqn7= -F_E+F_F+W_EF==0;
eqn8= cross(S4_E,-F_E)+cross(S4_F, F_F)==0;
eqn9= -F_F+F_G+W_FGH+AppliedLoad==0;
eqn10= cross(S5_F,-F_F)+ cross(S5_G,F_G)+ cross(S5_H, AppliedLoad)==0;

solution = solve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10], ...
    [FAx, FAy, FBx, FBy, FCx, FCy, FDx, FDy, FEx, FEy, FFx, FFy, FGx, FGy, T]); 

% VELOCITIES 
w_AB = [0 0 w_AB_val];

syms wBC wDCE
w_BC=[0 0 wBC]; w_DCE=[0 0 wDCE];

eqn11 = cross(w_AB, B-A)+ cross(w_BC, C-B) + cross(w_DCE, D-C)==0;
solutionVelocities = solve(eqn11,[wBC,wDCE]);
omega_BC = double(solutionVelocities.wBC);
omega_DCE = double(solutionVelocities.wDCE);

o_BC=[0 0 omega_BC]; o_DCE=[0 0 omega_DCE];

syms wEF wFGH
w_EF=[0 0 wEF]; w_FGH=[0 0 wFGH];
eqn12 = cross(o_DCE, E-D)+ cross(w_EF, F-E) + cross(w_FGH, G-F) == 0;
solutionSecondLoop = solve(eqn12,[wEF,wFGH]);
omega_EF  = double(solutionSecondLoop.wEF);
omega_FGH = double(solutionSecondLoop.wFGH);

omega_EF_vector = [0 0 omega_EF];
omega_FGH_vector= [0 0 omega_FGH];

% Save Angular velocities
wBC_arr(theta+1)=omega_BC; wDCE_arr(theta+1)=omega_DCE;
wEF_arr(theta+1)=omega_EF; wFGH_arr(theta+1)=omega_FGH;

%  ACCELERATIONS (angular)
alphaAB=[0 0 0];
syms aDCE aBC
a_BC=[0 0 aBC]; a_DCE=[0 0 aDCE];

a_ba = cross(alphaAB, B-A) + cross(w_AB, cross(w_AB, B-A));
a_cb = cross(a_BC, C-B) + cross(o_BC, cross(o_BC, C-B));
a_dc = cross(a_DCE, D-C) + cross(o_DCE, cross(o_DCE, D-C));
eqn13 = a_ba + a_cb + a_dc == 0;

solutionAcceleration = solve(eqn13,[aBC,aDCE]);
alpha_BC  = double(solutionAcceleration.aBC);
alpha_DCE = double(solutionAcceleration.aDCE);
alphaBC=[0 0 alpha_BC]; alphaDCE=[0 0 alpha_DCE];

syms aEF aFGH
a_EF=[0 0 aEF]; a_FGH=[0 0 aFGH];
a_ed = cross(alphaDCE, E-D)+ cross(o_DCE, cross(o_DCE, E-D));
a_fe = cross(a_EF, F-E)  + cross(omega_EF_vector,  cross(omega_EF_vector,  F-E));
a_gf = cross(a_FGH, G-F) + cross(omega_FGH_vector, cross(omega_FGH_vector, G-F));
eqn14 = a_ed + a_fe + a_gf == 0;

solutionAcceleration_secondLoop = solve(eqn14,[aEF,aFGH]);
alphaEF  = double(solutionAcceleration_secondLoop.aEF);
alphaFGH = double(solutionAcceleration_secondLoop.aFGH);

% Save angular accelerations
aBC_arr(theta+1)=alpha_BC; aDCE_arr(theta+1)=alpha_DCE;
aEF_arr(theta+1)=alphaEF;  aFGH_arr(theta+1)=alphaFGH;

% COM accelerations (magnitudes for plots / PMKS comparison) 
alphaEF_vector=[0 0 alphaEF]; alphaFGH_vector=[0 0 alphaFGH];
Accln_S1_A = cross(alphaAB, S1-A)+ cross(w_AB, cross(w_AB, S1-A));
Accln_S2_A = cross(alphaAB, B-A)+ cross(w_AB, cross(w_AB, B-A)) + cross([0 0 alpha_BC], S2-B) + cross(o_BC, cross(o_BC, S2-B));
Accln_S3_D = cross([0 0 alpha_DCE], S3-D)+ cross(o_DCE, cross(o_DCE, S3-D));
Accln_S4_D = cross([0 0 alpha_DCE], E-D) + cross(o_DCE, cross(o_DCE, E-D)) + cross(alphaEF_vector, S4-E) + cross(omega_EF_vector, cross(omega_EF_vector, S4-E));
Accln_S5_G = cross([0 0 alphaFGH], S5-G) + cross(omega_FGH_vector, cross(omega_FGH_vector, S5-G));

aS1(theta+1)=norm(Accln_S1_A);
aS2(theta+1)=norm(Accln_S2_A);
aS3(theta+1)=norm(Accln_S3_D);
aS4(theta+1)=norm(Accln_S4_D);
aS5(theta+1)=norm(Accln_S5_G);

% DYNAMICS (Newton s 2nd law) 
MMoI_AB = 0.001; MMoI_BC = 0.18; MMoI_DCE = 0.06; MMoI_EF = 0.01; MMoI_FGH = 1.12;
syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NT
F_A=[NFAx NFAy 0]; F_B=[NFBx NFBy 0]; F_C=[NFCx NFCy 0]; F_D=[NFDx NFDy 0];
F_E=[NFEx NFEy 0]; F_F=[NFFx NFFy 0]; F_G=[NFGx NFGy 0]; F_T=[0 0 NT];

eq1 = F_A + F_B + W_AB == M_AB*Accln_S1_A;
eq2 = cross(S1_A, F_A) + cross(S1_B, F_B) + F_T == MMoI_AB*alphaAB;
eq3 = -F_B + F_C + W_BC == M_BC*Accln_S2_A;
eq4 = cross(S2_B, -F_B) + cross(S2_C, F_C) == MMoI_BC*[0 0 alpha_BC];
eq5 = -F_C + F_D + F_E + W_DCE == M_DCE*Accln_S3_D;
eq6 = cross(S3_C, -F_C) + cross(S3_D, F_D) + cross(S3_E, F_E) == MMoI_DCE*[0 0 alpha_DCE];
eq7 = -F_E + F_F + W_EF == M_EF*Accln_S4_D;
eq8 = cross(S4_E, -F_E) + cross(S4_F, F_F) == MMoI_DCE*[0 0 alphaEF];   % kept as you had it
eq9 = -F_F + F_G + W_FGH + AppliedLoad == M_FGH*Accln_S5_G;
eq10= cross(S5_F, -F_F) + cross(S5_G, F_G) == MMoI_FGH*[0 0 alphaFGH];

solutionDyn = solve([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10], ...
    [NFAx,NFAy,NFBx,NFBy,NFCx,NFCy,NFDx,NFDy,NFEx,NFEy,NFFx,NFFy,NFGx,NFGy,NT]);

NTorque(theta+1) = double(solutionDyn.NT);

end 

% for theta

% Joint linear speeds from paths 
dt = (deg2rad(1))/w_AB_val; % time step per 1°
Bx=newB_x(:); By=newB_y(:);
Cx=newC_x(:); Cy=newC_y(:);
Ex=newE_x(:); Ey=newE_y(:);
Fx=newF_x(:); Fy=newF_y(:);
V_B = [0; hypot(diff(Bx), diff(By))/dt];
V_C = [0; hypot(diff(Cx), diff(Cy))/dt];
V_E = [0; hypot(diff(Ex), diff(Ey))/dt];
V_F = [0; hypot(diff(Fx), diff(Fy))/dt];

% Position sweep table & PMKS-MATLAB export 
Tpos = table(theta_deg.', newB_x.',newB_y.', newC_x.',newC_y.', newE_x.',newE_y.', newF_x.',newF_y.', ...
    'VariableNames', {'theta_deg','Bx','By','Cx','Cy','Ex','Ey','Fx','Fy'});
writetable(Tpos,'PositionSweep_BC_EF.xlsx');

Tout = table(theta_deg.', ...
    wBC_arr.', wDCE_arr.', wEF_arr.', wFGH_arr.', ...
    aBC_arr.', aDCE_arr.', aEF_arr.', aFGH_arr.', ...
    V_B, V_C, V_E, V_F, ...
    aS1.', aS2.', aS3.', aS4.', aS5.', ...
    NTorque.', ...
    'VariableNames', {'theta_deg','w_BC','w_DCE','w_EF','w_FGH','alpha_BC','alpha_DCE','alpha_EF','alpha_FGH', ...
                      'V_B','V_C','V_E','V_F','aS1','aS2','aS3','aS4','aS5','TorqueN'});
writetable(Tout,'Kinematics_For_Comparison.csv');

% Plots (labeled for rubric) 
% Paths
figure;
subplot(2,2,1); plot(newB_x,newB_y); grid on; axis equal; title('Joint B Path'); xlabel('x (m)'); ylabel('y (m)');
subplot(2,2,2); plot(newC_x,newC_y); grid on; axis equal; title('Joint C Path'); xlabel('x (m)'); ylabel('y (m)');
subplot(2,2,3); plot(newE_x,newE_y); grid on; axis equal; title('Joint E Path'); xlabel('x (m)'); ylabel('y (m)');
subplot(2,2,4); plot(newF_x,newF_y); grid on; axis equal; title('Joint F Path'); xlabel('x (m)'); ylabel('y (m)');

% Input torque vs theta
figure; plot(theta_deg, NTorque); grid on; xlim([0 360]);
title('Input Torque vs \theta (Dynamic)'); xlabel('\theta (deg)'); ylabel('Torque (N·m)');

% Angular velocities
figure; plot(theta_deg,wBC_arr, theta_deg,wDCE_arr, theta_deg,wEF_arr, theta_deg,wFGH_arr);
grid on; xlim([0 360]);
title('Angular Velocity of Links'); xlabel('\theta (deg)'); ylabel('\omega (rad/s)');
legend({'BC','DCE','EF','FGH'},'Location','best');

% Angular accelerations
figure; plot(theta_deg,aBC_arr, theta_deg,aDCE_arr, theta_deg,aEF_arr, theta_deg,aFGH_arr);
grid on; xlim([0 360]);
title('Angular Acceleration of Links'); xlabel('\theta (deg)'); ylabel('\alpha (rad/s^2)');
legend({'BC','DCE','EF','FGH'},'Location','best');

% Joint linear speeds
figure; plot(theta_deg,V_B, theta_deg,V_C, theta_deg,V_E, theta_deg,V_F);
grid on; xlim([0 360]);
title('Joint Linear Speed'); xlabel('\theta (deg)'); ylabel('Speed (m/s)');
legend({'B','C','E','F'},'Location','best');

% COM accelerations (magnitudes)
figure; plot(theta_deg,aS1, theta_deg,aS2, theta_deg,aS3, theta_deg,aS4, theta_deg,aS5);
grid on; xlim([0 360]);
title('Center-of-Mass Acceleration Magnitude'); xlabel('\theta (deg)'); ylabel('|a_{COM}| (m/s^2)');
legend({'S1(AB)','S2(BC)','S3(DCE)','S4(EF)','S5(FGH)'},'Location','best');

% Total computation and graph plotting time
elapsedTime = toc;
fprintf('Total computation and graph plotting time: %.2f seconds\n', elapsedTime);