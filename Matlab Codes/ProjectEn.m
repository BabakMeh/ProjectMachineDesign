%% Project Machine Design 2
% Babak Mehdizadeh & Nima mohammadi
% 995241050

clc;
clear;
%% Data and Inputs
% oil input temperature
Ti=60;
% radial force
W=3000;
% Estimated Delta T
delta_T=input('estimated delta_T(delta_T) at centigrade: ');
% rotational speed
N=8;
% pressure
pressure=0.469;
P=10^6*pressure;
% diameter
d=80;
%  radial lag
c=input('radial lag(c) at mm :');

%%
r=d/2;
% average tempreture in centigrade
Tf1=Ti+(delta_T)/2;
% TF1 at Fahrenheit
TF1=32+(Tf1*1.8);
disp(['Tf: ' , num2str(TF1) , 'Farenheit' ]);
% absolute viscosity
% with using diagram 12.14
mu1=input([' with using diagram 12.14 Enter the value of absolute viscosity at' num2str(TF1) ,'Fahrenheit: ' ]);
disp(['mu(viscosity): ' , num2str(mu1) , 'mureyn']);
mu_1=mu1*6.89*0.001; 
% Sommerfeld
S=((r/c)^2)*mu_1*N/P;
disp(['S(Sommerfeld): ' , num2str(S)]);


%find Temperature Rise Dimensionless Variable
% l/d=1

u=0.349109+S*6.0094+S*S*0.047467;
Delta_C=(10^-6)*P*u/0.12;

Delta_T_new=(delta_T+Delta_C)/2;
disp(['Delta_T_new: ' , num2str(Delta_T_new) ,'centigrade']);

% check mureyn to mpa.s convert
% check farenheit to centigrade convert

%% 
% now we should check the Delta_T_new that we fine is
% correct or not 
e=input('If Delta_T_new is near to Delta_T Enter 1 else Enter 0 :');

while e==0
    % TF2 at Fahrenheit
   Tf2=Ti+(Delta_T_new)/2; 
   TF2=32+(Tf2*9/5);
   % with using diagram 12.14
   mu2=input([' with using diagram 12.14 Enter the value of absolute viscosity at' num2str(TF2) , 'Fahrenheit: ']);
   mu_2=mu2*6.89*0.001; 
   S=((r/c)^2)*mu_2*N/P;
   disp(['S(Sommerfeld): ' , num2str(S)]);
   u=0.349109+S*6.0094+S*S*0.047467;
   Delta_C=(10^-6)*P*u/0.12;
   disp(['Delta_C :', num2str(Delta_C)]);
   Delta_T_new=(Delta_T_new+Delta_C)/2;
   disp(['Delta_T_new: ' , num2str(Delta_T_new)]);
   e=input('If Delta_T_new is near to Delta_T Enter 1 else Enter 0 :');
   
end
 disp(['final value of Delta_T :' , num2str(Delta_T_new)]);
%% 
% now after that find definitive value of Delta_T we have:

Tf_final=Ti+(Delta_T_new/2);
TF_Final=32+(Tf_final*9/5);
To=Ti+Delta_T_new;
  
mu_final=input(['Enter the value of absolute viscosity at' num2str(TF_Final) , 'Fahrenheit: ']);
mu_f=mu_final*6.89*0.001; 
S=((r/c)^2)*mu_f*N/P;

disp(['S(definitive value of Sommerfeld): ' , num2str(S)]);


% with using diagram 12.16 we read the value of h0/c 
% h0/c=h1
h1=input(' with using diagram 12.16  Enter the value of h0/c: ');
h0=h1*c;

% with using diagram 12.18 we read the value of rf/c
%rf/c=r1
r1=input('with using diagram 12.18 Enter the value of rf/c: ');
f=r1*c/r;

% with using diagram 12.19 we read the value of Q/rcNl
% Q1=Q/rcNl
Q1=input('with using diagram 12.19 Enter the value of Q/rcNl: ');
Q=Q1*r*c*N*d;

% with using diagram 12.20 we read the value of Qs/Q
% Q2=Qs/Q
Q2=input('with using diagram 12.20 Enter the value of Qs/Q: ');
Qs=Q*Q2;
% with using diagram 12.21 we read the value of P/Pmax
% P1= P/Pmax
P1=input('with using diagram 12.21  Enter the value of P/Pmax: ');
P_max=(P/P1)*10^(-6);

%H_loss
H_loss=f*W*d*pi*N/1000;


%% 
disp(['S(definitive value of Sommerfeld): ' , num2str(S)]);
disp(['Tf = ',num2str(Tf_final) , 'centigrade']);
disp(['To = ',num2str(To) , 'centigrade']);
disp(['h0 = ',num2str(h0) , 'm']);
disp(['f = ',num2str(f)]);
disp(['H_loss = ',num2str(H_loss) ,'Watt']);
disp(['Q = ',num2str(Q) , ' mm^3/s']);
disp(['Qs = ',num2str(Qs) , 'mm^3/s']);
disp(['pmax = ',num2str(P_max),' (Mpa)']);