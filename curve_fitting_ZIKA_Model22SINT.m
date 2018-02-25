function curve_fitting_ZIKA_Model22SINT

clear all
close all
clc


% load ZikaFlorida.txt
% tdata = ZikaFlorida(:,1);
% ydata = ZikaFlorida(:,2);

% load ZikaLocal.txt
% tdata = ZikaLocal(:,1);
% ydata = ZikaLocal(:,2);

load ZikaLocal2.txt
tdata = ZikaLocal2(:,1);
ydata = ZikaLocal2(:,2);

tforward = [tdata(1):0.1:tdata(end)+2]';

format long

Lambda_v = 1/10;
mu_v = 1/10;

mu = 1/(79*365);
Lambda = 19380000*mu;
% Lambda = 20000000*mu;

omega =0.0;
phi=0.2;
q=1;



% 
% k = 1.0e+03 * [0.000000028227958   0.000000005328385   
%                0.000026657369222   8.327830807908896] %res 310, best fit so far
           
% k = 1.0e+03 * [0.000000028227958     0.000000005328385   
%                0.00026657369222    0.00008327830807908896]
% 
% k = [0.000015472837978   0.000005456979318...
%      0.14                0.0013298604935864];


%k = [ 0.000084022423772   0.000006181789114   0.060000008322600   0.010000317188820  0.000000000003361]; %res3.410719528197382e+02
% initial_cond = [10^6 10^2 19380000 3 0 1];
%k = [ 0.000084022423772   0.000006181789114   0.060000008322600     0.000000000003361];

%k=[   0.000145738206859   0.000005178955972   0.100000329960743   0.000000000002683  0.0000001]; % 3.718343395113619e+02 lb gamma=1/10
% k =[0.64071787447   0.34170627996   0.100004932397014   0.105355950394680 ];
%k = 1.0e+03 *[2.792794440474008   0.493028184910855   0.000100113410368   0.00015610690177802]; % SI
% k = 1.0e+02 *[9.345581686266019   4.440596271170688 ...
%               0.001000024352960   0.001805296476230]
k = 1.0e+02 *[7.013391469297598   4.817397302890158 ...
              0.001000212579542   0.001310810966903];
%k=[0.000438514438704   0.000018825657362   0.998635777840970   3.981738121266249]% res 2.071936876458121e+02

initial_cond = [0.99 0.01 19380000 3 0 1  10];




function dy = model_zikav1(y,k)

dy = zeros(7,1);



beta = k(1);
beta_v = k(2);
gamma = k(3);
gammaA = k(4);




Sv = y(1);
Iv = y(2);
S = y(3);
I = y(4);
R = y(5);
C = y(6);

A = y(7);

dy(1) = Lambda_v - beta* (I+q*A)* Sv/(S+I+A+R) - mu_v*Sv;
dy(2) = beta*(I+q*A)*Sv/(S+I+A+R) - mu_v*Iv;
dy(3) = Lambda - beta_v* Iv*S/(S+I+A+R)  - mu*S + omega*R;
dy(7) = (1-phi)*beta_v*Iv*S/(S+I+A+R)  - (mu + gammaA)*A;
dy(4) = phi* beta_v*Iv*S/(S+I+A+R)  - (mu + gamma)*I;
dy(5) = gamma*I + gammaA * A - (mu + omega)*R;
dy(6) = phi* beta_v* Iv*S/(S+I+A+R);
end

function error_in_data = err_in_data(k)


 [t,y] = ode45(@(t,y)model_zikav1(y,k),tdata,initial_cond);
 
 CI = y(:,6); %cumulative incidences
  
 %Incidences = k(2)*y(:,2)*y(:,3); 
 error_in_data = sum((CI - ydata).^2) ;
 
 

end


lb = [0      0      1/10  1/10  ]; 
ub = [10000  10000  0.5   0.5 ];

 for j=1:5 
[k,fval] =  fminsearchbnd(@err_in_data,k,lb,ub,optimset('Display','iter','MaxIter',200,'MaxFunEvals',200)) 
% k = fmincon(@err_in_data,k,[],[],[],[],lb,ub,[],...
%     optimoptions('fmincon','Display','iter'))
 end


r1 = (k(1)*k(2)*Lambda_v*mu)
r2 = (mu_v^2  * (mu + k(3)) *Lambda)
r3 = (mu_v^2 * (mu + k(4))*Lambda)
Reproduction_Number = phi * r1/r2 + (1-phi)*q* r1/r3

 figure(1)

 [t y] = ode45(@(t,y)model_zikav1(y,k),tforward,initial_cond);
 plot(tdata, ydata, 'o', tforward, y(:,6), '-r')
 
 figure(2)

 tprojection = 0:0.1:2000;
%  [t y] = ode45(@(t,y)model_zikav1(y,k),tprojection,initial_cond);
%  plot(tdata, ydata, 'o', tprojection, y(:,6), '-r')
 
 

   [t y] = ode45(@(t,y)model_zikav1(y,k),tdata,initial_cond);
residuals = (ydata - y(:,6))./y(:,6);
mn = mean(residuals)*ones(1,length(ydata));
  plot(tdata, residuals, 'or', tdata, mn,'-b')
  
 figure(3)

 [t y] = ode45(@(t,y)model_zikav1(y,k),tprojection,initial_cond);
 plot(tprojection, y(:,4), '-r')
% % 
 end