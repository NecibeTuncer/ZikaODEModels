function curve_fitting_ZIKA_Model24SINT

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



% k=[0.000090075113591   0.000019152204649   0.100000000000000   0.379079535409867 ...
%                  0   0.140867306273061]; res 393.687
 % k=[ 0.000088711852099   0.000019374781996   0.100000000207413   0.378322961661689 ...
 %                  0   0.140867306273061]; %res 389.7
               
% k =[0.000087542690288   0.000020013604417   0.100043643298911   0.398746186199319 ...
%                  0   0.140867306273061]; %res 389.6
%   k= [87.411816098   200.80500712   0.100000236351274   0.400936005694349 ...
%                                 0   0.140867306273061];
% 
%   k = 1.0e+03 *[1.728482807660942   0.283077269452301 ...
%                 0.000100011477283   0.000499180621172 ...
%                 0.000122100174875   0.000683221062881];                          
                            %k=[0.000438514438704   0.000018825657362   0.998635777840970   3.981738121266249]% res 2.071936876458121e+02
% k =[64.801841203403484  58.117625279474062...
% 0.100001086084792   0.496358810242125...
% 0.340635923078425   0.805718653133023]

% k =[43.787049415002222  53.770680027218631...
%     0.100000013188961   0.498934984148567...
%      0.335406739235267   0.843767953339670];

k =1.0e+03*[2.825038695844068   0.100390800239187...
            0.000100000000446   0.000100000000446... 
            0.000000000938470   1];
         
initial_cond = [0.99 0.01 19380000 3 0 1  12];




function dy = model_zikav1(y,k)

dy = zeros(7,1);



beta = k(1);
beta_v = k(2);
gamma = k(3);
gammaA = k(4);
beta_d = k(5);
q_A = k(6);




Sv = y(1);
Iv = y(2);
S = y(3);
I = y(4);
R = y(5);
C = y(6);

A = y(7);

dy(1) = Lambda_v - beta* (I+q*A)* Sv/(S+I+R+A) - mu_v*Sv;
dy(2) = beta*(I+q*A)*Sv/(S+I+R+A) - mu_v*Iv;
dy(3) = Lambda - (beta_v* Iv*S + beta_d * S *(I+q_A * A))/(S+I+R+A) - mu*S + omega*R;
dy(7) = (1-phi)*( beta_v*Iv*S + beta_d * S *(I+q_A * A))/(S+I+R+A) - (mu + gammaA)*A;
dy(4) = phi* (beta_v*Iv*S + beta_d * S *(I+q_A * A))/(S+I+R+A) - (mu + gamma)*I;
dy(5) = gamma*I + gammaA * A - (mu + omega)*R;
dy(6) = phi*( beta_v* Iv*S + beta_d * S *(I+q_A * A))/(S+I+R+A);
end

function error_in_data = err_in_data(k)


 [t,y] = ode45(@(t,y)model_zikav1(y,k),tdata,initial_cond);
 
 CI = y(:,6); %cumulative incidences
  
 %Incidences = k(2)*y(:,2)*y(:,3); 
 error_in_data = sum((CI - ydata).^2) ;
 
 

end


lb = [0      0      0.1   0.1   0      1]; 
ub = [10000  10000  0.5   0.5   10000  1];

 for j=1:5
[k,fval] =  fminsearchbnd(@err_in_data,k,lb,ub,optimset('Display','iter','MaxIter',1000,'MaxFunEvals',1000)) 
% k = fmincon(@err_in_data,k,[],[],[],[],lb,ub,[],...
%   optimoptions('fmincon','Display','iter'))
 end


r1 = (k(1)*k(2)*mu*Lambda_v)
r2 = (mu_v^2 * Lambda * (mu + k(3)))
r3 = (mu_v^2 * Lambda * (mu + k(4)))
Reproduction_Number = phi * r1/r2 + (1-phi)*q* r1/r3 + phi*k(5)/((mu+k(3))) ...
                           +(1-phi)*k(6)*k(5)/((mu+k(4)))

 figure(1)

 [t y] = ode45(@(t,y)model_zikav1(y,k),tforward,initial_cond);
 plot(tdata, ydata, 'o', tforward, y(:,6), '-r')
 
%  figure(2)
% 
%  tprojection = 0:0.1:2000;
% %  [t y] = ode45(@(t,y)model_zikav1(y,k),tprojection,initial_cond);
% %  plot(tdata, ydata, 'o', tprojection, y(:,6), '-r')
%  
%  
% 
%    [t y] = ode45(@(t,y)model_zikav1(y,k),tdata,initial_cond);
% residuals = (ydata - y(:,6))./y(:,6);
% mn = mean(residuals)*ones(1,length(ydata));
%   plot(tdata, residuals, 'or', tdata, mn,'-b')
%   
%  figure(3)
% 
%  [t y] = ode45(@(t,y)model_zikav1(y,k),tprojection,initial_cond);
%  plot(tprojection, y(:,4), '-r')
% % % 
 end