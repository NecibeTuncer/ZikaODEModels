function  total_ARE = MonteCarlo_ZikaEpiModel6_R1

clear all
close all
clc

numiter = 1000; 


X = zeros(5,numiter);

Lambda_v = 1/10;
mu_v = 1/10;

mu = 1/(79*365);
Lambda = 19380000*mu;
omega = 0.0;
xi = 0.0000007;


tdata = [1 3 9 14 15 18 21 22 23 24 25 30 31 32 36 37 39 42 43 44 50 56 57 58 59 60 63 64 66 67 71 72 73];

tforward = (1:0.1:73);

tmeasure = [1 21 81 131 141 171 201 211 221 231 241 291  301 311 351 361 381 411 421 431 491 551 561 571 581 591 ...
            621 631 651 661 701 711 721];

tpmeasure = [21 291 311 421 431 631];

true_params = 1.0e+03 *[2.617723435993325   0.091377380457317   0.000100000000001   0.399007853530409...
              0.000020000000003];
   
initial_cond = [0.99 0.01 19380000 3 0 1 387600 1 1];

[~,y_trp] = ode45(@(t,y)(model_zikav1(y,true_params)),tforward,initial_cond);





function dy = model_zikav1(y,k)

dy = zeros(8,1);



beta = k(1);
beta_v = k(2);
gamma = k(3);
beta_vp = k(4);
gammap = k(5);


Sv = y(1);
Iv = y(2);
S = y(3);
I = y(4);
R = y(5);
C = y(6);

Sp = y(7);
Ip = y(8);
Cp = y(9);

dy(1) = Lambda_v - beta* (I+ Ip)* Sv /(S+I+R+Sp +Ip) - mu_v*Sv;
dy(2) = beta*(I+Ip)*Sv/(S+I+R+Sp +Ip) - mu_v*Iv;
dy(3) = Lambda - beta_v* Iv*S/(S+I+R+Sp +Ip) - (mu+xi)*S + omega*R + (1/280)*Sp;
dy(7) = xi*S - beta_vp * Iv * Sp/(S+I+R+Sp +Ip)  - mu*Sp - (1/280)*Sp;
dy(4) = beta_v*Iv*S/(S+I+R+Sp +Ip)  - (mu + gamma)*I;
dy(8) = beta_vp * Iv * Sp/(S+I+R+Sp +Ip)  -(mu+gammap) *Ip;
dy(5) = gamma*I +gammap * Ip - (mu + omega)*R;
dy(6) = beta_v* Iv*S/(S+I+R+Sp+Ip) +beta_vp * Iv * Sp/(S+I+R+Sp+Ip) ;
dy(9) = beta_vp * Iv * Sp/(S+I+R+Sp+Ip);
end

function error_in_data = err_in_data(k)

 [~,y] = ode45(@(t,y)model_zikav1(y,k),tdata,initial_cond);
 [~,yp] = ode45(@(t,y)model_zikav1(y,k),tforward,initial_cond);
 CI = y(:,6)'; %cumulative incidences
 CIP = yp(tpmeasure(:),9)'; 

 error_in_data = sum((CI - ydata).^2) + sum((CIP - ypdata).^2) ;
 
 

end


noiselevel = [0, 0.01, 0.05, 0.1, 0.2, 0.3];
total_ARE =  zeros(length(noiselevel), length(true_params));

for noisei = 1:6
    
rng default
noiselev = noiselevel(noisei)


for i= 1:numiter
    i
  ytrp = y_trp(tmeasure(:),6);
  ytrpreg = y_trp(tpmeasure(:),9);
  
  ydata =  (noiselev*(ytrp').*randn(1,size(tmeasure,2))) + ytrp';
  ypdata =  (noiselev*(ytrpreg').*randn(1,size(tpmeasure,2))) + ytrpreg';

k = true_params;

lb = [0        0       1/10     0            1/50 ];
%ub = [10000  10000   0.5  0.5  10000];
 
k =  fminsearchbnd(@err_in_data,k,lb,[],optimset('MaxIter',2000,'MaxFunEvals',2000)) 
                           
X(:,i) = k';
%conv(noisei,i) = exitflag;
 

% plot(tdata,ydata,'o',tforward,y_trp(:,6),'LineWidth',1.5)
% title('Cumulative New Cases Generated Data')

end



arescore = zeros(1,length(true_params));
    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    total_ARE(noisei,:) = arescore;



end

save('Model6MonteCarlo','total_ARE')
end