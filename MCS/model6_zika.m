function dy = model6_zika(y,k)

global Lambda_v mu_v Lambda mu omega gamma gammap xi

dy = zeros(9,1);



beta = k(1);
beta_v = k(2);
beta_vp = k(3);



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
dy(6) = beta_v* Iv*S/(S+I+R+Sp+Ip) + beta_vp * Iv * Sp/(S+I+R+Sp+Ip);
dy(9) = beta_vp * Iv * Sp/(S+I+R+Sp+Ip);


end