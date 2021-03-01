function [s_path, vol_path, V0] = Heston(S0, vol0, T, K, Theta, v, r, rho, NSteps, MPaths)

mu = [0,0];
sigma = [1,0;0,1];
dt = T/NSteps;

s_path = zeros(NSteps+1, MPaths);
s_path(1,:) = S0;     % First row(step) is the initial price
vol_path = zeros(NSteps+1, MPaths);
vol_path(1,:) = vol0;  % First row(step) is the initial volatility

k0 = -rho*K*Theta/v * dt;
k1 = 0.5*dt*(K*rho/v - 0.5) - rho/v;
k2 = 0.5*dt*(K*rho/v - 0.5) + rho/v;
k3 = 0.5*dt*(1-rho^2);
k4 = 0.5*dt*(1-rho^2);

for i = 1:MPaths
    for j = 1:NSteps
        pair_nm = mvnrnd(mu,sigma,1);   % generate a pair of two-dimensional norm
        vol_plus = max(vol_path(j,i),0);    % Full truncation scheme given by Lord, R., R. Koekkoek and D. van Dijk (2006)
        vol_path(j+1,i) = vol_path(j,i) + K*(Theta - vol_plus)*dt + v*(vol_plus^.5)*pair_nm(1,2)*(dt^.5);
        % vol_path(j+1,i) = max(vol_path(j,i),0);
        log_s = log(s_path(j,i)) + k0 + k1*vol_path(j,i) + k2*vol_path(j+1,i) + sqrt( k3*max(vol_path(j,i),0) + k4*max(vol_path(j+1,i),0) )*pair_nm(1,1);
        s_path(j+1,i) = exp(log_s);    % taking log ensure stock price not going negative

    end
end

% Calculate the price of the option
end_values = max(s_path(NSteps+1,:)-mean(s_path), 0);
V0_track = zeros(1, MPaths);
for i = 1:MPaths
    V0_track(1,i) = exp(-r*T) * mean(end_values(1,1:i));   % Calculate a mean after each new path, so that we keep trakc of estimation
end

V0 = V0_track(1,MPaths);
%V0 = exp(-r*T) * mean(end_values);

figure(1)
subplot(2,1,1);
plot(0:T/NSteps:T,s_path);
title('sample paths of Stock Price','fontsize',14);
ylabel('price','fontsize',14);
set(gca,'fontsize',14,'FontWeight','bold');
xlabel('time','fontsize',14);

subplot(2,1,2);
plot(0:T/NSteps:T,vol_path);
title('sample paths of Stochastic Volatility','fontsize',14);
ylabel('volatility','fontsize',14);
set(gca,'fontsize',14,'FontWeight','bold');
xlabel('time','fontsize',14);

figure(2)
plot(1:MPaths,V0_track);
title('Track of Option Price Estimation','fontsize',14);
ylabel('price','fontsize',14);
set(gca,'fontsize',14,'FontWeight','bold');
xlabel('mth_path','fontsize',14);



