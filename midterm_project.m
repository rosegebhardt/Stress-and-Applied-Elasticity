%% Rose Gebhardt & Hannah Quirk -- December 16, 2019 -- Stress and Applied Elasticity Mini-Project
clear all; close all; clc;

%% Equations Used

% Stress in the radial direction (Equation 3.55a)
sigma_r = @(ratio,theta) 0.5*((1-ratio.^-2)+((1+3*(ratio.^-4)-4*(ratio.^-2)).*cos(2*theta)));

% Stress in theta direction (Equation 3.55b)
sigma_t = @(ratio,theta) 0.5*((1+ratio.^-2)-((1+3*(ratio.^-4)).*cos(2*theta)));

% Shear in the radial-theta direction (Equation 3.55c)
tau_rt = @(ratio,theta) -0.5*(1-3*(ratio.^-4)+2*(ratio.^-2)).*sin(2*theta);

% von-Mises yield criterion
sigma_e = @(sigma_r,sigma_t,tau_rt) (sigma_r.^2 - sigma_r.*sigma_t + sigma_t.^2 + 3*(tau_rt.^2)).^0.5;

%% Task 1 -- Plots 

% Ratio of radius to hole radius (r/a)
domain_1 = linspace(1,4,1000);

% Plot yield criterion again radius ratio for theta = 0,30,60,90 degrees
figure(1)
plot(domain_1,sigma_e(sigma_r(domain_1,0), sigma_t(domain_1,0), tau_rt(domain_1,0)),...
    domain_1,sigma_e(sigma_r(domain_1,pi/6), sigma_t(domain_1,pi/6), tau_rt(domain_1,pi/6)),...
    domain_1,sigma_e(sigma_r(domain_1,pi/3), sigma_t(domain_1,pi/3), tau_rt(domain_1,pi/3)),...
    domain_1,sigma_e(sigma_r(domain_1,pi/2), sigma_t(domain_1,pi/2), tau_rt(domain_1,pi/2)))
xlim([0,4])
title('Yield Criterion as a Function of Radius')
xlabel('Ratio of Radius to Hole Radius'); ylabel('Ratio to Yield Criterion to Applied Distributed Load')
legend('0 degrees', '30 degrees', '60 degrees', '90 degrees')

%% Task 2 -- Stress Over Region

% Define domain
x_domain = linspace(-4,4,1000); y_domain = linspace(-4,4,1000);
[x,y] = meshgrid(x_domain,y_domain);
XY = [x(:),y(:)];

% Excluse data inside hole radius
for i = 1:1000000
    if XY(i,1)^2 + XY(i,2)^2 < 1
        XY(i,:) = [inf,inf];
    end
end

% Convert from Cartesian to radial coordinates
R = (XY(:,1).^2 + XY(:,2).^2).^0.5;
TH = atan(XY(:,2)./XY(:,1));

% Find yield condition at each point on domain
S_0 = abs(sigma_e(sigma_r(R,TH),sigma_t(R,TH),tau_rt(R,TH)));
S = reshape(S_0,length(x),length(y));

% Plot surface
figure(2)
surf(x,y,S,'EdgeColor','interp')
title('Yield Condition on an Infinite Plate with Circular Hole')
xlabel('$$\frac{x}{a}$$','interpreter','latex'); 
ylabel('$$\frac{y}{a}$$','interpreter','latex'); 
zlabel('$$\frac{\sigma_e}{s}$$','interpreter','latex');
colorbar

% Create contour plot
figure(3)
contour(x,y,S,20)
title('Contour Plot of Yield Condition')
xlabel('$$\frac{x}{a}$$','interpreter','latex'); 
ylabel('$$\frac{y}{a}$$','interpreter','latex'); 
zlabel('$$\frac{\sigma_e}{s}$$','interpreter','latex');
colorbar

%% Task 3 -- Maximum Yield Stress

% Get maximum yield criterion over domain and index number in array  
[M,I] = max(S_0);

% Use index to get xy-coordinate of maximum stress
max_pos = XY(I,:);
