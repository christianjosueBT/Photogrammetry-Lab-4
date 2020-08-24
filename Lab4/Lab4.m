clear all;
close all;
clc;
format longg;

% loading image coordinates and residuals 
image_3314 = load("Coords_3314.txt");
image_3316 = load("Coords_3316.txt");
vhat_3314 = load("Residuals_3314.txt");
vhat_3316 = load("Residuals_3316.txt");

%% 3314 image

% residuals
vhat_plot = zeros(14,2);
for i=1:length(vhat_plot)
	vhat_plot(i,1) = vhat_3314(i*2 - 1,1);
	vhat_plot(i,2) = vhat_3314(i*2,1);
end
% rmse of residuals
rmseX_3314 = rms(vhat_plot(:,1))
rmseY_3314 = rms(vhat_plot(:,2))
rmse_3314 = rms(vhat_3314(:))
% plotting residuals
figure;
quiver(image_3314(:,2), image_3314(:,3),vhat_plot(:,1),vhat_plot(:,2));
title({'3314 Image Residual Plot'});
xlabel('X (m)');
ylabel('Y (m)');
grid on;

%% 3316 image

% residuals
vhat_plot = zeros(14,2);
for i=1:length(vhat_plot)
	vhat_plot(i,1) = vhat_3316(i*2 - 1,1);
	vhat_plot(i,2) = vhat_3316(i*2,1);
end
% rmse of residuals
rmseX_3316 = rms(vhat_plot(:,1))
rmseY_3316 = rms(vhat_plot(:,2))
rmse3316 = rms(vhat_3316(:))
% plotting residuals 
figure;
quiver(image_3316(:,2), image_3316(:,3),vhat_plot(:,1),vhat_plot(:,2));
title({'3316 Image Residual Plot'});
xlabel('X (m)');
ylabel('Y (m)');
grid on;