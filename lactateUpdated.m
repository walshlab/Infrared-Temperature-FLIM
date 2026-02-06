%% Import Data 
%Importing data given from model results xlsx file 
temp1= importdata('temp1.txt'); %Temp 1 
temp2= importdata('temp2.txt'); %Temp 2
msec1= importdata('time1ms.txt'); % Time 1 
msec2= importdata('time2ms.txt'); %Time 2
%% Plot only 

% Exponential Decay Fitting
% Define the exponential decay model
expModel = fittype('A*exp(-x/tau) + C', 'independent', 'x', 'coefficients', {'A', 'tau', 'C'});

% Fit for temp1
[fitResult1, msec1_sub, temp1_sub] = fitExponentialDecay(msec1, temp1, expModel);

% Fit for temp2
[fitResult2, msec2_sub, temp2_sub] = fitExponentialDecay(msec2, temp2, expModel);


% Display fitted parameters
disp('Fit Parameters for Temp1:');
disp(fitResult1);
disp('Fit Parameters for Temp2:');
disp(fitResult2);

% Plot Results in a Single Figure
figure;
hold on;
% Raw Data
plot(msec2, temp2,'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6); % Temp2 data
plot(msec1, temp1,'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6); % Temp1 data

% Exponential Fits (Smooth Lines)
plot(msec2_sub, fitResult2(msec2_sub), 'k-', 'LineWidth', 3); % Fit for Temp2
plot(msec1_sub, fitResult1(msec1_sub), 'k:.', 'LineWidth', 3); % Fit for Temp1



 %Labels and Legends
xlabel('Time (ms)', 'FontSize', 16);
ylabel('Temperature (°C)', 'FontSize', 16);

% Define Equations for the Fits (LaTeX format)
eq1 = sprintf('$y = %.2f e^{-x/%.2f} + %.2f$', fitResult1.A, fitResult1.tau, fitResult1.C);
eq2 = sprintf('$y = %.2f e^{-x/%.2f} + %.2f$', fitResult2.A, fitResult2.tau, fitResult2.C);

legend({'0.96 J/cm$^2$ (5.01 ms)', ...
        '0.49 J/cm$^2$ (2.73 ms)', ...
        eq2, eq1}, ...
       'Interpreter', 'latex', 'Location', 'Northeast', 'FontSize', 16);

% Aesthetic Adjustments
set(gca, 'LineWidth', 1);           % Axis line thickness
set(gca, 'TickLength', [0.01, 0.01]);  % Tick mark length
box on 

hold off;


% Exponential Decay Fitting for Thermal Gradients 
function fitResult = fitSingleExponentialDecay(msec, temp, expModel, color)
    % Find index of maximum temperature
    [maxTemp, maxIdx] = max(temp);  

    % Select only data from the maximum temperature onward
    msec_sub = msec(maxIdx:end);
    temp_sub = temp(maxIdx:end);

    % Set initial parameter guesses
    Tmin = min(temp_sub); % Baseline (final temperature)
    A_init = maxTemp - Tmin; % Initial amplitude
    tau_init = (max(msec_sub) - min(msec_sub)) / 3; % Rough estimate for decay rate

    startPoints = [A_init, tau_init, Tmin];

    % Fit the model to the subset data
    fitResult = fit(msec_sub, temp_sub, expModel, 'Start', startPoints);

    % Plot results
    hold on;
    plot(msec_sub, fitResult(msec_sub), color, 'LineWidth', 2); % Fitted curve
    hold off;
end
%% Establishing Frequency and activation energy for lactate dehydrogenase 
act= 10393; %activation energy  
freq= 4.2/exp(-act/(8.314*(293))) %Finding frequency factor 

%% 1st Arrhenius Theory k=Ae^(-Ea/RT) for lactate dehydrogenase
%Temperature indexes 0-29 hold the linear region (exposure) 
%Temperature indexes 30-328 hold the exponential region
temp1k= temp1+273; %Conversion from celsius to kelvin 
RT1= 8.314*(temp1k); %Tempature multiplied by R Constant
ep1= -(act.*ones(length(temp1k),1)); %Activation energy 
div1= rdivide(ep1, RT1); %Division of Ea/RT
k1= freq*exp(div1); %Calculation of R constant 

figure(2) %Figure for Optica 
plot3(msec1,temp1k, k1)
grid on
fontSize= 16;
xlabel('Time (ms)', 'fontSize', fontSize)
ylabel('Temp(K)', 'fontSize', fontSize)
zlabel('Rate Constant, k', 'fontSize', fontSize)
title('LDH Rate Constant, k- 0.49 J/cm^2 (2.73 ms)', 'fontSize', fontSize)
legend('0.49 J/cm^2 (2.73 ms)', '0.96 J/cm^2 (5.01 ms)','fontSize', fontSize)
set(gca, 'LineWidth', 1);           % Axis line thickness
set(gca, 'TickLength', [0.01, 0.01]);  % Tick mark length

hold off

%% 2nd Arrhenius Theory k=Ae^(-Ea/RT) for lactate dehydrogenase
%Temperature indexes 0-29 hold the linear region (exposure) 
%Temperature indexes 30-328 hold the exponential region
temp2k= temp2+273; %Conversion from celsius to kelvin 
RT2= 8.314*(temp2k); %Tempature multiplied by R Constant
ep2= -(act.*ones(length(temp2k),1)); %Activation energy 
div2= rdivide(ep2, RT2); %Division of Ea/RT
k2= freq*exp(div2); %Calculation of R constant 

figure(3)
plot3(msec2,temp2k, k2)
xlabel('Time (ms)', 'fontSize', fontSize)
ylabel('Temp(K)', 'fontSize', fontSize)
zlabel('Rate Constant, k', 'fontSize', fontSize)
title('LDH Rate Constant, k- 0.96 J/cm^2 (5.01 ms)', 'fontSize', fontSize)
% lol= annotation('textbox', [0.45, 0.6, 0.1, 0.1], 'String', "Freq= 299.32")
% lol.FontSize= 13;
% lol2= annotation('textbox', [0.45, 0.5, 0.1, 0.1], 'String', "Ea= 10,393")
% lol2.FontSize=13; 


% Exponential Decay Fitting
% Define the exponential decay model
expModel = fittype('A*exp(-x/tau) + C', 'independent', 'x', 'coefficients', {'A', 'tau', 'C'});

% Fit for k1
[fitResult1, msec1_sub, k1_sub] = fitExponentialDecay(msec1, k1, expModel);

% Fit for k2
[fitResult2, msec2_sub, k2_sub] = fitExponentialDecay(msec2, k2, expModel);


% Display fitted parameters
disp('Fit Parameters for Temp1:');
disp(fitResult1);
disp('Fit Parameters for Temp2:');
disp(fitResult2);

% Plot Results in a Single Figure
figure;
hold on;
% Raw Data
plot(msec1, k1,'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6); % Temp1 data
plot(msec2, k2,'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6); % Temp2 data

% Exponential Fits (Smooth Lines)
plot(msec1_sub, fitResult1(msec1_sub), 'k--', 'LineWidth', 2); % Fit for Temp1
plot(msec2_sub, fitResult2(msec2_sub), 'k-', 'LineWidth', 2); % Fit for Temp2

 %Labels and Legends
xlabel('Time (ms)', 'FontSize', 16);
ylabel('Rate  Constant (k)', 'FontSize',16);
%ylabel('Temperature (°C)', 'FontSize', 16);
% ---- Add Second X-Axis ----
ax1 = gca; % Get current axis
ax2 = axes('Position', ax1.Position, 'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', 'Color', 'none'); % Create secondary axis

% ---- Convert Time (ms) to Temperature (°C) ----
% Assuming a known conversion between time and temperature (example below)
temp2_xaxis = temp2; % Get temperature values from the fit

% Set the secondary x-axis ticks and labels
set(ax2, 'XColor', 'k', 'YColor', 'none'); % Hide secondary y-axis
ax2.XLabel.String = 'Temperature (°C)';
ax2.XLabel.FontSize = 16;

% Match the second x-axis ticks with the first x-axis but with temperature values
set(ax2, 'XTick', get(ax1, 'XTick')); 
set(ax2, 'XTickLabel', num2str(temp2_xaxis(get(ax1, 'XTick'))', '%.1f'));

% ---- Add Legends ----
eq1 = sprintf('$y = %.2f e^{-x/%.2f} + %.2f$', fitResult1.A, fitResult1.tau, fitResult1.C);
eq2 = sprintf('$y = %.2f e^{-x/%.2f} + %.2f$', fitResult2.A, fitResult2.tau, fitResult2.C);

legend({'0.49 J/cm$^2$ (2.73 ms)', ...
        '0.96 J/cm$^2$ (5.01 ms)', ...
        eq1, eq2}, ...
       'Interpreter', 'latex', 'Location', 'Northeast', 'FontSize', 16);

% ---- Adjust Aesthetics ----
set(gca, 'LineWidth', 1);  
set(gca, 'TickLength', [0.01, 0.01]);  
box on;

hold off;

%% 3D Rate Constant Plots 

plot3(msec2,temp2k,k2, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6,'DisplayName', '0.96 J/cm^2 (5.01 ms)' ,'LineWidth', 1.5);
hold on 
plot3(msec1,temp1k,k1,'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'DisplayName', '0.49 J/cm^2 (2.73 ms)', 'LineWidth', 1.5);

fontSize= 16;
%title('Rate Constants')
xlabel('Time (ms)', 'fontSize', fontSize)
ylabel('Temp (K)', 'fontSize', fontSize)
zlabel('Rate Constant (k)', 'fontSize', fontSize)
legend( 'Location', 'northeast')
set(gca, 'LineWidth', 1);           % Axis line thickness
set(gca, 'TickLength', [0.01, 0.01]);  % Tick mark length
grid off 
box on 

hold off



%% Setting up parameters for relationship 
enzyme= 14.079; %Enzyme concentration in um 
nadh= 50; %NAD(P)H concentration in um 
s= 4; %Integer number of binding sites per enzyme molecule 

%% Alpha1 vs k rate constant data 1 
%Quadratic equation set up 
k= k1;
a= nadh;
b= (nadh*s*enzyme)-(nadh.^2)+k;
c= -k.*(nadh);

a1_k1= (-b+ sqrt((b.^2)-4*a*c))/(2*a);


figure(6)
plot(k, a1_k1)
xlabel('k rate constant')
ylabel('alpha 1')
title('k vs alpha 1 lactate dehydrogenase Data 1')

figure(7)
line(msec1, a1_k1, 'color', 'r')
xlabel('Time (ms)') 
ylabel('alpha 1')
ax1= gca; %current axes 
ax1_pos= ax1.Position; %Position of first axes 
ax2= axes('Position', ax1_pos, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');


%% Alpha1 vs k rate constant data 2
%Quadratic equation set up 
k3= k2; 
a2= nadh;
b2= (nadh*(s)*(enzyme))-(nadh.^2)+k3;
c2= -k3.*(nadh);

a1_k2= (-b2+ sqrt((b2.^2)-4*a2*c2))/(2*a2);


figure(8)
plot(k3, a1_k2)
xlabel('k rate constant')
ylabel('alpha 1')
title('k vs alpha 1 lactate dehydrogenase Data 2')



figure(9)
line(msec2, a1_k2, 'color', 'r')
xlabel('Time (ms)') 
ylabel('alpha 1')
ax1= gca; %current axes 
ax1_pos= ax1.Position; %Position of first axes 
ax2= axes('Position', ax1_pos, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');
%line(temp2k,a1_k2, 'Parent',ax2, 'color', 'k')
%xlabel('Temp (K)')
%ylabel('alpha 1')
%annotation('textbox', [0.25, 0.6, 0.1, 0.1], 'String', "lactate dehydrogenase Data 2")
%annotation('textbox', [0.2, 0.1, 0.1, 0.1], 'String', "alpha1= desired fraction of unbound NAD(P)H")



%% 3D Plots Alpha1
subplot(1,2,1)
plot3(msec1,temp1k, a1_k1)
xlabel('Time (ms)')
ylabel('Temp(K)')
zlabel('Free NADH, alpha 1')
title('0.49 J/cm^2 (2.73 ms)', 'fontSize', fontSize)
set(gca,'FontSize',18);

subplot(1,2,2)
plot3(msec2,temp2k, a1_k2)
xlabel('Time (ms)')
ylabel('Temp(K)')
zlabel('Free NADH, alpha 1')
title('0.96 J/cm^2 (5.01 ms)', 'fontSize', fontSize)
set(gca,'FontSize',18);

%% 3D Plot Alpha 1 
 a1_k1=a1_k1*100;
 a1_k2=a1_k2*100;

plot3(msec1,temp1k, a1_k1,'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6,'DisplayName', '0.49 J/cm^2 (2.73 ms)',  'LineWidth', 1.5);
hold on
plot3(msec2,temp2k, a1_k2,'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', '0.96 J/cm^2 (5.01 ms)', 'LineWidth', 1.5);
grid off


fontSize= 16;
xlabel('Time (ms)', 'fontSize', fontSize)
ylabel('Temp (K)', 'fontSize', fontSize)
zlabel('Free NADH (%)', 'fontSize', fontSize)
legend( 'Location', 'northeast')
set(gca, 'LineWidth', 1);           % Axis line thickness
set(gca, 'TickLength', [0.01, 0.01]);  % Tick mark length
hold off
box on 


%% Function to Perform Exponential Decay Fitting
function [fitResult, msec_sub, temp_sub] = fitExponentialDecay(msec, temp, expModel)
    % Find index of maximum temperature
    [maxTemp, maxIdx] = max(temp);  
    
    % Select only data from the maximum temperature onward
    msec_sub = msec(maxIdx:end);
    temp_sub = temp(maxIdx:end);
    
    % Set initial parameter guesses
    Tmin = min(temp_sub); % Baseline (final temperature)
    A_init = maxTemp - Tmin; % Initial amplitude
    tau_init = (max(msec_sub) - min(msec_sub)) / 3; % Rough estimate for decay rate
    
    startPoints = [A_init, tau_init, Tmin];

    % Fit the model to the subset data
    fitResult = fit(msec_sub, temp_sub, expModel, 'Start', startPoints);
end

