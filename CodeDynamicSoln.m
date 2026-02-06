%% Import and plot original image w/0 artifact at the end 

% Uploading File from SPCimage. (X spatial by Y spatial by Time)
tpsfImage = Convert_ASC_Image('F:\Paper 1\Dynamic Temp\102723_dypulse\0.5 ms\test_2023_10_27_16_04_binned_raw_data.asc');

% Artifact produced at the right-hand side, where no viable information is
% Sets artifact pixels as 0 (y positions from 249-256
tpsfImage(:,249:256,:)= 0; 

sumImage = sum(tpsfImage,3); 

minIntensity = min(sumImage(:));
maxIntensity = max(sumImage(:));

disp(['Min Intensity: ', num2str(minIntensity)]);
disp(['Max Intensity: ', num2str(maxIntensity)]);

% Normalize the image intensity to the range [0, 1]
normalizedImage = sumImage / maxIntensity;
imagesc(normalizedImage)
b= colorbar
colormap('parula')
ylabel(b, 'Fluorescence Intensity (a.u)', 'FontSize', 18)
axis image 
%set(gca, 'XTick', [], 'YTick', []); % This removes the tick marks and labels

%% Normalizes the IRF and produces the alpha 1 Values
% This is the instrument response function, which is something unique to 
% each system, can be obtained imaging urea crystals or YG beads

irfFile= readmatrix('IRF.txt'); %Loads IRF data
irfNorm = irfFile/max(irfFile);
 [a1,t1,a2,t2] = Fit_Deconvolution(newdata(1439,:),irfNorm(68:92));

alpha1= zeros([10,1]);%Single alpha1 array of nsize (8192 summed groups)
tao1= zeros([nsize, 1]);%Single tao1 array of nsize (8192 summed groups)
alpha2= zeros([nsize,1]);%Single alpha2 array of nsize (8192 summed groups)
tao2= zeros([nsize,1]); %Single tao2 array of nsize (8192 summed groups)
p= 1; %index p starts at 1

%Produces alpha1 values for all data points
 for p= 1:nsize
     [alpha1(p,1), tao1(p,1), alpha2(p,1), tao2(p,1)]= Fit_Deconvolution(newdata(p,:),irfNorm(68:92));
 end

 plot(entireTime(1:400), alpha1)
 xlabel('Time (ms)')
 ylabel('Free NADH')

%% Bins decays in groups of 8
nsize = 65536/8; %Dividing decay data into groups of 8
newdata = zeros([nsize,256]); %Empty array of nsize (8192 summed groups) x 256
time = zeros([nsize,1]); %Single time array of nsize (8192 summed groups) 
i = 1; %Index i starts at 1 
j = 1; %Index j starts at 1 
a=0;

while i < 65530
    newdata(j,:)=sum(test202310271604binnedrawdata(i:i+7,:),1); %Sums 16 X-val with corresponding y into a 1 dimensional array 
    if rem((i-256*a),241)==0; %Filters out the last bit of artifact on image 
        time(j,1)=time(j-1,1)+0.1*2*8; 
        i=i+8+8;
        a=a+1;
    else
        if j==1;
            time(j,1)=0;
        else
            time(j,1)=time(j-1,1)+0.1*8;  %In milliseconds
        end
        i=i+8; %Adds 16 to i

    end
    j=j+1; %Adds 1 to j 
end


entireTime=time; %Time for entire image (ms)

decayTime = 0:12.5/255:12.5; 


%% Calculating alpha1 moving average
binNumber= 8; 
alpha1= alpha1*100; 
alpha1_mean= movmean(alpha1, binNumber*3); 

%% Plotting Free NADH response over Time 
%artifact= 62;
target= 1561.6 %Time for first pixel of artifact 43
[~, index]= min(abs(entireTime-target));
laserOnIndex = index; % Index where the laser was turned on

% Calculate the response value (mean after laser is turned on)
responsePeriod = alpha1_mean(laserOnIndex:end); % Data after laser turned on

% Plot the time series data
figure;
p1= plot(entireTime, alpha1) %No smoothing
hold on 
%Highlight the baseline period
p2= plot(entireTime(1:laserOnIndex-1), alpha1_mean(1:laserOnIndex-1), 'b', 'LineWidth', 3); %Plots baseline period
% 
% %Highlight the response period
p3= plot(entireTime(laserOnIndex:end), alpha1_mean(laserOnIndex:end), 'r', 'LineWidth', 3); %Plots response period 
% 
% %Add a vertical line to indicate where the laser was turned on
xline(entireTime(laserOnIndex), '--k', 'LineWidth', 2); % Add a vertical dashed line for when laser starts 

xlabel('Time (ms)', 'FontSize', 28)
ylabel('Free NADH (%)', 'FontSize', 28)

xlim([0 7000])
xticks(0:2000:7000)

ylim([50 100])
yticks(40:15:100)
set(gca, 'FontSize', 23)
% Display the legend (automatically gathers DisplayNames from each plot)
lgd= legend({'Experiment Data', 'Baseline Period','Response Period'}, 'Location', 'northeast');
lgd.FontSize = 22;


%% Percent change plot use this one 

% startIndex= index_Artifact-10;
startIndex= 1; %We are starting at index 1 
index_Artifact= laserOnIndex; %This is the index for our artifact 
%baselineValue= mean(alpha1_mean); for control stuff 
baselineValue= mean(alpha1_mean(1:index_Artifact-1));
percentChange= ((alpha1_mean-baselineValue)/baselineValue)*100; %Computes percent change 

% Plot the percent change
figure;
plot(entireTime(startIndex:end), percentChange(startIndex:end));
hold on 
% Highlight the baseline period
plot(entireTime(startIndex:laserOnIndex-1), percentChange(startIndex:laserOnIndex-1), 'b', 'LineWidth', 2);

% Highlight the response period
plot(entireTime(laserOnIndex:end), percentChange(laserOnIndex:end), 'r', 'LineWidth', 2);

% Add a vertical line to indicate where the laser was turned on
xline(entireTime(laserOnIndex), '--k', 'LineWidth', 2); % Add a vertical dashed line

% title('Percent Change of Free NADH (\alpha_1)', 'FontSize', 18)
xlabel('Time (ms)', 'FontSize', 28)
ylabel('Percent Change (%)', 'FontSize', 28)
% Optionally, you can add a legend
lgd= legend('Percent Change', 'Baseline Period', 'Response Period', 'Location', 'northeast');
lgd.FontSize = 22;
% grid on;
xlim([0 7000])
xticks(0:2000:7000)

ylim([-10 70])
yticks(-10:25:70)
set(gca, 'FontSize', 23)
% max_Value= max(alpha1_mean(index_Artifact:index_Artifact+500))
max_Value= max(alpha1_mean(index_Artifact:end))
totalIncrease= max_Value- baselineValue
 percentChangeOverall= ((max_Value-baselineValue)/baselineValue)*100
 fprintf('Baseline value: %.3f\n', baselineValue)


%% Plotting Free NADH response over Time CONTROL

% Plot the time series data
figure;
p1= plot(entireTime, alpha1) %No smoothing
hold on 
p2= plot(entireTime, alpha1_mean,'b', 'LineWidth', 3) %alpha1  smoothing for CONTROL


xlabel('Time (ms)', 'FontSize', 28)
ylabel('Free NADH (%)', 'FontSize', 28)
xlim([0 7000])
xticks(0:2000:7000)

ylim([50 100])
yticks(40:15:100)
set(gca, 'FontSize', 23)
% Display the legend (automatically gathers DisplayNames from each plot)
%legend({'Experiment Data', 'Baseline Period','Response Period'}, 'Location', 'Best');


%% Percent change plot CONTROL 

baselineValue= mean(alpha1_mean); %for control stuff 
percentChange= ((alpha1_mean-baselineValue)/baselineValue)*100; %Computes percent change 

% Plot the percent change
figure;
plot(entireTime, percentChange,'b', 'LineWidth', 2);

% title('Percent Change of Free NADH (\alpha_1)', 'FontSize', 18)
xlabel('Time (ms)', 'FontSize', 28)
ylabel('Percent Change (%)', 'FontSize', 28)
% Optionally, you can add a legend
% legend('Percent Change', 'Baseline Period', 'Response Period', 'Location', 'Best');
% grid on;
xlim([0 7000])
xticks(0:2000:7000)

ylim([-10 70])
yticks(-10:25:70)
set(gca, 'FontSize', 23)
%% Function to deconvolute and fit data 
function [a1Norm,t1,a2Norm,t2] = Fit_Deconvolution(decayCurve,irf)
    timeFitting = [1:256]*12.5/256; %Nanoseconds 
    if max(decayCurve)> 0
            deconvDecay = deconvlucy(decayCurve',irf,10); %Deconvolutes data
            [~, peakPosition] = max(deconvDecay); %Finds the max peak of decay curve
            fittingPart = [deconvDecay(peakPosition:256)' zeros(1,peakPosition-1)]; 
            decayFit = fit(timeFitting',fittingPart','exp2','StartPoint',[0.6,-4, 1, -0.8],'Lower',[0.5,-5, 0.03,-1],'Upper',[80,-1.67,60,-0.3]); %Fits data to a 2 exponential form  
            fittingValues = coeffvalues(decayFit); %Array of all coefficients 
            a1 = fittingValues(1);
            a1Norm = fittingValues(1)/(fittingValues(1) + fittingValues(3) );
            t1 = -1/fittingValues(2);
            a2 = fittingValues(3);
            a2Norm = fittingValues(3)/(fittingValues(1) + fittingValues(3) );
            t2 = -1/fittingValues(4);
    else
            a1Norm = 0;
            t1 = 0;
            a2Norm = 0;
            t2 = 0;
    end
%                f = fittype('a*exp(b*x) + c*exp(d*x)')
%                c = cfit(f,a1,-1/t1,a2,-1/t2)
%                plot(c,timeFitting,fittingPart)
end


%% Function for opening up file type in Matlab
function photonDecayImage3D = Convert_ASC_Image(imageDirectory) %Function to read ASC File 
PhotonRawDecayImage = importdata(imageDirectory,' ',11);
PhotonDecayImage = PhotonRawDecayImage.data;
photonDecayImage3D =  reshape(PhotonDecayImage,[256 256 256]);
photonDecayImage3D = permute(photonDecayImage3D,[2 1 3]);
end
