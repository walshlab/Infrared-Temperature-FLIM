%% Import and plot original image w/0 artifact at the end 

% Uploading File from SPCimage. (X spatial by Y spatial by Time)
tpsfImage = Convert_ASC_Image('F:\Paper 1\Dynamic Temp\MCF-7\2 MCF7\test_2024_12_26_18_51_binned_raw_data.asc');
irfFile= readmatrix('IRF.txt'); %Loads IRF data for deconvolution 

% Artifact produced at the right-hand side, where no viable information is
% Sets artifact pixels as 0 (y positions from 249-256)
tpsfImage(:,249:256,:)= 0; 
%tpsfImage(:,196:256,:)= 0; 

sumImage = sum(tpsfImage,3); %Sums up dimensions so the intensity image can be created
imagesc(sumImage) %Plots intensities 


%% Normalize the data for visualization in imageSegmenter and get the first mask 
% Check the intensity range of the summed image
minIntensity = min(sumImage(:));
maxIntensity = max(sumImage(:));

disp(['Min Intensity: ', num2str(minIntensity)]);
disp(['Max Intensity: ', num2str(maxIntensity)]);

% Normalize the image intensity to the range [0, 1]
normalizedImage = sumImage / maxIntensity;
imageSegmenter(normalizedImage);  

%% Plotting Mask 
croppedImage = maskedImage(:, [1:248]);
imagesc(croppedImage);
%imagesc(maskedImage)
b= colorbar
colormap('parula')
ylabel(b, 'Fluorescence Intensity (a.u)', 'FontSize', 18)
axis image 
set(gca, 'XTick', [], 'YTick', []); % This removes the tick marks and labels

clear b 


%% Colorbar only 
figure;
imagesc([]);            % creates an empty axes
colormap('parula');
b = colorbar;
ylabel(b, 'Fluorescence Intensity (a.u)', 'FontSize', 18);
axis off;               % hides the axes box and ticks
%% Generate 2D array for timing of laser dwell time (Anna added)
% sumImage= normalizedImage *maxIntensity; %Unnormalizes the intensities 
scan_time=[]; 
time=1; %timing starts at 1 

for r=1:height(BW) % Going through each row
    for c=1:width(BW) % Going through each column
            scan_time(r,c)=time* 0.100; %microsecond dwell time of 100 microseconds 
         time=time+1;
    end
    % % flipping row, b/c laser scanning snakes through image (BIDIRECTIONAL) 
    % if rem(r,2)==0 %if remainder is two (for every alt. row, to flip)
    %     scan_time(r,:)=flip(scan_time(r,:));
    % 
    % end
end 

% Pinpoint timing of region of interest (ROI)
masked_timing= scan_time.*BW; %Filters to have only the time for when there is a cell 


clear r c time



% %% Taking into account the bidirectional gathering of data 
% % Iterate through each slice in the third dimension
% for k = 1:size(tpsfImage, 3)
%     for r = 2:2:size(tpsfImage, 1) % Only even rows
%         tpsfImage(r, :, k) = flip(tpsfImage(r, :, k)); % Flip each even row
%     end
% end

%% Finding the row and column values for positions= 1 (BW mask) Jocelyn 
% Preallocate a cell array to store indices for each row
valOnes  = cell(256, 1);  % Makes cell array for each row 

for r = 1:256 % Loop through each row of the array
    columnindex = find(BW(r, :) ~= 0);  % Find the column indices of non-zero elements in the current row
    % Store the row and corresponding column indices
    valOnes{r} = [repmat(r, length(columnindex), 1), columnindex'];  % Create a matrix of (row, col) pairs
end

clear columnindex r 

%% Getting decay information and time stored (NO BINNING)
% First two columns for row and column indices, next column for information

noBinDecay = []; % Matrix to store decay data without binning 

noBinTime = []; %Matrix to store time corresponding to decays without binning 

for r = 1:256
    if ~isempty(valOnes{r})  % Check if valOnes contains any non-zero elements for this row
        rowColPairs = valOnes{r};  % Get the row-column pairs from valOnes
        
        for i = 1:size(rowColPairs, 1)
            row = rowColPairs(i, 1);  % Extract row index
            col = rowColPairs(i, 2);  % Extract column index
            
            % Scan_time for each(row, col)
            noBinTime = [noBinTime; row, col, scan_time(row, col)];  
            % Third dimension decay for each (row, col)
            thirdDimData = squeeze(tpsfImage(row, col, :))';  % 1x256 vector
            
            % Matrix with row, col, and decay information
            noBinDecay = [noBinDecay; row, col, thirdDimData];  % Row, Col, decay data
        end
    end
end

clear col i r row rowColPairs thirdDimData


%% Binning decay
binNumber = input('How many pixels do you want to average? \n Please type an integer: ');

% Extract the row indices from the first column of noBinDecay
rowIndices = noBinDecay(:, 1);

% Extract the column indices from the second column of noBinDecay
colIndices = noBinDecay(:, 2);

% Find unique row indices
uniqueRowIndices = unique(rowIndices);

% Initialize array to store the final binned data
binDecayTotal = [];

% Loop through each unique row index
for i = 1:length(uniqueRowIndices)
    % Find rows in noBinDecay that match the current unique row index
    matchingRows = noBinDecay(rowIndices == uniqueRowIndices(i), :);
    
    % Get the number of rows for the current unique row index
    numRows = size(matchingRows, 1);
    
    % Check if there are enough rows for the specified binNumber
    if numRows >= binNumber
        % Determine how many full groups of binNumber can be formed
        numGroups = floor(numRows / binNumber);  % Number of full groups of binNumber
        
        % Loop through each group
        for j = 1:numGroups
            % Define the rows for the current group of binNumber within matchingRows
            groupRows = (1 + (j-1) * binNumber):(j * binNumber);
            
            % Extract the current column indices for this group
            currentColIndices = matchingRows(groupRows, 2);  % Get the column indices from the second column
            
            % Sum the information columns for the current group (excluding the first two columns)
            summedInfo = sum(matchingRows(groupRows, 3:end), 1);  % Summing columns 3 to end  
            
            % Store the row index, summed information, and the last column index
            binDecayTotal = [binDecayTotal; uniqueRowIndices(i), summedInfo, currentColIndices(end)];
        end
    end
end

binnedDecay= binDecayTotal(:, 2:257);

indexBinnedDecay=  binDecayTotal(:, [1, 258]);

clear colIndices currentColIndices groupRows i j matchingRows numGroups numRows rowIndices summedInfo uniqueRowIndices
clear binDecayTotal

%% Binning time data 
binNumber = input('How many pixels do you want to average? \n Please type an integer: ');

% Extract the row indices from the first column of noBinTime
rowIndices = noBinTime(:, 1);

% Extract the column indices from the second column of noBinTime
colIndices = noBinTime(:, 2);

% Find unique row indices
uniqueRowIndices = unique(rowIndices);

% Initialize an array to store the final binned data and corresponding times
binTimeTotal = [];

% Loop through each unique row index
for i = 1:length(uniqueRowIndices)
    % Find rows in noBinTime that match the current unique row index
    matchingRows = noBinTime(rowIndices == uniqueRowIndices(i), :);
    
    % Get the number of rows for the current unique row index
    numRows = size(matchingRows, 1);
    
    % Check if there are enough rows for the specified binNumber
    if numRows >= binNumber
        % Determine how many full groups of binNumber can be formed
        numGroups = floor(numRows / binNumber);  % Number of full groups of binNumber
        
        % Loop through each group
        for j = 1:numGroups
            % Define the rows for the current group of binNumber
            groupRows = (1 + (j-1) * binNumber):(j * binNumber);

            % Extract the current column indices for this group
            currentColIndices = matchingRows(groupRows, 2);  % Get the column indices from the second column
            
            % Extract the last row for this group
            lastRow = matchingRows(groupRows(end), :);  % Get the last row of the group
            
            % Store the row index and the last values of the group
            binTimeTotal = [binTimeTotal; lastRow, currentColIndices(end)];
        end
    end
end
binnedTime= binTimeTotal(:, 3);

indexBinnedTime=  binTimeTotal(:, [1,2]);

clear colIndices currentColIndices groupRows i j lastRow matchingRows numGroups numRows rowIndices uniqueRowIndices


%% Normalizes the IRF and produces the alpha 1 Values 
irfNorm = irfFile/max(irfFile); %Normalizes the IRF File 

size2DbinnedDecay= size(binnedDecay);
nsize= size2DbinnedDecay(1);

% [a1,t1,a2,t2] = Fit_Deconvolution(newdata(1439,:),irfNorm(68:92));

alpha1= zeros([nsize,1]);%Single alpha1 array of nsize
tao1= zeros([nsize, 1]);%Single tao1 array of nsize
alpha2= zeros([nsize,1]);%Single alpha2 array of nsize 
tao2= zeros([nsize,1]); %Single tao2 array of nsize 
p= 1; %index p starts at 1

%Produces alpha1 values for all data points
 for p= 1:nsize
     [alpha1(p,1), tao1(p,1), alpha2(p,1), tao2(p,1)]= Fit_Deconvolution(binnedDecay(p,:),irfNorm(68:92));
 end
 
 clear p size2DbinnedDecay nsize

%% Calculating alpha1 moving average
binNumber= 8; 
alpha1= alpha1*100; 
alpha1_mean= movmean(alpha1, binNumber*3); 

%% Arranging alpha1 for Heat Map 
indexalpha1= [indexBinnedDecay, alpha1_mean];
data = indexalpha1;  % Using indexalpha1 

% Create a 256x256 array initialized with zeros
heat_alpha1= zeros(256, 256); %Empty

% Loop through each row of the 216x3 array to populate resultArray
for i = 1:size(data, 1)
    rowIndex = data(i, 1);    % Row index from the first column
    colIndex = data(i, 2);    % Column index from the second column
    infoValue = data(i, 3);    % Information from the third column
    
    % Store the information in the corresponding position in the result array
    heat_alpha1(rowIndex, colIndex) = infoValue;
    
    % Copy the value to previous columns based on binNumber - 1
    for b = 1:(binNumber - 1)
        previousCol = colIndex - b;  % Calculate the previous column index
        % Check if the previous column index is valid
        if previousCol >= 1
            heat_alpha1(rowIndex, previousCol) = infoValue;
        end
    end
end

clear data rowIndex b i previousCol 

%% Heat map of alpha1 PLOT

heat_alpha1(heat_alpha1 ==0 )= NaN; %Sets values of 0 to NaN
figure;
imagesc(heat_alpha1);
c= colorbar; % Add a color bar to indicate the scale
ylabel(c, 'Free NADH (\alpha_1)', 'FontSize', 18)

caxis([min(heat_alpha1(:)) max(heat_alpha1(:))]);
axis image
% set(gca, 'XTick', [], 'YTick', []); % This removes the tick marks and labels
% Adjust colormap
colormap('parula'); % Choose your desired colormap
% Highlight a specific row (e.g., row 5) in red
artifact= 80; 
 highlightRow = artifact;
hold on;
plot(1:size(heat_alpha1, 2), highlightRow * ones(1, size(heat_alpha1, 2)), 'w', 'LineWidth', 1);  % Plot a red line across the row

clear highlightRow c

%% Plotting Free NADH response over Time 
artifact= 52;
index = find(binTimeTotal(:,1) == artifact); %Finds the index where the artifact starts 
index_Artifact= index(1); %First occurence of the artifact 
[max_Value, max_Index ]= max(alpha1_mean(index_Artifact:end)); %Max value and index of where free NADH is highest


laserOnIndex = index_Artifact; % Index where the laser was turned on

% Calculate the response value (mean after laser is turned on)
responsePeriod = alpha1_mean(laserOnIndex:end); % Data after laser turned on

% Plot the time series data
figure;
p1= plot(binnedTime, alpha1) %No smoothing
hold on 
%p2= plot(binnedTime, alpha1_mean,'b', 'LineWidth', 3) %alpha1  smoothing for CONTROL

%Highlight the baseline period
p2= plot(binnedTime(1:laserOnIndex-1), alpha1_mean(1:laserOnIndex-1), 'b', 'LineWidth', 3); %Plots baseline period
% 
% %Highlight the response period
p3= plot(binnedTime(laserOnIndex:end), alpha1_mean(laserOnIndex:end), 'r', 'LineWidth', 3); %Plots response period 
% 
% %Add a vertical line to indicate where the laser was turned on
xline(binnedTime(laserOnIndex), '--k', 'LineWidth', 2); % Add a vertical dashed line for when laser starts 

xlabel('Time (ms)', 'FontSize', 28)
ylabel('Free NAD(P)H (%)', 'FontSize', 28)
xlim([0 7000])
xticks(0:2000:7000)

ylim([50 100])
yticks(40:15:100)
set(gca, 'FontSize', 23)
% Display the legend (automatically gathers DisplayNames from each plot)
lgd= legend({'Experiment Data', 'Baseline Period','Response Period'}, 'Location', 'northeast');
lgd.FontSize = 25;

%% Plotting Free NADH response over Time CONTROL

% Plot the time series data
figure;
p1= plot(binnedTime, alpha1) %No smoothing
hold on 
p2= plot(binnedTime, alpha1_mean,'b', 'LineWidth', 3) %alpha1  smoothing for CONTROL


xlabel('Time (ms)', 'FontSize', 28)
ylabel('Free NAD(P)H (%)', 'FontSize', 28)

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
plot(binnedTime, percentChange,'b', 'LineWidth', 2);

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


%% Percent change plot use this one 

% startIndex= index_Artifact-10;
startIndex= 1; %We are starting at index 1 
index_Artifact= laserOnIndex; %This is the index for our artifact 
%baselineValue= mean(alpha1_mean); for control stuff 
baselineValue= mean(alpha1_mean(1:index_Artifact-1));
percentChange= ((alpha1_mean-baselineValue)/baselineValue)*100; %Computes percent change 

% Plot the percent change
figure;
plot(binnedTime(startIndex:end), percentChange(startIndex:end));
hold on 
% Highlight the baseline period
plot(binnedTime(startIndex:laserOnIndex-1), percentChange(startIndex:laserOnIndex-1), 'b', 'LineWidth', 2);

% Highlight the response period
plot(binnedTime(laserOnIndex:end), percentChange(laserOnIndex:end), 'r', 'LineWidth', 2);

% Add a vertical line to indicate where the laser was turned on
xline(binnedTime(laserOnIndex), '--k', 'LineWidth', 2); % Add a vertical dashed line

% title('Percent Change of Free NADH (\alpha_1)', 'FontSize', 18)
xlabel('Time (ms)', 'FontSize', 28)
ylabel('Percent Change (%)', 'FontSize', 28)
% Optionally, you can add a legend
lgd= legend('Percent Change', 'Baseline Period', 'Response Period', 'Location', 'northeast');
lgd.FontSize = 25;
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

%% Bins decays in groups of 8 for solutions 
nsize = 65536/8; %Dividing decay data into groups of 8
newdata = zeros([nsize,256]); %Empty array of nsize (8192 summed groups) x 256
time = zeros([nsize,1]); %Single time array of nsize (8192 summed groups) 
i = 1; %Index i starts at 1 
j = 1; %Index j starts at 1 
a=0;

while i < 65530
    newdata(j,:)=sum(test202212051506binnedrawdata(i:i+7,:),1); %Sums 16 X-val with corresponding y into a 1 dimensional array 
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

decayTime = 0:12.5/255:12.5; %In nanoseconds for 1 decay curve 
 
%% Normalizes the IRF and produces the alpha 1 Values
% This is the instrument response function, which is something unique to 
% each system, can be obtained imaging urea crystals or YG beads

irfFile= readmatrix('IRF.txt'); %Loads IRF data
irfNorm = irfFile/max(irfFile);

size2DbinnedDecay= size(binnedDecay); %For cells do not need for soln 
nsize= size2DbinnedDecay(1); %For cells do not need for soln

alpha1= zeros([nsize,1]);%Single alpha1 array of nsize (8192 summed groups)
tao1= zeros([nsize, 1]);%Single tao1 array of nsize (8192 summed groups)
alpha2= zeros([nsize,1]);%Single alpha2 array of nsize (8192 summed groups)
tao2= zeros([nsize,1]); %Single tao2 array of nsize (8192 summed groups)
p= 1; %index p starts at 1
 

% for p= 1:nsize
%      [alpha1(p,1), tao1(p,1), alpha2(p,1), tao2(p,1)]= Fit_Deconvolution_v1(newdata(p,:),irfNorm(68:92));
%  end


 for p= 1:nsize
     [alpha1(p,1), tao1(p,1), alpha2(p,1), tao2(p,1)]= Fit_Deconvolution_v1(binnedDecay(p,:),irfNorm(68:92));
 end


%% Function to deconvolute and fit data/ Function to open up asc file 
function [a1Norm,t1,a2Norm,t2] = Fit_Deconvolution(decayCurve,irf)
    timeFitting = [1:256]*12.5/256; %Nanoseconds 
    if max(decayCurve)> 0
            deconvDecay = deconvlucy(decayCurve',irf,10); %Deconvolutes data
            [~, peakPosition] = max(deconvDecay); %Finds the max peak of decay curve
            fittingPart = [deconvDecay(peakPosition:256)' zeros(1,peakPosition-1)]; 
            
            %Conditions are needed to ensure alpha1 is never above 1 
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
               % f = fittype('a*exp(b*x) + c*exp(d*x)')
               % c = cfit(f,a1,-1/t1,a2,-1/t2)
               % plot(c,timeFitting,fittingPart)
               %  xlabel('Time (ns)');
               %  ylabel('Decay Amplitude');
               %  title('Double Exponential Fit');
               %  legend('Fitted Curve', 'Decay Curve');
               %  grid on;
end


%Function to open up asc file 
function photonDecayImage3D = Convert_ASC_Image(imageDirectory) %Function to read ASC File 
PhotonRawDecayImage = importdata(imageDirectory,' ',11);
PhotonDecayImage = PhotonRawDecayImage.data;
photonDecayImage3D =  reshape(PhotonDecayImage,[256 256 256]);
photonDecayImage3D = permute(photonDecayImage3D,[2 1 3]);
end

