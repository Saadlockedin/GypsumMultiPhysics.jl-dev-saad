% MATLAB script to analyze TGA data and plot weight (%) vs. temperature
clear; close all; clc;

% Define the path to the directory containing the text files
folderPath = 'C:\Users\saada\Desktop\GSU\GypsumMultiPhysics.jl-dev-saad\TGAdata\Bl_Moisture and Mold Resistant Drywall';

% Define the pattern for the files you want to read
filePattern = fullfile(folderPath, 'Gyp_*_BI_T1_Nitrogen.txt');

% Get a list of all files matching the pattern
fileList = dir(filePattern);

% Extract the numeric part from the filenames and sort the files
fileNumbers = arrayfun(@(x) sscanf(x.name, 'Gyp_%d_GW_T2_Nitrogen.txt'), fileList);
[~, sortedIdx] = sort(fileNumbers);
sortedFileList = fileList(sortedIdx);

    
% Initialize a cell array to store the legend entries
legendEntries = {};

% Initialize a figure
figure(1);
hold on;
figure(2);
hold on;

% Define colors for the plots (you can customize this as needed)
colors = lines(length(sortedFileList));
%k=1
for k = 1:length(sortedFileList)
    % Extract the filename without path and extension for legend and image name
    [~, name, ~] = fileparts(sortedFileList(k).name);

    filename = fullfile(folderPath, sortedFileList(k).name);
    
    % Read the data from the text file
    opts = detectImportOptions(filename, 'FileType', 'text', 'Delimiter', '\t');
    data = readtable(filename, opts);
    
    % Extract values from row 43 to end, and columns 1 to 6 only
    data = data(42:end, 1:6);
    
    % Label the columns
    data.Properties.VariableNames = {'Time_min', 'Temperature_C', 'Weight_mg', 'Balance_Purge_Flow_mL_min', 'Sample_Purge_Flow_mL_min', 'Deriv_Weight_percent_C'};
    
    % Extract relevant columns
    temperature = data.Temperature_C;
    weight_mg = data.Weight_mg;
    Deriv_Weight = data.Deriv_Weight_percent_C;

    % Calculate the initial weight and the weight percentage
    initial_weight = weight_mg(1); % assuming the initial weight is the first entry
    % weight_percentage = weight_mg
    weight_percentage = (weight_mg / initial_weight) * 100;
    
    % Plot weight (%) vs. temperature
    % yyaxis left; % Use the left y-axis
    
    % Add the current file name to the legend entries
    legendEntries{end + 1} = name;

    figure(1);
    plot(temperature, weight_percentage, 'Color', colors(k, :), 'LineWidth', 1.5, 'MarkerSize',0.1);
    ylabel('Weight (%)');
    xlabel('Temperature (°C)');
    legend(legendEntries, 'Interpreter', 'none', 'Location', 'northeast');
    grid on;
    
    figure(2);
    % yyaxis right; % Use the right y-axis
    plot(temperature, Deriv_Weight, 'Color', colors(k, :), 'LineWidth', 1.5, 'MarkerSize',0.1);
    ylabel('Deriv. Weight (%/°C)');
    xlabel('Temperature (°C)');
    legend(legendEntries, 'Interpreter', 'none', 'Location', 'northeast');
    grid on;
    
end

hold off;

% Save the plot as an image (optional, uncomment to use)
% saveas(gcf, 'Combined_TGA_Plot.png');
