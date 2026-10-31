%% Plot Transmission by division of envelopes
 
%% clean up
clear all; close all; clc;

%% Font, LineWidth and MarkerSize
FontSize = 16;
MarkerSize = 3;
LineWidth = 1;

%% Collect data
n = 2; % number of columns

%% Read in Reference Data
numFileBase = 1; % number of base files to be averaged
cdExtBase = ''; % file extention for base
%cdExtBase = 'Air_Ref_90\'; % file extention for base

for i=1:numFileBase
    txtNameRUN = [cdExtBase, 'air_0001.txt'];% name of THz file(ref)
    fibRUN = fopen(txtNameRUN, 'r');% opens file
        fscanf(fibRUN, 'Frequency (GHz)	THz Photocurrent (nA)');% first line of file
        [dataRUN, amountRUN] = fscanf(fibRUN,'%f', [n,inf]);% scans in data   
    if diff(dataRUN(1,:)) > 0
        xRUNbase(i,:) = dataRUN(1,:);
        yRUNbase(i,:) = dataRUN(2,:);
    else
        xRUNbase(i,:) = flip(dataRUN(1,:));
        yRUNbase(i,:) = flip(dataRUN(2,:));
    end
end

%% Read in Sample Data
numFileSamp = 1; % number of sample files to be averaged
cdExtSamp = ''; % file extention for sample
%cdExtSamp = 'Ga2O3_90flat_TE\'; % file extention for sample

for i=1:numFileSamp     
    txtNameRUN=[cdExtSamp,'200-1000 GHz-Shadi-SRR-GaO-backplane-0deg-annealed_0001_0001.txt'];
    fibRUN = fopen(txtNameRUN, 'r');% opens file
        fscanf(fibRUN,'Frequency (GHz)	THz Photocurrent (nA)');
        [dataRUN,amountRUN]=fscanf(fibRUN,'%f', [n,inf]);% scans in data  
    if diff(dataRUN(1,:)) > 0
        xRUNsamp(i,:) = dataRUN(1,:);% used 
        yRUNsamp(i,:) = dataRUN(2,:);% used 
    else
        xRUNsamp(i,:) = flip(dataRUN(1,:));% used 
        yRUNsamp(i,:) = flip(dataRUN(2,:));% used 
    end
end
        
fclose('all');% closes file
    
%% Average Multiple Scans
 xRUNbaseAve = xRUNbase(1,:);
 yRUNbaseAve = mean(yRUNbase,1);
 xRUNsampAve = xRUNsamp(1,:);
 yRUNsampAve = mean(yRUNsamp,1);
      
%% Verify Averaged Data
% figure
% hold on
% for i=1:numFileBase
%     plot(xRUNbase(i,:),yRUNbase(i,:))
% end
% plot(xRUNbaseAve,yRUNbaseAve)
% hold off

%% Interpolating the photocurrents for phase estimation
x_base_int = xRUNbaseAve(1):1E-3:xRUNbaseAve(length(xRUNbaseAve));
x_samp_int = xRUNsampAve(1):1E-3:xRUNsampAve(length(xRUNsampAve));
y_base_int = spline(xRUNbaseAve,yRUNbaseAve,x_base_int);
y_samp_int = spline(xRUNsampAve,yRUNsampAve,x_samp_int);

figure(1)
plot(x_base_int,y_base_int)
hold on
plot(x_samp_int,y_samp_int)

%% No interpolation prior to finding envelope 
% Finding envelope of sample and reference scan
[est_base_high, est_base_low] = envelope(yRUNbaseAve,30,'rms');
[est_samp_high, est_samp_low] = envelope(yRUNsampAve,30,'rms');

% Plotting envelopes of sample and reference
figure(2)
hold on
    plot(xRUNbaseAve,est_base_high,'.',xRUNbaseAve,est_base_low,'.','MarkerSize',3)
    plot(xRUNsampAve,est_samp_high,'.',xRUNsampAve,est_samp_low,'.','MarkerSize',3)
    title('Calculated Envelopes of THz Photocurrent','FontSize',FontSize)
    xlabel(' Frequency (GHz)','FontSize',FontSize)
    ylabel(' THz Photocurrent (nA)','FontSize',FontSize)
    set(gca,'FontSize',FontSize)
    %plot(xENVbaseAve,movmean(xENVsampAve./yRUNbaseAve,50),'b')
hold off

%% Plotting measured photocurrents with envelopes
% Plotting photocurrents and their envelopes overlaid
figure(3)
hold on
    plot(xRUNbaseAve,yRUNbaseAve,'b',xRUNbaseAve,est_base_high,'k')
    plot(xRUNsampAve,yRUNsampAve,'r',xRUNsampAve,est_samp_high,'k')
    legend('Reference PC','Reference Envelope','Sample PC','Sample Envelope')
    title('Measured THz Photocurrents and Calculated Envelopes','FontSize',FontSize)
    xlabel('Frequency (GHz)','FontSize', FontSize)
    ylabel('THz Photocurrent (nA)','FontSize', FontSize)
    set(gca,'FontSize',FontSize)
hold off

%% Finding THz-Transmission from envelopes 
TxBase = xRUNbaseAve;
TyBase = est_base_high;
TxSamp = xRUNsampAve;
TySamp = est_samp_high;

if length(TxBase) > length(TxSamp)
    xtot = length(TxSamp);
else
    xtot = length(TxBase);
end

Tx = ((TxBase(1:xtot)+TxSamp(1:xtot))./2);
Ty = ((TySamp(1:xtot)./TyBase(1:xtot))); % Square the calculated transmission % equation 5 of paper

TyMean = movmean(Ty,500); 

%% Error Bars

%to make error bars and scatter points
c=[56/255,61/255,150/255; %blue
   175/255,54/255,60/255; %red
   52/255,52/255,52/255; %black
   70/255,148/255,73/255; %green
   94/255,60/255,108/255; %purple
   115/255,82/255,68/255; %brown
   8/255,133/255,161/255; %cyan
   87/255,108/255,67/255; %dark green
   % 214/255,126/255,44/255; %orange
   % 187/255,86/255,149/255; %magenta
    ]; 

%% Plot THz data
m=1;
x = Tx;
xLimit = [min(min(x)), max(max(x))];
%xLimit = [0, 50];
xTitle = 'Frequency (GHz)';

y = Ty;
% yLimit = [min(min(y)), max(max(y))];
yLimit = [0, 1.5];
yTitle = 'Transmission';

Legend={'Calculated';
        'Moving Mean';
        };

% plot
figure1=figure(4);
axes1 = axes('Parent',figure1,'FontSize',FontSize); % creates axes in arbitrary positions
xlim(axes1,xLimit); % sets a limit to x-axis
ylim(axes1,yLimit); % sets a limit to x-axis

hold on
for k = 1:m
    plot(x(k,(1:1:end)),(y(k,(1:1:end))),'.-','MarkerSize',MarkerSize,'Color', [c(k,1),c(k,2),c(k,3)]); % plot
   
end
hold off
%% GSO(110)\_Substrat\_Theta100, TOO\_36\_17\_DSP\_MAS\_C\_1\_Theta100
Title = 'Ga2O3-Small-Sample\_90degree';
%% Labeling
title(Title,'FontSize',FontSize)
ylabel(yTitle,'FontSize',FontSize)
xlabel(xTitle,'FontSize',FontSize)
   
%% Extras
box on
set(gca,'XMinorTick','on','YMinorTick','on')

save 0deg x y



% Export transmission to CSV
T = table(x(:), y(:), 'VariableNames', {'Frequency_GHz', 'Transmission'});
writetable(T, '0deg_transmission.csv');