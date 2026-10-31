%% Measured THz transmission of a patterned sample (as deposited and annealed) against an HFSS conductivity sweep
% 200-800 GHz, 0 deg = TE. Edit the data folder and file names below.
% No new air scan with the annealed run, so both use the Oct 19 air reference.
clear all; close all; clc;
D = [pwd filesep];   % folder with the THz scans and the HFSS export
OUT = D;
FontSize = 16; LineWidth = 1;
c=[56/255,61/255,150/255;   %blue
   175/255,54/255,60/255;   %red
   52/255,52/255,52/255;    %black
   70/255,148/255,73/255;   %green
   94/255,60/255,108/255;   %purple
   8/255,133/255,161/255];  %cyan

%% measured: envelope ratio, same as plots_oct8.m
[f, air] = readscan([D '200-800 GHz-Shadi-air_0001_0001.txt']);
[f0, d0] = readscan([D '200-800 GHz-Shadi-SRR-front-GaO-0deg_0001_0001.txt']);
[f1, d1] = readscan([D '200-800 GHz-Shadi-SRR-front-GaO-0deg-annealled_0001_0001.txt']);
d0 = interp1(f0, d0, f); d1 = interp1(f1, d1, f);
[ea,~] = envelope(air,30,'rms'); [e0,~] = envelope(d0,30,'rms'); [e1,~] = envelope(d1,30,'rms');
T0 = e0./ea; T1 = e1./ea;

%% HFSS TE-T sweep (H = 686 um)
H = readmatrix([D 'TE-T-different cond-200-1000-front.csv']);
sig = [94000 150000 210000];
lab = {'94,000','150,000','210,000'};
sty = {'--',':','-.'};
col = [c(2,:); c(5,:); c(4,:)];
hf = cell(1,3); ht = cell(1,3);
for i = 1:3
    k = H(:,1) == sig(i);
    hf{i} = H(k,2)*1000; ht{i} = H(k,3);
end

%% 1) annealed vs as-deposited
figure('Position',[100 100 860 500]); hold on
plot(f, movmean(T0,50), 'Color', c(1,:), 'LineWidth', LineWidth)
plot(f, movmean(T1,50), 'Color', c(2,:), 'LineWidth', LineWidth)
hold off; box on; grid on; set(gca,'FontSize',FontSize,'XMinorTick','on','YMinorTick','on')
xlim([200 800]); ylim([0.2 0.85])
xlabel('Frequency (GHz)'); ylabel('Transmission')
legend('as-deposited 10/19','annealed 10/20','Location','northeast','FontSize',12)
exportgraphics(gcf,[OUT '200-800GHz_annealed_vs_asdeposited_Oct20.png'],'Resolution',250)

%% 2) annealed vs HFSS 210,000 S/m (best of the sweep)
figure('Position',[100 100 860 500]); hold on
plot(f, movmean(T1,50), 'Color', c(1,:), 'LineWidth', LineWidth)
plot(hf{3}, ht{3}, '--', 'Color', c(2,:), 'LineWidth', 1.8)
hold off; box on; grid on; set(gca,'FontSize',FontSize,'XMinorTick','on','YMinorTick','on')
xlim([200 800]); ylim([0.2 0.85])
xlabel('Frequency (GHz)'); ylabel('Transmission')
legend('annealed 10/20','HFSS, ITO 210,000 S/m','Location','northeast','FontSize',12)
exportgraphics(gcf,[OUT '200-800GHz_annealed_vs_HFSS_210k_Oct20.png'],'Resolution',250)

%% 3) fringe-averaged: annealed vs HFSS conductivities
nm = round(63/mean(diff(f)));
T1l = movmean(T1,nm); k = f>=240 & f<=760;
figure('Position',[100 100 860 500]); hold on
plot(f(k), T1l(k), 'Color', c(1,:), 'LineWidth', 1.8)
for i = 1:3
    nh = round(63/mean(diff(hf{i})));
    hl = movmean(ht{i},nh); kh = hf{i}>=240 & hf{i}<=760;
    plot(hf{i}(kh), hl(kh), sty{i}, 'Color', col(i,:), 'LineWidth', 1.8)
end
hold off; box on; grid on; set(gca,'FontSize',FontSize,'XMinorTick','on','YMinorTick','on')
xlim([200 800]); ylim([0.25 0.65])
xlabel('Frequency (GHz)'); ylabel('Transmission, fringe-averaged')
legend(['annealed 10/20'; strcat({'HFSS, ITO '}, lab', {' S/m'})], 'Location','northeast','FontSize',12)
exportgraphics(gcf,[OUT '200-800GHz_annealed_fringe-averaged_vs_HFSS_Oct20.png'],'Resolution',250)

%% 4) before and after anneal, each with its best HFSS conductivity
figure('Position',[100 100 860 500]); hold on
plot(f, movmean(T0,50), 'Color', c(1,:), 'LineWidth', LineWidth)
plot(f, movmean(T1,50), 'Color', c(2,:), 'LineWidth', LineWidth)
plot(hf{1}, ht{1}, '--', 'Color', c(1,:), 'LineWidth', 1.8)
plot(hf{3}, ht{3}, '--', 'Color', c(2,:), 'LineWidth', 1.8)
hold off; box on; grid on; set(gca,'FontSize',FontSize,'XMinorTick','on','YMinorTick','on')
xlim([200 800]); ylim([0.2 0.9])
xlabel('Frequency (GHz)'); ylabel('Transmission')
legend('before anneal (10/19)','after anneal (10/20)','HFSS, ITO 94,000 S/m','HFSS, ITO 210,000 S/m','Location','northeast','FontSize',12)
exportgraphics(gcf,[OUT '200-800GHz_before_after_anneal_vs_HFSS_Oct20.png'],'Resolution',250)

%% numbers for the slides
T0s = movmean(T0,50); T1s = movmean(T1,50);
bands = [200 300; 300 400; 400 500; 500 600; 600 700; 700 800];
for b = 1:6
    kb = f>=bands(b,1) & f<bands(b,2);
    fprintf('%d-%d GHz: as-deposited %.3f, annealed %.3f, change %+.3f\n', bands(b,1), bands(b,2), mean(T0s(kb)), mean(T1s(kb)), mean(T1s(kb))-mean(T0s(kb)));
end
for i = 1:3
    kk = f>=210 & f<=790;
    r = sqrt(mean((T1s(kk) - interp1(hf{i},ht{i},f(kk))).^2));
    nh = round(63/mean(diff(hf{i}))); hl = movmean(ht{i},nh);
    kl = f>=240 & f<=760; dd = T1l(kl) - interp1(hf{i},hl,f(kl));
    fprintf('annealed vs ITO %6d S/m: rms 210-790 GHz = %.3f | fringe-averaged mean(measured - HFSS) = %+.3f\n', sig(i), r, mean(dd));
end

function [x,y] = readscan(fn)
    fid = fopen(fn,'r'); fgetl(fid); d = fscanf(fid,'%f',[2 inf]); fclose(fid);
    if d(1,2) < d(1,1), d = fliplr(d); end
    x = d(1,:)'; y = d(2,:)';
end
