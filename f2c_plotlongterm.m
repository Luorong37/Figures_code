%% load movie and ROI
data_path = 'I:\1_Data\3f. Long-term imaging\20230810' ;
load('I:\1_Data\3f. Long-term imaging\20230810-161506recordSCN_Analysis\1_raw_ROI.mat')

data_dir = dir(data_path);
dir_name = {data_dir([data_dir.isdir]).name};
dir_name = dir_name(~ismember(dir_name,{'.','..'}));
traces_original = [];
traces_bgfitcorr = [];
background = zeros([size(bwmask),length(dir_name)]);
rois = struct();
rois.bwmask = bwmask;

% load traces
for i = 1:length(dir_name)
    fprintf('%d for %d\n',i,length(dir_name))
    current_folder = fullfile(data_path,dir_name{i});
    [movie, ncols, nrows, nframes] = load_movie(current_folder);
    % load trace
    % [~, current_traces] = select_ROI(movie, nrows, ncols, bwmask, []); 
    % traces = [traces;current_traces];

    % load background
    [current_background, current_background_fitted, current_traces, current_traces_bgfitcorr, background_mask] = remove_background(movie, ncols, nrows, rois);
    traces_original = [traces_original;current_traces];
    traces_bgfitcorr = [traces_bgfitcorr;current_traces_bgfitcorr];

    % save background mask
    background(:,:,i) = background_mask;
    close();
end

cutindex = 722981;
traces_original_cutted = traces_original(1:cutindex,:);
traces_bgfitcorr_cutted = traces_bgfitcorr(1:cutindex,:);
clear movie;
save(fullfile(data_path,'longterm_traces.mat'))


% load('I:\0_Code\Figure\traces.mat')
%% fit
% [traces_corrected, fitted_curves] = fit_exp2(traces);
% 计算数据的总长度，并将其分为7段（假设数据长度可被7整除）
% partlength = size(traces,1)/7;
% 计算用于填充的长度（数据长度的5%），用于避免滤波时的边界效应
traces = traces_bgfitcorr_cutted ;
padlength = round(0.1*size(traces,1));
padded_traces = [repmat(traces(1, :), padlength, 1); traces; repmat(traces(end, :), padlength, 1)];
filtered_traces = zeros(size(padded_traces));
Fs = 500; % 采样频率 (Hz)，请根据实际情况调整
FcL = 50;  % 截止频率 (Hz)
% 设置高通滤波器的截止频率。这里的`1/t(end)`是数据的周期，因此高频截止频率为每个周期的反向值转换为每分钟。
t = size(traces,1)/Fs;
FcH = 1/t(end)*60;%截止频率 (Hz)

% 对每一列进行高通滤波校正漂白
for col = 1:size(traces, 2)
    fprintf('Processing roi %d\n',col)
    filtered_traces(:, col) = highpass(padded_traces(:, col), FcH, Fs);
end
traces_dFF = filtered_traces./mean(traces_bgfitcorr_cutted);
% %对每一列进行低通滤波降
% for col = 1:size(traces, 2)
%     fprintf('Processing roi %d\n',col)
%     denoised_traces(:, col) = lowpass(traces_dFF(:, col), FcL, Fs);
% end
% 去除填充部分，得到最终的校正后的信号
traces_corrected = traces_dFF(padlength+1:end-padlength, :);
% plot
fig = figure();
set(fig,'Position',get(0,'Screensize'));
fit_axe = subplot(2,1,1);
fited_axe = subplot(2,1,2);
t = (1:size(traces,1)).*0.002/60;
% plot
for i = 1: size(traces_corrected,2)
    plot(t,traces(:,i),'Parent',fit_axe);
    hold(fit_axe, 'on');
    plot(t,traces_corrected(:,i),'Parent',fited_axe);
    hold(fited_axe, 'on');
end
hold off;

% note
title(fit_axe, 'Original and Fitted Curves');
title(fited_axe, 'Corrected Traces');
legend(fit_axe, 'Original Trace', 'Fitted Curve');
save(fullfile(data_path,'2_fitted_traces.mat'),'traces_corrected','traces','Fs','FcL');
saveas(gcf,fullfile(data_path,'2_fitted_traces.png'));

fprintf('Finished expotential fit\n')
 %% Peak finding(manual)
% 
parts = 4;
MinPeakProminence_factor = 0.5; % 定义最小峰值显著性因子
nrois = size(traces,2);
peaks_polarity = zeros(parts,nrois);
peaks_threshold = zeros(parts,nrois);
peaks_amplitude = cell(parts, nrois);
peaks_index = cell(parts, nrois);
peaks_sensitivity = cell(parts, nrois);

% calculate part index
step_frame = round(length(traces_corrected)/parts) ;
partindex = {};
start_frame = 0;
end_frame = 0;
for i = 1:parts
    start_frame = 1 + end_frame;
    end_frame = round(min(start_frame + step_frame ,length(traces_corrected)));
    partindex{i} = start_frame:end_frame;
end

% peakfinding

for i = 1:parts
    [current_peaks_polarity, current_peaks_threshold, current_peaks_index, current_peaks_amplitude, current_peaks_sensitivity] = ...
        peak_finding(traces_corrected(partindex{i},:),MinPeakProminence_factor)
    peaks_polarity(i,:) = current_peaks_polarity;
    peaks_threshold(i,:) = current_peaks_threshold;
    for j = 1:nrois
        peaks_amplitude{i,j} = current_peaks_amplitude{j};
        peaks_index{i,j} = current_peaks_index{j};
        peaks_sensitivity{i,j} = current_peaks_sensitivity{j};
    end
end
save(fullfile(data_path,'partpeakfinding.mat'), 'traces_corrected', 'peaks_polarity', 'peaks_threshold', 'peaks_index', 'peaks_amplitude', 'peaks_sensitivity');
peaks_index_intergrated = cell(1,nrois);
for i = 1:parts 
    for j = 1:nrois
    peaks_index_intergrated{j} = [peaks_index_intergrated{j} ;peaks_index{i,j}+(i-1)*(step_frame+1)];
    end
end
peaks_index=peaks_index_intergrated ;
pffolder = sprintf('%d partial peakfinding at manual',parts);
mkdir (fullfile(data_path,pffolder))
save(fullfile(data_path,pffolder,'partpeakfinding.mat'), 'traces_corrected', 'peaks_polarity', 'peaks_threshold', 'peaks_index', 'peaks_index_intergrated','peaks_amplitude', 'peaks_sensitivity');


%% peak finding (auto)

parts = 1;

nrois = size(traces,2);
peaks_polarity = zeros(parts,nrois);
peaks_threshold = zeros(parts,nrois);
peaks_amplitude = cell(parts, nrois);
peaks_index = cell(parts, nrois);
peaks_sensitivity = cell(parts, nrois);

% calculate part index
step_frame = round(length(traces_corrected)/parts) ;
partindex = {};
start_frame = 0;
end_frame = 0;
for i = 1:parts
    start_frame = 1 + end_frame;
    end_frame = round(min(start_frame + step_frame -1 ,length(traces_corrected)));
    partindex{i} = start_frame:end_frame;
end
MinPeakProminence_factor = 0.5; % 定义最小峰值显著性因子
threshold_factor = 5;
outlierPercentage = 0.10;
pffolder = sprintf('%d partial peakfinding at threshold = %.1f,MinPeakProminence_factor = %.2f, outlierPercentage = %.2f',parts,threshold_factor,MinPeakProminence_factor,outlierPercentage);
mkdir (fullfile(data_path,pffolder))
%peakfinding
for p = 1:parts
    fprintf('Processing part %d',p)
    for i = 1:nrois
        current_trace = traces_corrected(partindex{p},i);

        % 确定峰值极性
        peaks_polarity(p,i) = -1;
        polarity = 'Negative';

        % set threshold
        meanData = mean(current_trace);
        stdData = std(current_trace);
        zScores = abs((current_trace - meanData) / stdData);
        
        % 计算离群值的数量

        numOutliers = ceil(length(current_trace) * outlierPercentage);
        [~, sortedIndices] = sort(zScores, 'descend');
        inlierIndices = sortedIndices(numOutliers+1:end);

        baseline = current_trace(inlierIndices);
        avg_baseline = mean(baseline);

        % Calculate std from the baseline
        noise = std(baseline);

        % peaks_threshold(p,i) = meanData +threshold_factor*stdData;
        peaks_threshold(p,i) = avg_baseline +threshold_factor*noise;


        plot_trace = current_trace * peaks_polarity(p,i) ;

        % 寻找峰值
        MinPeakProminence = (max(current_trace)-min(current_trace)) * MinPeakProminence_factor;
        [peak_y, peak_x] = findpeaks(plot_trace, 'MinPeakProminence', ...
            MinPeakProminence, 'MinPeakHeight', peaks_threshold(p,i));

        peaks_amplitude{p,i}  = peak_y;
        peaks_index{p,i} = peak_x+start_frame;
        peaks_sensitivity{p,i} = current_trace(peak_x)-mean(current_trace);
        peaks_amplitude{p,i} = peaks_amplitude{p,i} * peaks_polarity(p,i);


        figure()
        plot(plot_trace); hold on;
        plot(ones(1,size(current_trace,1)).*peaks_threshold(p,i),'LineWidth',2);
        if numel(peak_x) ~= 0
        plot(peak_x, (peaks_amplitude{p,i}* peaks_polarity(p,i)),'v');
        end
        title(sprintf('peakfinding of p = %d, roi %d.png',p,i));
        hold off;
        saveas(gcf,fullfile(data_path,pffolder,sprintf('peakfinding of roi %d, p = %d.png',i,p)));
        saveas(gcf,fullfile(data_path,pffolder,sprintf('peakfinding of roi %d, p = %d.fig',i,p)));
        close(gcf);
        if numel(peak_x) == 0
            peaks_polarity(p,i) = 0;
        end
        fprintf('%d peaks found in ROI %d, Polarity : %s\n',numel(peak_x),i,polarity)
    end
    fprintf('All peaks found.\n')
end
save(fullfile(data_path,pffolder,'peakfinding.mat'), 'MinPeakProminence_factor','parts','threshold_factor', ...
    'peaks_polarity', 'peaks_threshold', 'peaks_index', 'peaks_amplitude', 'peaks_sensitivity');



% for i = 1:nrois
%     step_frame = round(length(traces_corrected))/parts ;
% 
%     start_frame = 1;
%     end_frame = 0;
%     for p = 1:parts
% 
% 
%         if p ==1
%             start_frame = 1;
%         else
%             start_frame = end_frame+1 ;
%         end
% 
%         if start_frame + step_frame < length(t)
%             end_frame = start_frame + step_frame;
%         else
%             end_frame = length(t);
%         end
% 
%         current_trace = traces_corrected(start_frame:end_frame,i);
% 
%         % 确定峰值极性
%         % abs(min(current_trace) - mean(current_trace)) < max(abs(current_trace) - mean(current_trace))
% 
%             peaks_polarity(p,i) = -1;
%             polarity = 'Negative';
% 
%         % set threshold
%         peaks_threshold(p,i) =  mean(current_trace)+3*std(current_trace);
%         plot_trace = current_trace * peaks_polarity(p,i) ;
% 
%         % 寻找峰值
%         MinPeakProminence = (max(current_trace)-min(current_trace)) * MinPeakProminence_factor;
%         [peak_y, peak_x] = findpeaks(plot_trace, 'MinPeakProminence', ...
%             MinPeakProminence, 'MinPeakHeight', peaks_threshold(p,i));
% 
%         peaks_amplitude{p,i}  = peak_y;
%         % plot peak
%         % figure()
%         % plot(plot_trace); hold on;
%         % plot(ones(1,size(current_trace,1)).*peaks_threshold(p,i),'LineWidth',2);
%         % plot(peak_x, (peaks_amplitude{p,i}),'v','MarkerFaceColor',colors(i,:));
%         % hold off;
% 
% 
%         % peaks_threshold(i) = (peaks_threshold(i) + peaks_polarity(i) -1)/peaks_polarity(i);
%         peaks_threshold(p,i) = peaks_threshold(p,i)*peaks_polarity(p,i);
%         peaks_index{p,i} = peak_x+start_frame;
%         % current_trace = current_trace;
%         peaks_sensitivity{p,i} = current_trace(peak_x)-mean(current_trace);
% 
% 
% 
%         % peaks_amplitude{i} = (peaks_amplitude{i} + peaks_polarity(i) -1)/peaks_polarity(i);
%         peaks_amplitude{p,i} = peaks_amplitude{p,i} * peaks_polarity(p,i);
% 
%         if numel(peak_x) == 0
%             peaks_polarity(p,i) = 0;
%         end
%         fprintf('%d peaks found in ROI %d, Polarity : %s\n',numel(peak_x),i,polarity)
%     end
% end
% fprintf('All peaks found.\n')
% 
% peaks_polarity =  peaks_polarity(1,1:nrois);
% peaks_amplitude_all = cell(1, nrois);
% peaks_index_all = cell(1, nrois);
% peaks_sensitivity_all = cell(1, nrois);
% for i = 1:nrois
%     for p = 1:parts
%         peaks_amplitude_all{i} = cat(1,peaks_amplitude_all{i},peaks_amplitude{p,i});
%         peaks_index_all{i} = cat(1,peaks_index_all{i},peaks_index{p,i});
%         peaks_sensitivity_all{i} = cat(1,peaks_sensitivity_all{i},peaks_sensitivity{p,i});
% 
%     end
% end
% 
% peaks_amplitude = peaks_amplitude_all;
% peaks_index = peaks_index_all;
% peaks_sensitivity = peaks_sensitivity_all;
% 
% 
% snr_threshold = 3;  % 定义 SNR 阈值
% baselines = mean(traces_corrected);
% % 遍历每个 ROI 和部分，筛选 SNR > 3 的峰值
% for i = 1:nrois
%         % 获取当前部分的峰值
%         logical_index = false(size(traces_corrected, 1), 1);
%         indexs = peaks_index{i};
%         for j = 1:length(indexs)
%             if isempty(j)
%                 continue
%             else
%                 logical_index(round(indexs(j))) = true;
%             end
% 
%         end
%         current_amplitude = traces_corrected(logical_index,i);
%         current_sensitivity = peaks_sensitivity{i};
% 
%         % 计算背景噪声标准差（标准差可从原信号的残差计算）
%         noise_std = std(traces_corrected(:, i)); % 背景噪声
%         if isempty(current_amplitude)
%             continue;
%         end
% 
%         % 计算 SNR
%         snr_values = (current_amplitude-baselines(i)) / noise_std;
% 
%         % 筛选 SNR >= 3 的峰值
%         valid_indices = abs(snr_values) >= snr_threshold;
%         peaks_amplitude{i} = current_amplitude(valid_indices);
%         peaks_index{ i} = peaks_index{i}(valid_indices);
%         peaks_sensitivity{i} = current_sensitivity(valid_indices);
% 
% end
% 
% 
% 
% 
% %  find peaks
% % also peak_finding(traces,MinPeakProminence_factor) for diy factor
% % [peaks_polarity, peaks_threshold, peaks_index, peaks_amplitude, peaks_sensitivity] = ...
% %                                         peak_finding(traces_corrected,0.3);
% traces_corr_flipped = traces_corrected.*peaks_polarity + 1 - peaks_polarity;
% 
% % plot trace
% fig = figure();
% set(fig,'Position',get(0,'Screensize'));
% offset_array = offset_plot(traces_corr_flipped,t);
% xlim tight
% sgtitle('Peak finding');
% % plot label
% for i = 1:nrois
%     % for p = 1:parts
%     % % plot threshold
%     % plot(t,ones(1,numel(t)/parts).*(peaks_threshold(p,i)*peaks_polarity(i)+ 1 -peaks_polarity(i)) ...
%     %     + offset_array(i),'LineWidth',2); hold on;
%     % end
%     % plot peak
%     dt = t(2)-t(1);
%     plot(peaks_index{i}.*dt, (peaks_amplitude{i}.*peaks_polarity(i)+ 1 -peaks_polarity(i)) ...
%         + offset_array(i),'v'); hold on;
% end
% 
% 
% 
% % fig_filename = fullfile('3_peak_finding.fig');
% % png_filename = fullfile('3_peak_finding.png');
% % peak_filename = fullfile( '3_peak_finding.mat');
% % 
% % saveas(gcf, fig_filename, 'fig');
% % saveas(gcf, png_filename, 'png');

%% bin spike heatmap
peaks_index_intergrated = cell(1,nrois);
for i = 1:parts 
    for j = 1:nrois
    peaks_index_intergrated{j} = [peaks_index_intergrated{j} ;peaks_index{i,j}+(i-1)*(step_frame+1)];
    end
end

binsize = 500;
fullpeaksarray = zeros(size(traces_corrected));
for i = 1:size(peaks_index_intergrated,2)
    currentpeaks = peaks_index_intergrated{i};
    for j = 1:size(currentpeaks,1)
        fullpeaksarray(currentpeaks(j),i) = 1;
    end
end
binnum = floor(size(traces_corrected,1)/binsize);
binarray = fullpeaksarray(1:binnum*binsize,:);
binarray = reshape(binarray,binsize,binnum,size(traces_corrected,2));
binarray = sum(binarray,1);
binarray = squeeze(binarray)';

% 生成示例矩阵（包含0和正整数）

% 显示矩阵
imagesc(binarray);

% 自定义颜色映射
cmap = [
    0.95 0.95 0.95;    % 
    0.9 0.75 0.75;    % 
    0.9 0.2 0.2;    % 
    ];

% 使用插值的方式生成更平滑的颜色梯度
cmap_interp = interp1([0, 1, max(binarray(:))], cmap, linspace(0, max(binarray(:)), 256));

% 应用自定义的颜色映射
colormap(cmap_interp);

% 设置颜色轴范围，将0和正整数的颜色差异显现出来
clim([0 max(binarray(:))]);
set(gca, 'TickLength', [0 0]);
title(sprintf('Heat map of spikes at bin %d',binsize))
% 添加颜色条
colorbar;
%%
figure
SNR_traces = zeros(size(traces_corrected));
for i = 1: nrois
    SNR_traces(:,i) = calculate_SNR(traces_corrected(:,i),peaks_index{i},AP_window_width);
end
offset_plot(SNR_traces,t);

