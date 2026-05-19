%第二部分时域处理
%开始时域处理了，最常用的，发sci足够了


%数据导入并分段
clear all; clc; close all
data_path = 'path_to_your_data';%预处理完数据路径
cd(data_path)
files = dir('*.set');
fn = {files.name};
Cond = {  '1'  '2'  }; %% condition name  mark名称
for i = 1:length(fn)
     setname = fn{i}; %% filename of set file
     EEG = pop_loadset('filename',setname,'filepath',data_path); 
     EEG = eeg_checkset( EEG );
    for j = 1:length(Cond)%对各Mark重新分段
        EEG_new = pop_epoch( EEG, Cond(j), [-0.2  1], 'newname', 'Merged datasets pruned with ICA', 'epochinfo', 'yes'); %% epoch by conditions, input to EEG_new
        EEG_new = eeg_checkset( EEG_new );
        EEG_new = pop_rmbase( EEG_new, [-200     0]); %% baseline correction for EEG_new
        EEG_new = eeg_checkset( EEG_new );
        EEG_avg(i,j,:,:) = squeeze(mean(EEG_new.data,3));  %% average across trials for EEG_new, EEG_avg dimension: subj*cond*channel*time
    end 
end
EEG.times = EEG.times(:, 401:1000);%更新一下时间数据方便后面画图
%%  单因素找显著时间窗
data_test = squeeze(EEG_avg(:,:,28,:)); %% 选电极找显著, data_test: subj*cond*time
%对于每个时间点
for i = 1:size(data_test,3)
    %提取所有被试 所有条件 第i个时间点的数据
    %data_anova sub * cond
    data_anova = squeeze(data_test(:,:,i)); %% select the data at time point i
    %在第i个时间点下进行重复测量方差分析
    %注意 在使用anova_rm函数时 要保证函数在当前路径下 或者set path中
    %此函数只能做单因素重复测量方差分析
    %要求输入的数据组织形式是 被试为行 条件为列，off关闭弹窗
    %输出变量p中的第一个值是条件的主效应
    %table是方差分析的表
    [p, table] = anova_rm(data_anova,'off');  %% perform repeated measures ANOVA
    %汇总每次统计下的p值
    P_anova(i) = p(1); %% save the data from ANOVA
end

mean_data = squeeze(mean(data_test,1)); %% dimension: cond*time
figure; 
subplot(211);plot(EEG.times, mean_data,'linewidth', 1.5); %% waveform for different condition 
set(gca,'YDir','reverse');
axis([-200 1000 -35 25]);
subplot(212);plot(EEG.times,P_anova); axis([-200 1000 0 0.05]); %% plot the p values from ANOVA
%%
%% point-by-point repeated measures of ANOVA across channels
%找到感兴趣时间范围内的采样点的位置信息
test_idx = find((EEG.times>=80)&(EEG.times<=120)); %% define the intervals
%sub * con * ch
data_test = squeeze(mean(EEG_avg(:,:,:,test_idx),4)); %% select the data in [197 217]ms, subj*cond*channel
%对于每个通道
for i = 1:size(data_test,3)
    data_anova = squeeze(data_test(:,:,i)); %% select the data at channel i
    %在第i个时间点下进行重复测量方差分析
    %注意 在使用anova_rm函数时 要保证函数在当前路径下 或者set path中
    %此函数只能做单因素重复测量方差分析
    %要求输入的数据组织形式是 被试为行 条件为列，off关闭弹窗
    %输出变量p中的第一个值是条件的主效应

    [p, table] = anova_rm(data_anova,'off');  %% perform repeated measures ANOVA
    %汇总每次统计下的p值
    P_anova111(i) = p(1); %% save tANOVA
    
end

figure;
for i = 1:2
    subplot(1,3,i); 
    topoplot(squeeze(mean(data_test(:,i,:),1)),EEG.chanlocs,'maplimits',[-20 20]); 
end
subplot(1,3,3); topoplot( P_anova111,EEG.chanlocs,'maplimits',[0 0.05]);


%% 提取脑电统计的数据，选择通道，时间在上环节提取，这些数据放spss里
clear all; clc; close all
data_path = 'path_to_your_data';
cd(data_path)
files = dir('*.set');
fn = {files.name};
%Cond = {  'S 21'  'S 22'  'S 23'  'S 24'  };
Cond = {  '1'  '2'  '3'  '4' }; %% condition name
for i = 1:length(fn)
     setname = fn{i}; %% filename of set file
     EEG = pop_loadset('filename',setname,'filepath',data_path); 
     EEG = eeg_checkset( EEG );
    for j = 1:length(Cond)
        EEG_new = pop_epoch( EEG, Cond(j), [-0.2  1], 'newname', 'Merged datasets pruned with ICA', 'epochinfo', 'yes'); %% epoch by conditions, input to EEG_new
        EEG_new = eeg_checkset( EEG_new );
        EEG_new = pop_rmbase( EEG_new, [-200     0]); %% baseline correction for EEG_new
        EEG_new = eeg_checkset( EEG_new );
        EEG_avg(i,j,:,:) = squeeze(mean(EEG_new.data,3));  %% average across trials for EEG_new, EEG_avg dimension: subj*cond*channel*time
    end 
end
EEG.times = EEG.times(:, 401:1000);%更新一下时间数据方便后面画图

data_test = squeeze(EEG_avg(:,:,46,:)); %% 选电极找显著, data_test: subj*cond*time
mean_data = squeeze(mean(data_test,1)); %%可提取为画图数据 可以直接MATLAB画图，也可以复制黏贴出来放origion画图  dimension: cond*time



test_idx = find((EEG.times>=80)&(EEG.times<120)); %
data_test = squeeze(mean(EEG_avg(:,:,:,test_idx),4));
data_djtongji= squeeze(data_test(:,:,29));%选择想
% define the intervals
%sub * con * ch提取电极

%% 画地形图
N2_interval=find((EEG.times>=80)&(EEG.times<=120)); %% N2 interval
P2_interval=find((EEG.times>=400)&(EEG.times<=800)); %% P2 interval
%EEG_avg 4D sub * con * ch * t
%提取所有被试 所有条件 所有通道 感兴趣时间范围内的数据
%沿着时间做平均 沿着被试做平均
%-> con * ch
N2_amplitude=squeeze(mean(mean(EEG_avg(:,:,:,N2_interval),4),1)); %% N2 amplitude for each subject, condition, and channels
P2_amplitude=squeeze(mean(mean(EEG_avg(:,:,:,P2_interval),4),1)); %% P2 amplitude for each subject, condition, and channels
%将两个条件的数据进行纵向拼接
all_amp = [N2_amplitude;P2_amplitude];
%将数据取绝对值
abs_all_amp = abs(all_amp);
%找最大值,把并向下取整
max_val = floor(0.7*max(abs_all_amp(:)));
%生成一个空画布
figure; %% divide the panel into 4 rows and 2 colums
%对于每个条件
for i = 1:4
    N2_data = N2_amplitude(i,:); %% average across subjects
    %绘制2行四列的第i张
    subplot(2,4,i); 
    topoplot(N2_data,EEG.chanlocs,'maplimits',[-1*max_val max_val]); title(['L',num2str(i),' N2 Amplitude'],'fontsize',16); %% plot N2 scalp map (group-level)
    P2_data = P2_amplitude(i,:); %% average across subjets
    %绘制两行四列的i+4张
    subplot(2,4,i+4); 
    topoplot(P2_data,EEG.chanlocs,'maplimits',[-1*max_val max_val]); title(['L',num2str(i),' P2 Amplitude'],'fontsize',16); %% plot P2 scamp map (group-level)
end
