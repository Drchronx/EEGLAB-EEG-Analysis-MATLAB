%第三部分时频分析，可以结合时域结果一起补充撰写，可解释数据多一维度，有时候有新发现新惊喜。
clear all; clc;
data_path = 'path_to_your_data';%定义预处理完的数据路径，我这就不该了按自己数据把路径更改
cd(data_path)
files = dir('*.set');
fn = {files.name};
%定义条件信息
Cond={  '1'  '2'  '3'  '4'  };

%% time-frequency analysis for multiple conditions

% 对于所有被试 
for i=1:length(fn)
    setname = fn{i}; %% filename of set file
    EEG = pop_loadset('filename',setname,'filepath',data_path); 
    EEG= eeg_checkset( EEG );
    for j=1:length(Cond)
        %已重新分段的方式提取j条件的数据
        EEG_new = pop_epoch( EEG, Cond(j), [-1  2], 'newname', 'Merged datasets pruned with ICA   epochs epochs', 'epochinfo', 'yes'); 
        EEG_new = eeg_checkset( EEG_new );
        EEG_new = pop_rmbase( EEG_new, [-1000     0]); 
        EEG_new = eeg_checkset( EEG_new );
        %对每个通道
        for nchan=1:size(EEG_new.data,1)
            %提取第nchan通道下的所有时间分段的数据
            %时间 *分段
            x = squeeze(EEG_new.data(nchan,:,:)); 
            %定义时间轴信息
            xtimes=EEG_new.times/1000; 
            %定义在什么时间点加窗
            t=EEG_new.times/1000;
            %定义频率轴
            f=1:1:30; 
            %定义采样率
            Fs = EEG.srate;
            %定义时间窗
            winsize = 0.400; 
            %进行STFT
            %
            [S, P, F, U] = sub_stft(x, xtimes, t, f, Fs, winsize); 
            % 被试 * 条件 * 通道 * 频率 * 时间点
            %进行个体水平叠加平均（分段）
            %P_DATA 5维  被试*通道*频率*时间
            P_data(i,j,nchan,:,:)=squeeze(mean(P,3)); %%P_data (without baseline correciton):  subj*cond*chan*f*time
        end
    end    
end
%%
%save('P_data.mat','P_data')
%load('P_data.mat');
%% baseline correction 
%p_data 5维 被试、条件、通道、频率、时间
% 对于所有被试、条件、通道、频率，都做基线校正
%定义基线范围
t_pre_idx=find((t>=-0.8)&(t<=-0.2));
for i=1:size(P_data,1)
    %对每个被试
    for j=1:size(P_data,2)
        %对于每个通道
        for ii=1:size(P_data,3)
            %对于每个频率点
            for jj=1:size(P_data,4)
                %对于时间序列
                temp_data=squeeze(P_data(i,j,ii,jj,:));
                %进行基线校正
                P_BC(i,j,ii,jj,:)=temp_data-mean(temp_data(t_pre_idx));
            end
        end
    end
end
%%
save('P_BC.mat','P_BC')
%load('P_BC.mat');
%% %%多因素找显著
%% 

% 定义时间 (毫秒)，并转化为 秒 (s) 以匹配你的 ROI [0.5 1.0]
% 你的数据是 -1000ms 到 1998ms，间隔 2ms
%t_ms = -1000:2:1998; 
%t = t_ms / 1000;  % 【关键】转换成秒！也就是 -1.0s 到 1.998s

% 定义频率 f (假设是1-30Hz，需根据你 data_test 的第3维度修改)
% 请务必确认 size(data_test, 3) 是多少，如果是 30，就写 1:30
%f = 1:30;   % <--- 如果你的 f 还没定义，请取消注释并修改它




data_test=squeeze(P_BC(:,:,19,:,:));
for i=1:size(data_test,3)
    %对于每个时间点
     for j=1:size(data_test,4)
        %提取每个被试每个条件ij时评点下的数据
        data_anova=squeeze(data_test(:,:,i,j)); %% select the data at time-frequency point
        %进行重复测量方差分析
        sub = 26; % 被试数量
        F1_con = 2; %F1自变量水平数
        F2_con = 2; %F2自变量水平数
        total_con = F1_con * F2_con ;
        S = repmat([1:sub], 1, total_con).';
         % F1和F2 为自变量的水平
        F1 = sort(repmat([1:F1_con], 1, sub*F2_con)).';
        F2 = repmat(sort(repmat([1:F2_con], 1, sub)), 1, F1_con).';
        %FACTNAMES为变量名称
        FACTNAMES = {'一条件', '二条件'};
        %重复测量方差分析
        tab = rm_anova2(data_anova , S, F1, F2, FACTNAMES);
         %输出p值
        pzzz = tab{3, 6};%第二条件主效应
        pzz = tab{2, 6}%第一条件主效应
        p = tab{4, 6}%交互作用
        %汇总每次统计下的p值
        P_anovazz(i,j) = pzz(1); %% save the data from ANOVA
        P_anovazzz(i,j) = pzzz(1); 
        P_anova(i,j) = p(1); 
     end
end

%figure; imagesc(t,f,P_anova); axis xy;
%% 1. 处理数据：将不显著的值设为 NaN
P_plot = P_anova; 
P_plot(P_anova >= 0.05) = NaN;

% 2. 绘图
figure; 
imagesc(t, f, P_plot); 
axis xy; 
colorbar;

% 3. 调整色标范围 (旧版本使用 caxis)
caxis([0 0.05]); 

% 4. 设置背景颜色为白色（让NaN部分显示为白色）
set(gca, 'color', 'w'); 

% 可选：反转色图，让P值越小（越显著）颜色越深/红
colormap(flipud(jet));
%%
%% 0. 画图没有框
% 假设你关注的是第19通道
chan_idx = 19; 

% 分别计算4个条件的群组平均 (Grand Average)
% 假设 P_BC 的第2维是条件 (1, 2, 3, 4)
GrandAvg_C1 = squeeze(mean(P_BC(:, 1, chan_idx, :, :), 1)); 
GrandAvg_C2 = squeeze(mean(P_BC(:, 2, chan_idx, :, :), 1)); 
GrandAvg_C3 = squeeze(mean(P_BC(:, 3, chan_idx, :, :), 1)); 
GrandAvg_C4 = squeeze(mean(P_BC(:, 4, chan_idx, :, :), 1)); 

%% 1. 计算统一的色标范围 (关键步骤)
% 将所有4个条件的数据拼在一起找最大最小值，确保颜色对比真实有效
all_values = [GrandAvg_C1(:); GrandAvg_C2(:); GrandAvg_C3(:); GrandAvg_C4(:)];

% 方案A：自动范围 (从全局最小到全局最大)
c_lims = [min(all_values) max(all_values)]; 

% 方案B：对称范围 (如果你看的是ERD/ERS或相减后的差波，推荐用这个)
% limit = max(abs(all_values));
% c_lims = [-limit limit]; 

%% 2. 绘图 (采用 2x2 布局)
% 设置画布大小，使其接近正方形以容纳 2x2 的图
figure('Color', 'w', 'Position', [100, 100, 1000, 800]); 

% 为了代码简洁，我们将数据放入一个元胞数组(Cell Array)中进行循环绘图
data_list = {GrandAvg_C1, GrandAvg_C2, GrandAvg_C3, GrandAvg_C4};
title_list = {'Condition 1', 'Condition 2', 'Condition 3', 'Condition 4'};

for i = 1:4
    subplot(2, 2, i); % 创建 2行2列 的子图
    
    % 绘图核心
    imagesc(xtimes, f, data_list{i}); 
    axis xy;                % 翻转Y轴
    caxis(c_lims);          % 【重要】锁定统一色标范围 (新版MATLAB可用 clim(c_lims))
    
    % 装饰与标注
    title([title_list{i} ' (Ch ' num2str(chan_idx) ')'], 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    % 绘制0时刻竖线
    hold on;
    xline(0, '--k', 'LineWidth', 1.5); 
    hold off;
    
    % 色条与配色
    colorbar;               % 显示色条
    colormap(jet);          % 设定色系 (也可以换成 'parula' 或 'turbo')
end

% 这是一个可选步骤：添加一个总标题
sgtitle(['Time-Frequency Representations - Channel ' num2str(chan_idx)]);



%% 0. 准备工作：提取4个条件的数据
%% 0. 准备工作：提取4个条件的数据
% 假设你关注的是第19通道
chan_idx = 19; 

% 分别计算4个条件的群组平均
GrandAvg_C1 = squeeze(mean(P_BC(:, 1, chan_idx, :, :), 1)); 
GrandAvg_C2 = squeeze(mean(P_BC(:, 2, chan_idx, :, :), 1)); 
GrandAvg_C3 = squeeze(mean(P_BC(:, 3, chan_idx, :, :), 1)); 
GrandAvg_C4 = squeeze(mean(P_BC(:, 4, chan_idx, :, :), 1)); 

%% 1. 计算统一的色标范围
all_values = [GrandAvg_C1(:); GrandAvg_C2(:); GrandAvg_C3(:); GrandAvg_C4(:)];
c_lims = [min(all_values) max(all_values)]; 

%% 2. 绘图 (采用 2x2 布局)
figure('Color', 'w', 'Position', [100, 100, 1000, 800]); 

data_list = {GrandAvg_C1, GrandAvg_C2, GrandAvg_C3, GrandAvg_C4};
title_list = {'Condition 1', 'Condition 2', 'Condition 3', 'Condition 4'};

% === 【设置框框参数：700-800ms, 4-7Hz】 ===
% 注意：请确保你的 xtimes 单位是毫秒(ms)
% 如果是秒(s)，请改成 [0.7 0.8]
ROI_time = [0.600 0.800];  
ROI_freq = [4 7];

% 计算矩形的 [x起点, y起点, 宽度, 高度]
rect_pos = [ROI_time(1), ROI_freq(1), ROI_time(2)-ROI_time(1), ROI_freq(2)-ROI_freq(1)];

for i = 1:4
    subplot(2, 2, i); 
    
    % 1. 绘制时频图
    imagesc(xtimes, f, data_list{i}); 
    axis xy;                
    caxis(c_lims);          
    
    hold on; 
    
    % 2. 【绘制实线长方形框】
    rectangle('Position', rect_pos, ...
              'EdgeColor', 'k', ...     % 边框颜色：黑色('k')，看不清可改为白色('w')
              'LineWidth', 2, ...       % 线宽
              'LineStyle', '-');        % 【修改处】：'-' 代表实线
          
    % 3. 绘制0时刻竖线
    xline(0, '--k', 'LineWidth', 1.5); 
    
    hold off;
    
    % 4. 装饰与标注
    title([title_list{i} ' (Ch ' num2str(chan_idx) ')'], 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Time (ms)'); 
    ylabel('Frequency (Hz)');
    
    colorbar;               
    colormap(jet);          
end

sgtitle(['Time-Frequency Representations - Channel ' num2str(chan_idx)]);


%% fdr correction to account for multiple comparisons
%使用eeglab工具箱的fdr函数进行fdr校正
%p_fdr是说原来的p只需要小于等于多少才能通过校正
%p_masked显示该结果是否过校正 1过 0没过
[p_fdr1, p_masked] = fdr(P_ttest, 0.05); %% fdr correction for p values from ttest
figure; imagesc(t,f,P_ttest); axis xy; caxis([0 p_fdr1]); 

% 对方差分析的结果做FDR 校正
[p_fdr2, p_masked] = fdr(P_anova, 0.05);%% fdr correction for p values from ANOVA
figure; imagesc(t,f,P_anova); axis xy; caxis([0 p_fdr2]); 
%figure; imagesc(t,f,p_masked); axis xy; caxis([0 1]); 

% 对方差分析的结果做FDR 校正
[p_fdr2, p_masked] = fdr(P_anova, 0.05);%% fdr correction for p values from ANOVA
figure; imagesc(t,f,P_anova); axis xy; caxis([0 p_fdr2]); 
%figure; imagesc(t,f,p_masked); axis xy; caxis([0 1]); 

%%
%画出几号通道几条件的情况的组平均时频能量图
clc;close all

aim = '4';
for i=1:length(Cond)
   if Cond{i} == aim
      cond_ind = i ;
      break
   end
end

channel_plot=19; 
figure;  
imagesc(t*1000,f,squeeze(mean(P_BC(:,cond_ind,channel_plot,:,:),1))); 
axis xy; colorbar; 
xlabel('Time (ms)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12); 
title('Baseline-corrected TFR','fontsize',15);

%%
%请在一张图上画出四种情况的0.2s-0.5s的频段的组平均地形图
limits= [-0.5 0.5];
ROI_t=[0.60 0.80]; ROI_f=[4 7]; 
ROI_t_idx=find((t>=ROI_t(1))&(t<=ROI_t(2)));
ROI_f_idx=find((f>=ROI_f(1))&(f<=ROI_f(2)));

TFD_huatu=squeeze(mean(mean(mean(P_BC(:,:,:,ROI_f_idx,ROI_t_idx),1),4),5)); 
figure;
for j=1:length(Cond)
    plot_data = TFD_huatu(j,:);
    subplot(2,3,j);
    topoplot(plot_data,EEG.chanlocs,'maplimits',limits);colorbar; 
    title([Cond{j},' TFD\_频段'],'fontsize',10);
   
end


%% 提取数据，选择通道，选择时间，选择频率

data_test=squeeze(P_BC(:,:,20,:,:));
%被试×条件×频率×时间
ROI1_t=[0.60 0.80]; ROI1_f=[4 7];
ROI1_t_idx=find((t>=ROI1_t(1))&(t<=ROI1_t(2)));
ROI1_f_idx=find((f>=ROI1_f(1))&(f<=ROI1_f(2)));
data_test1=squeeze(mean(mean(data_test(:,:,ROI1_f_idx,ROI1_t_idx),3),4))
