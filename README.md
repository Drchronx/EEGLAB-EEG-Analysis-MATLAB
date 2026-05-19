# EEGLAB-EEG-Analysis-MATLAB

基于 MATLAB 和 EEGLAB 的脑电数据分析脚本集合，覆盖 EEG 预处理、ERP、ERSP/时频分析、静息态与任务态功能连接、ROI 统计和全脑连接扫描。项目适用于心理学、认知神经科学、人机交互、神经营销等实验场景中的 EEG 数据处理与探索性分析。

This repository provides MATLAB scripts for EEG analysis based on EEGLAB. It includes workflows for preprocessing, ERP analysis, ERSP/time-frequency analysis, resting-state and task-based functional connectivity, ROI-level statistics, and whole-brain connectivity scanning.

## Features

- EEG preprocessing with EEGLAB, including import, resampling, channel selection, re-referencing, filtering, epoching, and baseline correction
- ERP analysis for single-factor, within-subject, and mixed experimental designs
- ERSP/time-frequency analysis with STFT-based workflows
- Resting-state EEG segmentation and analysis
- ROI-level and whole-brain functional connectivity analysis
- Mass-univariate statistical testing and FDR correction
- Export support for visualization tools such as BrainNet Viewer
- Scripts for extracting data for SPSS/Excel-based follow-up statistics

## Repository Structure

```text
eeg_pingyu/   EEG preprocessing and frequency-domain analysis examples
erp/          ERP preprocessing and statistical analysis scripts
ersp/         ERSP and time-frequency analysis scripts
gnlj/         Functional connectivity analysis and BrainNet Viewer export scripts
```

## Requirements

- MATLAB
- EEGLAB
- EEGLAB-compatible data import plugins as required by your raw data format, such as Curry/Neuroscan import support
- BrainNet Viewer, only for scripts that export or render `.node` and `.edge` files

The scripts were written as research-oriented MATLAB workflows. Before running them, update paths, event markers, channel indices, time windows, frequency bands, and experimental design settings according to your own dataset.

## Usage

1. Install MATLAB and EEGLAB.
2. Add EEGLAB and this repository to the MATLAB path.
3. Replace example values such as `group_dir`, `data_path`, `save_path`, event codes, channel indices, ROI definitions, and condition names.
4. Run the relevant script for your analysis target:
   - ERP analysis: scripts under `erp/`
   - ERSP/time-frequency analysis: scripts under `ersp/`
   - Functional connectivity analysis: scripts under `gnlj/`
   - EEG preprocessing/frequency analysis examples: scripts under `eeg_pingyu/`

## Data and Privacy

This repository is intended to contain source code only. Do not upload raw EEG data, processed participant data, identifying information, local experiment folders, or unpublished project materials. The `.gitignore` file excludes common EEG data formats, MATLAB outputs, and intermediate processing folders by default.

Some scripts contain placeholder paths such as `path_to_your_data`. Replace these placeholders with paths on your own machine before running the scripts.

## Suggested GitHub About

MATLAB scripts for EEG analysis based on EEGLAB, covering preprocessing, ERP, ERSP/time-frequency analysis, resting-state and task-based functional connectivity, ROI-level statistics, and whole-brain connectivity analysis.

## License

This project is released under the PolyForm Noncommercial License 1.0.0. See [LICENSE](LICENSE) for details.

Noncommercial use is permitted. Commercial use is not permitted without separate written permission from the repository owner. See [COMMERCIAL.md](COMMERCIAL.md).

Because this license restricts commercial use, this repository should be described as source-available for noncommercial use rather than OSI-approved open source.
