
WIP for standardized pre-processing pipline for data collected at NatMEG.

# Scripts
### opm_analysis_main.m
Main script to run analysis. Specify paths to data and parameters (e.g. filter cutoff, pre- and post-stim time window) here.

### ica_MEG.m
Function to run ICA on OPM or SQUID data (depending on params.ch) with option to manually select components to remove from data.

### read_osMEG.m
Read in On Scalp MEG data from OPM sensors, including filtering, downsampling, epoching and identification of bad channels (using `opm_badchannels.m`).

### opm_badchannels.m
Function to detect bad channels in OPM data. Return `badchs` with:
- `badchs_flat `: Flat channels, i.e. channels with no signal.
- `badchs_std`: Channels with excessive standard deviation, likely contaminated by high levels of noise or artifacts.
- `badchs_neighbors`: Channels that have poor correlation with their spatial neighbors.
- `badchs_outlier`: Channels that show outlier behavior in their frequency spectrum.

### read_cvMEG.m
- Read MEG data from conventional SQUID sensors, including filtering and epoching.

### timelock.m
Create time-locked epochs for OPM or SQUID data. Baseline/demean data and save butterfly plot.

### insidePointcloud.m
Not in use

# To do
- Read data from BIDS, loop over subjects from CSV/list
- Read parameters from separate config file
- Separate OPM and SQUID analysis
- Don't specify params.manual_ica more than once (now in params, function and when calling function?)
    


