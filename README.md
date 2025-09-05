The code in this repository include all data preprocessing, quality control and figure plots for the orographic enhancement paper (2025).
Due to the license of original data sources (ERA5-Land, MODIS, and Soils data) and the storage limit of our local resources, these raw data need to be downloaded (~20TB) prior to use the scripts in this repository.
The core of this paper lies in the calculation of rainfall uncertainty. The execution steps are as follows:

    1. Use Prepare_inputs_summary.m script to prepare inputs including: interpolation of raw datasets, mapping to basins, quality control.
    2. Use Run_world_900m_IRC.sh to conduct a rainfall uncertainty estimation process (IRC, or inverse rainfall correction). 
        Since calculations in this work are mostly done using matrix fashion, Matlab was used for faster performance. Therefore, the matlab scripts in this repository need to be adjusted to user's working path.
        The IRC procedure includes the following scripts: Pre_rising_WORLD_30day.m; rainupv2_WORLD_30day.m; Slow_recession_v2_WORLD_30day.m with each model time step requiring an execution of particle tracking using Backtrack_IRC_WORLD_30days.m and transition_WORLD_30day.m.
        An executable of the hydrological model is provided in this repository. Detailed user guide/instructions regarding the DCHM model will be provided upon request since we are working on licensing it. 
    3. Use plot_IRC_results.m script to generate figures presented in the paper.
    4. Use AI_transfer_learning.m script to construct Random Forest Regressor for Cross-regional model validataion and generate corresponding figures. 

Note: outputs from the IRC process are very large due to the use of Lagrangian tracking at high resolution (5 minutes scale), the total raw data outputs is >80TB in binary format. These raw outputs can be downloaded upon request.
