The code in this repository include all data preprocessing, quality control and figure plots for the orographic enhancement paper (2025).
Due to the license of original data sources (ERA5-Land, MODIS, and Soils data) and the storage limit of our local resources, these raw data need to be downloaded prior to use the scripts in this repository.
The core of this paper lies in the calculation of rainfall uncertainty. The execution steps are as follows:

    Use Prepare_inputs_summary.m script to prepare inputs including: interpolation of raw datasets, mapping to basins, quality control.
    Use Run_world_900m_IRC.sh to conduct a rainfall uncertainty estimation process (IRC, or inverse rainfall correction). 
    Since the calculations are mainly done using matrix, Matlab was used for faster performance. Therefore, the matlab scripts in this repository need to be adjusted to user's working path.
    Use plot_IRC_results.m script to generate figures presented in the paper.
    Use AI_transfer_learning.m script to construct feedforward neural networks for transfer learning and generate corresponding figures. 

Note: outputs from the IRC process are very large due to the use of Lagrangian tracking at high resolution (5 minutes scale), the total raw data outputs is >80TB in binary format. These raw outputs can be downloaded upon requests.
