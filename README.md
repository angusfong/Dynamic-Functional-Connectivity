## Code supporting our paper "Dynamic functional connectivity during task performance and rest predicts individual differences in attention across studies", NeuroImage (https://doi.org/10.1016/j.neuroimage.2018.11.057) in December 2018.
### Authors: Angus Ho Ching Fong, Kwangsun Yoo, Monica D. Rosenberg, Sheng Zhang, Chiang-Shan R. Li, Dustin Scheinost, R. Todd Constable, Marvin M. Chun

Abstract:
Dynamic functional connectivity (DFC) aims to maximize resolvable information from functional brain scans by considering temporal changes in network structure. Recent work has demonstrated that static, i.e. time-invariant resting-state and task-based FC predicts individual differences in behavior, including attention. Here, we show that DFC predicts attention performance across individuals. Sliding-window FC matrices were generated from fMRI data collected during rest and attention task performance by calculating Pearson's r between every pair of nodes of a whole-brain atlas within overlapping 10–60s time segments. Next, variance in r values across windows was taken to quantify temporal variability in the strength of each connection, resulting in a DFC connectome for each individual. In a leave-one-subject-out-cross-validation approach, partial-least-square-regression (PLSR) models were then trained to predict attention task performance from DFC matrices. Predicted and observed attention scores were significantly correlated, indicating successful out-of-sample predictions across rest and task conditions. Combining DFC and static FC features numerically improves predictions over either model alone, but the improvement was not statistically significant. Moreover, dynamic and combined models generalized to two independent data sets (participants performing the Attention Network Task and the stop-signal task). Edges with significant PLSR coefficients concentrated in visual, motor, and executive-control brain networks; moreover, most of these coefficients were negative. Thus, better attention may rely on more stable, i.e. less variable, information flow between brain regions.

Files for main results: <br>
combModel_nested.R: runs internal validation, calls comb_predict.R <br>
comb_predict.R: runs one iteration of LOOCV<br>
tune_w.R: optimizes parameter w<br>
