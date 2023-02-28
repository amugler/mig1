Matlab code:

To plot Fig 2D, 3A:
- Run fig23_process1.m: fits activator model to control data (fig234_data_ctrl.txt)
- Run fig23_process2.m: fits repressor model to control data (fig234_data_ctrl.txt)
- Run fig23_plot: plots figure panels

To plot Fig 4E,F:
- Run fig4_process1.m: finds division time from control data (fig234_data_ctrl.txt)
- Run fig4_process2.m: finds mean size ratio from wild-type data (fig4_data_ratio.txt)
- Run fig4_process3.m: fits activator model with division to control data (fig234_data_ctrl.txt)
- Run fig4_plot: plots figure panels

To compute statistics for Fig 5B:
- Run fig5_process: runs ANOVA on data (fig5_data_ctrl.txt, fig5_data_gof.txt, fig5_data_lof.txt)

To plot supplementary feedback models:
- Run auto_ctrl.m: fits auto-regulation model to control data (fig5_data_ctrl.txt)
- Run auto_lof.m: fits auto-regulation model to bar-1(ga80) data (fig5_data_lof.txt)
- Run auto_gof.m: fits auto-regulation model to ΔN-BAR-1 data (fig5_data_gof.txt)
- Run auto_noise.m: computes CV for auto-regulation model
- Run both_ctrl.m: fits dual-feedback model to control data (fig5_data_ctrl.txt)
- Run both_lof.m: fits dual-feedback model to bar-1(ga80) data (fig5_data_lof.txt)
- Run both_gof.m: fits dual-feedback model to ΔN-BAR-1 data (fig5_data_gof.txt)
- Run both_noise.m: computes CV for dual-feedback model
- Run feedback_noise.m: plots figure panels
