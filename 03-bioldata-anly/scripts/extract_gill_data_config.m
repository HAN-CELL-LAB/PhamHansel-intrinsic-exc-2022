%% Stimulus info 
stim_info = struct; 
stim_info.pulse_on_time         = 100; % ms 
stim_info.pulse_off_time        = 600; % ms 
stim_info.hyperpol_on_time      = 1000; % ms
stim_info.hyperpol_off_time     = 1100; % ms
stim_info.hyperpol_Iinj_amp     = 100; % pA

%% Defining parameters 
spike_params = struct; 
spike_params.peak_prom          = 40; 
spike_params.min_dist           = 1; 
spike_params.min_height         = -10; 
spike_params.min_width          = 0.05;  
spike_params.analyze_in_stim    = true; 

AP_params = struct; 
AP_params.max_tpre              = 1; 
AP_params.min_dV_dt             = 8;  
AP_params.scale_min_dV_dt       = false; 
% AP_params.min_dV_dt             = 0.1;  
% AP_params.scale_min_dV_dt       = true; 
AP_params.min_pts_cross_thres   = 2; 
AP_params.max_tpost             = 10; 
AP_params.rise_start_APfactor   = 0.1;
AP_params.rise_stop_APfactor    = 0.9;
AP_params.rise_interp_factor    = 100;
AP_params.fall_start_APfactor   = 0.9;
AP_params.fall_stop_APfactor    = 0.1; 
AP_params.fall_interp_factor    = 100;

pass_params = struct; 
pass_params.Vbase_calc_time_window  = 20; 
pass_params.Rin_calc_time_window    = 20; 

%% Save all in one
configs = struct; 
configs.stim_info = stim_info; 
configs.spike_params = spike_params;
configs.AP_params = AP_params; 
configs.pass_params = pass_params;
