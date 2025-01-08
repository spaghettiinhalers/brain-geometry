%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculates accuracy and MSE for the reconstruction of a task fMRI spatial map using eigenmodes
%%%
%%% Original: James Pang, Monash University, 2022
%%%
%%% Adjusted: Jill Bay, MPI Leipzig, 2023
%%% Adjusted: Michelle Zhou, Max Planck Institute for Human Cognitive and Brain Sciences, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Central directory in which all project code and data are stored
project_dir = '/data/pt_02994/';

% Load relevant repository with MATLAB functions
addpath(genpath('functions_matlab'));

% Load surface files
surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = 'midthickness';

[vertices, faces] = read_vtk(sprintf('%sdata/template_surfaces_volumes/%s_%s-%s.vtk', project_dir, surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices;
surface_midthickness.faces = faces;

% Load cortex mask
cortex = dlmread(sprintf('%sdata/template_surfaces_volumes/%s_cortex-%s_mask.txt', project_dir, surface_interest, hemisphere));
cortex_ind = find(cortex);

num_vertices = length(cortex);
hemisphere = 'lh';

disp('Surface files loading complete.');



%% RECONSTRUCT A SINGLE-SUBJECT TASK FMRI SPATIAL MAP

disp('BEGINNING TASK FMRI SPATIAL MAP RECONSTRUCTION');

% =========================================================================
%                    Load eigenmodes and empirical data           
% =========================================================================

% Load relevant eigenmodes

% Geometric eigenmodes: 
% eigenmodes = load(sprintf('%sdata/eigenmodes/fsLR_32k_midthickness-%s_emode_200.txt', project_dir, hemisphere));

% Connectome eigenmodes:
% eigenmodes = load(sprintf('%sdata/eigenmodes/connectome_eigenmodes-%s_200.mat', project_dir, hemisphere));
% eigenmodes = eigenmodes.eig_vec;

% EDR eigenmodes:
% eigenmodes = load(sprintf('%sdata/network_eigenmodes/EDR_eigenmodes-%s_200.mat', project_dir, hemisphere));
% eigenmodes = eigenmodes.eig_vec;

% Load task fMRIdata
data = load(sprintf('/data/pt_02994/data/empirical/S255_tfMRI_ALLTASKS_raw_%s', hemisphere));
data_to_reconstruct = data.zstat;

num_modes = size(eigenmodes, 2);

field_names = fieldnames(data_to_reconstruct);
num_tasks = numel(field_names);

disp('Eigenmodes loaded.');

% =========================================================================
%       Calculate reconstruction beta coefficients for each subject
% =========================================================================

recon_beta_struct = struct(); % store reconstruction weights

recon_corr_vertex_struct = struct(); % store reconstruction accuracy on vertices
recon_corr_parc_struct = struct(); % store reconstruction accuracy on parcellations

mse_vertex_struct = struct();
mse_parc_struct = struct();

for task_id = 1:num_tasks
    task_name = field_names{task_id};
    task_data = data_to_reconstruct.(task_name);

    num_subjects = size(task_data, 2);

    recon_beta = zeros(num_modes, num_modes, num_subjects);
    
    % For each subject, decompose eigenmode data
    for subj = 1:num_subjects
        fprintf('Calculating coefficients for subject %s.', subj);
        for mode = 1:num_modes
            basis = eigenmodes(cortex_ind, 1:mode);
            recon_beta(1:mode, mode, subj) = calc_eigendecomposition(task_data(cortex_ind, subj), basis, 'matrix');
        end 
    end 
    
    % Store reconstruction beta coefficients of subjects for a task
    recon_beta_struct.(task_name) = recon_beta;

    disp('Reconstruction beta coefficients calculated.');

% =========================================================================
%     Calculate reconstruction accuracy for each subject   
% =========================================================================

    recon_corr_vertex = zeros(num_subjects, num_modes);
    recon_corr_parc = zeros(num_subjects, num_modes);

    mse_vertex = zeros(num_subjects, 1);
    mse_parc = zeros(num_subjects, 1);

    parc_name = 'Glasser360'; % other parcellations can be used as well
    parc = dlmread(sprintf('%sdata/parcellations/fsLR_32k_%s-%s.txt', project_dir, parc_name, hemisphere));

    for subj = 1:num_subjects
        fprintf('Calculating reconstruction accuracy for for subject %s.', subj);
        for mode = 1:num_modes

            recon_temp = eigenmodes(cortex_ind, 1:mode)*recon_beta(1:mode, mode, subj);
            recon_corr_vertex(subj, mode) = corr(task_data(cortex_ind, subj), recon_temp);

            if mode == num_modes
                mse_vertex(subj) = nanmean((task_data(cortex_ind, subj) - recon_temp).^2);
            end
    
            recon_temp = eigenmodes(:, 1:mode)*recon_beta(1:mode, mode, subj);
            recon_corr_parc(subj, mode) = corr(calc_parcellate(parc, task_data(:, subj)), calc_parcellate(parc, recon_temp));

            if mode == num_modes
                mse_parc(subj) = nanmean((calc_parcellate(parc, task_data(:, subj)) - calc_parcellate(parc, recon_temp)).^2);
            end

        end
    end
    
    recon_corr_vertex_struct.(task_name) = recon_corr_vertex;
    recon_corr_parc_struct.(task_name) = recon_corr_parc;

    mse_vertex_struct.(task_name) = mse_vertex;
    mse_parc_struct.(task_name) = mse_parc;

end

disp('Reconstruction accuracies calculated.');

% =========================================================================
%                               Calculate MSE                     
% =========================================================================

fprintf('Average MSE for vertex-level reconstruction calculated: %s', nanmean(mse_vertex_struct.(task_name)));
fprintf('Average MSE for parcellation-level reconstruction calculated: %s', nanmean(mse_parc_struct.(task_name)));

% =========================================================================
%                      Some visualizations of results                      
% =========================================================================

% Reconstruction accuracy vs number of modes at vertex and parcellated levels
figure('Name', 'tfMRI reconstruction - accuracy');
hold on;
plot(1:num_modes, recon_corr_vertex_struct.(task_name), 'k-', 'linewidth', 2, 'displayname', 'vertex')
plot(1:num_modes, recon_corr_parc_struct.(task_name), 'b-', 'linewidth', 2, 'displayname', 'parcellated')
hold off;
leg = legend('fontsize', 12, 'location', 'southeast', 'box', 'off');
set(gca, 'fontsize', 10, 'ticklength', [0.02 0.02], 'xlim', [1 num_modes], 'ylim', [0 1])
xlabel('number of modes', 'fontsize', 12)
ylabel('reconstruction accuracy', 'fontsize', 12)

% Reconstruct spatial map
N = num_modes;
surface_to_plot = surface_midthickness;
data_to_plot = eigenmodes(:, 1:N)*recon_beta(1:N,N);
medial_wall = find(cortex==0);
with_medial = 1;

fig = draw_surface_bluewhitered_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = sprintf('tfMRI reconstruction - surface map using %i modes', N);
