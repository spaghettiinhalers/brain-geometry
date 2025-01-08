%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculates high-resolution connectome-based eigenmodes for:
%%% (1) a high-resolution structural connectome
%%% (2) a synthetic stochastically-generated exponential distance rule (EDR)-based connectome
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
surface_to_analyze = surface_midthickness;
num_modes = 200;

disp('Surface files loading complete.');



%% CALCULATE CONNECTOME EIGENMODES

disp('BEGINNING CONNECTOME EIGENMODES CALCULATION');

% =========================================================================
%                   Generate surface local connectivity                    
% =========================================================================

surface_connectivity = calc_surface_connectivity(surface_to_analyze);

% Remove vertices corresponding to the medial wall
surface_connectivity = surface_connectivity(cortex_ind, cortex_ind);

disp('Surface local connectivity generated.');

% =========================================================================
%                      Generate synthetic connectome                       
% =========================================================================

% Replace the connectome variable below with your empirical data
% However, make sure that you remove the vertices corresponding to the medial wall first
connectome = load(sprintf('%sdata/empirical/S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat', project_dir));
connectome = connectome.avgSC_L;

% Threshold connectome
m = nnz(surface_connectivity);
n = nnz(connectome);
threshold = 4*m/n;
connectome_threshold = threshold_edges_proportional(connectome, threshold, 'weakest');

disp('Synthetic connectome generated.');

% =========================================================================
%            Combine surface local connectivity and connectome             
% =========================================================================

surface_with_connectome = surface_connectivity + (connectome_threshold>0);
surface_with_connectome(surface_with_connectome>0) = 1;

% Verify connection density: 
N = size(surface_with_connectome,1);
connection_density = nnz(surface_with_connectome) / (N * (N-1)) * 100; % percentage

fprintf('Surface local connectivity and connectome combined with a connection density of %s percent.', connection_density);

% =========================================================================
%                           Calculate the modes                            
% =========================================================================

[eig_vec_temp, eig_val, Laplacian] = calc_network_eigenmode(surface_with_connectome, num_modes);

save(sprintf('%sdata/laplacians/Laplacian_connectome_%s.mat', project_dir, hemisphere), 'Laplacian', '-v7.3');

% Bring back medial wall vertices with zero values
eig_vec = zeros(num_vertices, num_modes);
eig_vec(cortex_ind,:) = eig_vec_temp(:,1:num_modes);
save(sprintf('%sdata/eigenmodes/connectome_eigenmodes-%s_%i.mat', project_dir, hemisphere, num_modes), 'eig_val', 'eig_vec', '-v7.3')

disp('Connectome modes calculated.');



%% CALCULATE EDR EIGENMODES

disp('BEGINNING EDR EIGENMODES CALCULATION');

% =========================================================================
% Calculate connection lengths of surface vertices without the medial wall                    
% =========================================================================

euclidean_dist = squareform(pdist(surface_to_analyze.vertices(cortex_ind,:)));

% Change function according to line of best fit
sigmoid = @(p, x) p(1) * x + p(2) + p(3) ./ (1 + exp(-p(4) * (x - p(5))));
p = [1, 5, 45, 0.15, 35];

% Calculate connection length
surface_dist = sigmoid(p, euclidean_dist);

disp('Connection lengths calculated.');

% =========================================================================
%                    Generate synthetic EDR connectome                     
% =========================================================================

% Probability function
Pspace_func = @(scale, distance) exp(-scale*distance);

% Generate pseudorandom numbers to compare with probability function
rng(1) % set for reproducibility
rand_prob = rand(size(surface_dist));
rand_prob = triu(rand_prob,1) + triu(rand_prob,1)';
rand_prob(1:1+size(rand_prob,1):end) = 1;

% Calculate empirical scale exponent parameter based on structural connectivity data
connectome_weights = connectome(:);
connectome_distances = surface_dist(:);

exp_function = @(alpha, d) exp(-alpha * d);
error_function = @(alpha) sum((connectome_weights - exp_function(alpha, connectome_distances)).^2);

alpha_init = 0.120;
alpha_empirical = fminsearch(error_function, alpha_init);

% Calculate probability
scale = alpha_empirical;
Pspace = Pspace_func(scale, surface_dist);
Pspace(1:1+size(Pspace,1):end) = 0;
Pspace = Pspace/max(Pspace(:));
    
% Generate EDR connectome
connectome = double(rand_prob < Pspace);
connectome(1:1+size(connectome,1):end) = 0;
   
disp('Synthetic EDR connectome generated.');

% =========================================================================
%                           Calculate the modes                            
% =========================================================================

[eig_vec_temp, eig_val, Laplacian] = calc_network_eigenmode(connectome, num_modes);

% Bring back medial wall vertices with zero values
eig_vec = zeros(num_vertices, num_modes);
eig_vec(cortex_ind,:) = eig_vec_temp(:,1:num_modes);
save(sprintf('%sdata/eigenmodes/EDR_eigenmodes-%s_%i.mat', project_dir, hemisphere, num_modes), 'eig_val', 'eig_vec', '-v7.3')

disp('EDR modes calculated.');