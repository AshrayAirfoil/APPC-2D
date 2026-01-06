%% Angle of Attack Sweep
close all; clear; clc;

% 1. SETUP PATHS
% Adds all subfolders (Airfoils, colormaps, mesh2d, etc.) to MATLAB path
addpath(genpath(pwd)); 

% 2. CONFIGURATION
% Adjust these filenames to match your actual files
vecAng = 0;
% Note: Using fullfile() fixes slash direction issues (Windows vs Mac/Linux)
file_nacelle = fullfile('Airfoils', 'ONR-Coords', sprintf('nacelleVec%g.dat', vecAng));
file_main    = fullfile('Airfoils', 'ONR-Coords', sprintf('mainVec%g.dat', vecAng));
%file_krueger = fullfile('Airfoils', 'ONR-Coords', 'krueger.dat');

% Check if files exist before starting
if ~isfile(file_nacelle), error('File not found: %s', file_nacelle); end
if ~isfile(file_main),    error('File not found: %s', file_main); end

% Define Geometry Object
surfaceFiles(1).name = file_nacelle;
surfaceFiles(2).name = file_main;
%surfaceFiles(3).name = file_krueger;

CT = 10;                 % Fixed Thrust Coefficient
alpha_range = -5:1:20;   % Sweep Range (Degrees)

% 3. RUN SIMULATION
fprintf('Starting Alpha Sweep (CT=%g)...\n', CT);
fprintf('%-10s %-10s %-20s\n', 'Alpha', 'Cl', 'Status');
fprintf('%s\n', repmat('-',1,45));

results = [];

for i = 1:length(alpha_range)
    alpha = alpha_range(i);
    
    try
        % --- CALL SOLVER ---
        [surfaces, wake, iter, chord] = main(surfaceFiles, alpha, CT);
        
        % --- CALCULATE CL (Integrated from Cp) ---
        % Convert Alpha to Radians for rotation matrix
        alRad = alpha * pi / 180;
        R = [cos(alRad) -sin(alRad); sin(alRad) cos(alRad)];
        
        Cn = 0; Ca = 0; % Normal and Axial Force Accumulators
        
        % Iterate through all surfaces (Nacelle, Main, Krueger)
        for k = 1:length(surfaces)
            % Derotate geometry to body frame
            center = repmat([0.25;0], 1, length(surfaces(k).endPoints(:,1)));
            pts_body = (R*(surfaces(k).endPoints' - center) + center)';
            
            % Recalculate panel angles in body frame
            pt1 = pts_body(1:surfaces(k).m, :);
            pt2 = pts_body(2:surfaces(k).n, :);
            theta_body = atan2((pt2(:,2)-pt1(:,2)), (pt2(:,1)-pt1(:,1)));
            
            % Leading edge split index
            % Note: We use 'co' (collocation points) to find the LE
            co_body = pt1 + (pt2 - pt1)/2;
            target = find(co_body(:,1) == min(co_body(:,1)), 1);
            
            % Integrate Pressure
            for j = 1:length(surfaces(k).CP)
                cp = surfaces(k).CP(j);
                dl = abs(surfaces(k).DL(j));
                th = theta_body(j);
                
                % Standard convention: Upper surface adds to Normal, Lower subtracts
                % (or vice versa depending on panel orientation).
                % Based on designDriver logic:
                if j <= target
                    Cn = Cn + cp * dl * cos(th);
                    Ca = Ca - cp * dl * sin(th);
                else
                    Cn = Cn - cp * dl * cos(th);
                    Ca = Ca + cp * dl * sin(th);
                end
            end
        end
        
        % Rotate Forces to Wind Frame
        force_wind = R * [Cn; Ca];
        current_Cl = force_wind(1);
        
        % Store Results
        results(i).alpha = alpha;
        results(i).Cl = current_Cl;
        fprintf('%-10d %-10.4f %-20s\n', alpha, current_Cl, 'Success');
        
    catch ME
        % --- ERROR HANDLING ---
        % This will print the specific error message instead of just "Failed"
        fprintf('%-10d %-10s %-20s\n', alpha, 'NaN', ['Error: ' ME.message]);
        results(i).alpha = alpha;
        results(i).Cl = NaN;
    end
end

% 4. PLOT AND SUMMARY
valid_results = results(~isnan([results.Cl]));

if ~isempty(valid_results)
    [max_Cl, idx] = max([valid_results.Cl]);
    best_alpha = valid_results(idx).alpha;
    
    fprintf('\n=== RESULTS ===\n');
    fprintf('Best Alpha: %d deg | Max Cl: %.4f\n', best_alpha, max_Cl);
    
    figure('Color','w');
    plot([valid_results.alpha], [valid_results.Cl], '-o', 'LineWidth', 2);
    grid on;
    xlabel('Angle of Attack (deg)');
    ylabel('Lift Coefficient (C_l)');
    title(sprintf('Cl vs Alpha (CT=%g)', CT));
else
    fprintf('\nAll runs failed. Check error messages above.\n');
end
