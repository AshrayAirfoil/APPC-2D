%% CORRECTED SYMMETRIC NOZZLE OPTIMIZER (V2 - Fixed)
% -------------------------------------------------------------------------
% OBJECTIVE:
% Design a Symmetric Nozzle Housing with:
% 1. OUTSIDE: Smooth, continuous round curve (Suction Side).
% 2. INSIDE: Flat cowling channel (Pressure Side).
% 3. INLET: Sharp Edge (x=0).
% 4. EXIT: Round Edge (x=1).
% -------------------------------------------------------------------------

close all; clear; clc;

% =========================================================================
% 1. CONFIGURATION
% =========================================================================

% --- Physics ---
CT_design   = 10;
Alpha_deg   = 6;         % Fixed Angle of Attack
Re          = 1e6;

% --- Geometry Settings ---
inlet_gap   = 0.45;      % Gap between Sharp Leading Edges (Inlet)
chord_len   = 1.0;       % Normalized Chord

% --- Search Space ---
search_folder = './Airfoils/ONR-Coords'; 
% Pitch Angles: Rotate Top Foil DOWN (-) to converge the tail
pitch_angles = 0:-1:-15; 

% =========================================================================
% 2. INITIALIZATION
% =========================================================================
addpath(genpath(pwd));
temp_dir = './Temp_Correct_Opt';
if ~exist(temp_dir, 'dir'), mkdir(temp_dir); end

candidates = struct('name', {}, 'x', {}, 'y', {});
cnt = 0;

if isfolder(search_folder)
    files = dir(fullfile(search_folder, '*.dat'));
    fprintf('Processing %d candidates...\n', length(files));
    
    for i = 1:length(files)
        fname = files(i).name;
        [raw_x, raw_y] = read_coords_robust(fullfile(files(i).folder, fname));
        
        if isempty(raw_x), continue; end
        
        % --- A. STANDARDIZE ORIENTATION ---
        % 1. Normalize to 0..1
        raw_x = raw_x - min(raw_x);
        scale = 1.0 / max(raw_x);
        raw_x = raw_x * scale;
        raw_y = (raw_y - mean(raw_y)) * scale;
        
        % 2. Detect Geometry (Sharp Front vs Round Front)
        % We need the "Curved/Suction" side to be UP (Positive Y).
        % We need the "Sharp" edge to be at x=0 (Inlet).
        
        [is_standard_shape, ~, ~] = check_shape(raw_x, raw_y);
        
        if is_standard_shape
            % Standard = Round Front (0), Sharp Back (1).
            % FLIP X to make Sharp Front (0).
            final_x = 1.0 - raw_x;
            final_y = raw_y;
            
            % RE-ORDER POINTS to maintain loop continuity for solver
            final_x = flipud(final_x);
            final_y = -flipud(final_y);
        else
            % Already Sharp Front.
            final_x = raw_x;
            final_y = -raw_y;
        end
        
        cnt = cnt + 1;
        candidates(cnt).name = fname;
        candidates(cnt).x = final_x * chord_len;
        candidates(cnt).y = final_y * chord_len;
    end
end

fprintf('Candidates prepared: %d\n', cnt);
fprintf('%s\n', repmat('-',1,95));
fprintf('%-20s | %-8s | %-8s | %-8s | %-8s | %-12s\n', ...
    'Airfoil', 'Pitch', 'Cl', 'L/D', 'ExitGap', 'Status');
fprintf('%s\n', repmat('-',1,95));

results = struct([]);
sim_idx = 0;

% =========================================================================
% 3. MAIN OPTIMIZATION LOOP
% =========================================================================

for i = 1:length(candidates)
    cand = candidates(i);
    
    for theta = pitch_angles
        
        % --- A. TOP FOIL CONSTRUCTION ---
        % 1. Get Coordinates (Sharp @ 0, Round @ 1)
        nx = cand.x; 
        ny = cand.y;
        
        % 2. Pitch Rotation (Pivot at Sharp Inlet x=0)
        % Rotate DOWN (Negative) to converge tail
        rad = theta * pi / 180;
        [nx_rot, ny_rot] = rotate_coords(nx, ny, rad, 0, 0);
        
        % 3. Shift Vertical
        % Move Sharp Inlet to +Gap/2
        ny_top = ny_rot + (inlet_gap / 2);
        nx_top = nx_rot;
        
        % --- B. BOTTOM FOIL CONSTRUCTION ---
        % 1. MIRROR the Top Foil
        % This puts the Curved Surface (Top of Top Foil) onto the Bottom.
        nx_bot = nx_top;
        ny_bot = -ny_top;
        
        % Reverse point order for Bottom to maintain CCW/CW loop integrity
        nx_bot = flipud(nx_bot);
        ny_bot = flipud(ny_bot);
        
        % --- C. GEOMETRY VALIDATION ---
        
        % 1. Convergence Check (Exit Gap)
        [~, idx_top_te] = max(nx_top);
        [~, idx_bot_te] = max(nx_bot);
        
        y_exit_top = ny_top(idx_top_te);
        y_exit_bot = ny_bot(idx_bot_te);
        
        gap_exit = y_exit_top - y_exit_bot;
        
        % Constraint: Strictly Convergent
        if gap_exit >= inlet_gap
            continue; % Divergent geometry
        end
        
        % Constraint: Open Exit
        if gap_exit < 0.10
            continue; % Choked
        end
        
        % Constraint: Intersection
        if min(ny_top) <= 0
            continue; % Intersection
        end
        
        % --- D. SIMULATION ---
        try
            f_base = fullfile(temp_dir, 'base.dat');
            f_nac  = fullfile(temp_dir, 'nac.dat');
            write_file(f_base, nx_bot, ny_bot);
            write_file(f_nac, nx_top, ny_top);
            
            % Surface 1 = Top (Upper Wake)
            % Surface 2 = Bottom (Lower Wake)
            sFiles(1).name = f_nac;
            sFiles(2).name = f_base;
            if length(sFiles)>2, sFiles(3:end)=[]; end
            
            % Run Solver
            [~, ~] = evalc('main(sFiles, Alpha_deg, CT_design)');
            [surfaces, ~, ~, ~] = main(sFiles, Alpha_deg, CT_design);
            
            % Metrics
            Cl = get_cl(surfaces, Alpha_deg);
            
            % Total Drag Estimate (Sum of both)
            Cd = est_drag(nx_top, ny_top, Re) + est_drag(nx_bot, ny_bot, Re);
            
            sim_idx = sim_idx + 1;
            results(sim_idx).name   = cand.name;
            results(sim_idx).theta  = theta;
            results(sim_idx).Cl     = Cl;
            results(sim_idx).LD     = Cl/Cd;
            results(sim_idx).exit   = gap_exit;
            results(sim_idx).nx     = nx_top; results(sim_idx).ny = ny_top;
            results(sim_idx).bx     = nx_bot; results(sim_idx).by = ny_bot;
            
            if Cl > 1.2
                 fprintf('%-20s | %-8.1f | %-8.3f | %-8.2f | %-8.3f | OK\n', ...
                    cand.name, theta, Cl, Cl/Cd, gap_exit);
            end
            
        catch
            % Fail
        end
    end
end

rmdir(temp_dir, 's');

% =========================================================================
% 4. RESULTS & VISUALIZATION
% =========================================================================

if isempty(results)
    error('No valid geometries found. Try increasing Inlet Gap.');
else
    [~, idx] = max([results.Cl]);
    best = results(idx);
    
    fprintf('\n>>> BEST CORRECTED DESIGN <<<\n');
    fprintf('Airfoil:      %s\n', best.name);
    fprintf('Pitch Angle:  %.1f deg (Converging)\n', best.theta);
    fprintf('Max Cl:       %.4f\n', best.Cl);
    fprintf('Convergence:  %.2f (Inlet) -> %.3f (Exit)\n', inlet_gap, best.exit);
    
    % --- PLOT ---
    figure('Name', 'Corrected Nozzle Geometry', 'Color', 'w');
    hold on;
    
    % Fill Geometry
    fill(best.nx, best.ny, [0.2 0.2 0.2], 'EdgeColor', 'k'); % Top
    fill(best.bx, best.by, [0.2 0.2 0.2], 'EdgeColor', 'k'); % Bottom
    
    % Visual Guides
    yline(0, '-.k', 'Centerline');
    xline(0, ':r', 'Inlet (Sharp)');
    xline(1, ':b', 'Exit (Round)');
    
    % Flow Arrow
    quiver(-0.2, 0, 0.4, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    text(-0.25, 0.1, 'FLOW', 'Color', 'r', 'FontWeight', 'bold');
    
    % Annotation
    text(0.5, max(best.ny)+0.1, 'Curved Outside', 'HorizontalAlignment', 'center');
    text(0.5, 0.05, 'Cowling Inside', 'HorizontalAlignment', 'center');
    
    axis equal; grid on;
    xlabel('x/c'); ylabel('y/c');
    title(sprintf('Corrected Nozzle: Round Outer Shell, Flat Inner Channel (Cl=%.3f)', best.Cl));
end

% =========================================================================
% 5. HELPER FUNCTIONS (FIXED)
% =========================================================================

function [is_standard, t_front, t_back] = check_shape(x, y)
    % Compare thickness at 10% and 90%
    % Standard = Round Front (0), Sharp Back (1)
    
    [x_u, I] = unique(x); y_u = y(I);
    tol = 0.05;
    y_f = y(x > 0.1-tol & x < 0.1+tol);
    y_b = y(x > 0.9-tol & x < 0.9+tol);
    
    % FIX: Assign defaults if empty
    if isempty(y_f) || isempty(y_b)
        is_standard = true; 
        t_front = 0; t_back = 0;
        return; 
    end
    
    t_front = max(y_f) - min(y_f);
    t_back  = max(y_b) - min(y_b);
    
    % Assign to outputs
    % If Front is thicker than Back, it's Standard
    is_standard = (t_front > t_back);
end

function [xr,yr] = rotate_coords(x,y,r,cx,cy)
    x0=x-cx; y0=y-cy; xr=x0*cos(r)-y0*sin(r)+cx; yr=x0*sin(r)+y0*cos(r)+cy;
end

function Cl = get_cl(S,a)
    rad=a*pi/180; R=[cos(rad) -sin(rad); sin(rad) cos(rad)]; Cn=0; Ca=0;
    for k=1:length(S)
        c=repmat([0.25;0],1,length(S(k).endPoints)); pts=(R*(S(k).endPoints'-c)+c)';
        pt1=pts(1:end-1,:); pt2=pts(2:end,:); th=atan2(pt2(:,2)-pt1(:,2),pt2(:,1)-pt1(:,1));
        co=pt1+(pt2-pt1)/2; tgt=find(co(:,1)==min(co(:,1)),1);
        for j=1:length(S(k).CP)
            F=S(k).CP(j)*abs(S(k).DL(j));
            if j<=tgt, Cn=Cn+F*cos(th(j)); Ca=Ca-F*sin(th(j));
            else, Cn=Cn-F*cos(th(j)); Ca=Ca+F*sin(th(j)); end
        end
    end
    Fw=R*[Cn;Ca]; Cl=Fw(1);
end

function Cd = est_drag(x,y,Re)
    p=sum(sqrt(diff(x).^2+diff(y).^2)); c=max(x)-min(x); t=max(y)-min(y);
    Cd=(0.074/Re^0.2)*(p/c)*(1+2*(t/c)+60*(t/c)^4);
end

function [x,y] = read_coords_robust(f)
    try d=importdata(f); if isstruct(d),d=d.data; end; x=d(:,1); y=d(:,2); catch, x=[]; y=[]; end
end

function write_file(f,x,y)
    fid=fopen(f,'w'); for k=1:length(x), fprintf(fid,'%.6f %.6f\n',x(k),y(k)); end; fclose(fid);
end
