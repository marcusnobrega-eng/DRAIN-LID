%% =========================================================================
%                    Four-Case Hydrus vs. Model Comparison
%                 2 rows x 4 cols (per case: h-profile | theta-profile)
% =========================================================================
close all; clear; clc;

%% ---------------------- CASE DEFINITIONS (EDIT THESE) --------------------
% For each case provide: xlsPath (with sheets 'Head' and 'Moisture'),
% hydrusOut (.out file), label (for titles), and pdf output (optional)
cases(1).label     = 'Example 1';
cases(1).xlsPath   = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Modeling_Results\Example1_Infiltration\Data\SimulationResults.xlsx';
cases(1).hydrusOut = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Modeling_Results\Example1_Infiltration\Data\Example_1_Hydrus.out';

cases(2).label     = 'Example 2';
cases(2).xlsPath   = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Modeling_Results\Example2_Clay_Loam_Soil\Data\SimulationResults.xlsx';
cases(2).hydrusOut = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Modeling_Results\Example2_Clay_Loam_Soil\Data\Example_2_Hydrus.out';

cases(3).label     = 'Example 3';
cases(3).xlsPath   = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Modeling_Results\Example3_CapillaryRise\Data\SimulationResults.xlsx';
cases(3).hydrusOut = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Modeling_Results\Example3_CapillaryRise\Data\Example_3_CapillaryRise.out';

cases(4).label     = 'Example 4';
cases(4).xlsPath   = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Modeling_Results\Example4_TopNeumann_Sandy\Data\SimulationResults.xlsx';
cases(4).hydrusOut = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Modeling_Results\Example4_TopNeumann_Sandy\Data\Example4_TopNeumann_Sandy.out';

% Optional: output PDF for the whole 2x4 figure
pdfOutput = 'C:\Users\marcu\Documents\GitHub\LID_Tool\Modeling_Results\Hydrus_Comparison\FourCase_Comparison.pdf';

%% ------------------------- GLOBAL PLOTTING SETTINGS ----------------------
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontName','Montserrat');
set(groot,'defaultAxesFontSize',14);
set(groot,'defaultAxesLineWidth',1.5);

% Distinguishable 8-color palette (Bright)
hexColors = { '#003a7d', '#008dff', '#ff73b6', '#c701ff', ...
              '#4ecb8d', '#ff9d3a', '#f9e858', '#d83034' };

% Convert each hex code to an [R G B] row in 0–1 range
COL = zeros(numel(hexColors), 3); % preallocate for speed
for i = 1:numel(hexColors)
    COL(i,:) = hex2rgb(hexColors{i}); % each row is a valid RGB triplet
end


% Use 5 key instants (fractions of total simulation time)
key_fracs = [0.1, 0.25, 0.5, 0.75, 1.0];

% Line/marker styles
lw  = 2.5;              % thick borders
styHyd = '-';           % Hydrus
styMod = ':';           % Your model

%% === Axis limits (EDIT HERE) ===
xlim_h   = [-1.0 0.0];        % head x-axis limits
xlim_th  = [0.0 0.5];         % theta x-axis limits
ylim_all = [-1.0 0.0];        % depth y-axis limits

%% === Figure size (A4 width, half A4 height) ===
fig = figure('Color','w','Units','inches','Position',[1, 1, 11.69, 5.8]);


%% ------------------------------- FIGURE ---------------------------------

% A4 landscape
fig = figure('Color','w','Units','inches','Position',[1, 1, 11.69, 8.27]);

% We will fill subplots in order: 
% Row 1:  (1) h1 | (2) th1 | (3) h2 | (4) th2
% Row 2:  (5) h3 | (6) th3 | (7) h4 | (8) th4

nCases = numel(cases);
for c = 1:nCases
    %% ---------------------- Load model data for case c -------------------
    T_head   = readtable(cases(c).xlsPath, 'Sheet', 'Head');
    time_s   = T_head.Time_s;
    model_h  = table2array(T_head(:, startsWith(T_head.Properties.VariableNames,'h_Node'))); % [nt x nz]

    T_th     = readtable(cases(c).xlsPath, 'Sheet', 'Moisture');
    model_th = table2array(T_th(:, startsWith(T_th.Properties.VariableNames,'theta_Node'))); % [nt x nz]

    % If you have 'params' per case in MAT, load it here; else infer domain depth
    % load('Examples\Example1_Infiltration.mat'); % <-- adapt per case if needed
    % For depth vector z, we assume evenly spaced nodes over [-L, 0]
    % If L unknown, infer from your model's existing convention:
    L = 1.0;  % <-- EDIT if known; otherwise replace with your params.L
    nz = size(model_h,2);
    z  = linspace(-L, 0, nz);

    %% ------------------ Parse HYDRUS .out for h & theta ------------------
    [tH, H_hyd, TH_hyd] = parseHydrusOut(cases(c).hydrusOut);

    %% ------------------ Interpolate HYDRUS to model times ----------------
    % Keep matching number of vertical nodes with your model
    H_int  = interp1(tH, H_hyd(:,1:nz),  time_s, 'linear', 'extrap');   % [nt x nz]
    TH_int = interp1(tH, TH_hyd(:,1:nz), time_s, 'linear', 'extrap');   % [nt x nz]

    %% ------------------ Choose key indices (common style) ----------------
    key_idx = max(1, round(key_fracs * numel(time_s)));
    key_idx = unique(key_idx, 'stable');

    %% ------------------ Compute final-step RMSEs (profiles) --------------
    last = numel(time_s);
    rmse_h_final  = RMSE(model_h(last,:),  flip(H_int(last,:)));
    rmse_th_final = RMSE(model_th(last,:), flip(TH_int(last,:)));

    %% ----------------------------- PLOTS ---------------------------------
    % Column indices in the 2x4 grid for this case
    switch c
        case 1, colH = 1; colT = 2; row = 1;
        case 2, colH = 3; colT = 4; row = 1;
        case 3, colH = 1; colT = 2; row = 2;
        case 4, colH = 3; colT = 4; row = 2;
    end
    spH = subplot(2,4,(row-1)*4 + colH); hold(spH,'on');
    spT = subplot(2,4,(row-1)*4 + colT); hold(spT,'on');

    % Cycle through 5 time snapshots with distinct colors
    for k = 1:numel(key_idx)
        idx = key_idx(k);
        colk = COL(mod(k-1,size(COL,1))+1,:);

        % Pressure head
        plot(spH, flip(H_int(idx,:)),  z, styHyd, 'LineWidth', lw, 'Color', colk); % Hydrus
        plot(spH, model_h(idx,:),      z, styMod, 'LineWidth', lw, 'Color', colk); % Model

        % Moisture
        plot(spT, flip(TH_int(idx,:)), z, styHyd, 'LineWidth', lw, 'Color', colk);
        plot(spT, model_th(idx,:),     z, styMod, 'LineWidth', lw, 'Color', colk);
    end

    % Axis labels/titles
    xlabel(spH, '$h$ [m]');
    ylabel(spH, 'Depth [m]');
    % title(spH,  sprintf('\\textbf{%s: $h$ profiles}', cases(c).label));

    xlabel(spT, '$\theta$ [cm$^3$/cm$^3$]');
    ylabel(spT, 'Depth [m]');
    % title(spT,  sprintf('\\textbf{%s: $\\theta$ profiles}', cases(c).label));

    % Formatting (minor ticks inside, major outside, longer)
    formatAxes(spH);
    formatAxes(spT);
    

    

    % Annotate final-step RMSEs (only, per your request)
    txtH  = sprintf('\\textbf{RMSE$_{t_{end}}$ = %.3g m}', rmse_h_final);
    txtTh = sprintf('\\textbf{RMSE$_{t_{end}}$ = %.3g}',   rmse_th_final);
    annotationInAxes(spH, txtH);
    annotationInAxes(spT, txtTh);

    % Optional: lock x-lims per variable across all cases (uncomment & set)
    % xlim(spH, [-2 0.2]);     % example limits for h
    % xlim(spT, [0 0.5]);      % example limits for theta
end

% Tighten layout
set(gcf,'Renderer','painters'); % vector-friendly
tiledlayout = [];               % (placeholder if you later move to tiledlayout)

% === Force Montserrat on all text elements ===
set(findall(fig,'-property','FontName'), 'FontName','Montserrat');
set(findall(fig,'-property','FontUnits'), 'FontUnits','points');
set(findall(fig,'-property','FontSizeMode'), 'FontSizeMode','manual');

% ---------------------------- Export (optional) ---------------------------
if ~isempty(pdfOutput)
    exportgraphics(fig, pdfOutput, 'ContentType','vector', ...
        'BackgroundColor','white', 'Resolution',300);
end

%% === Add top legend for time ratio mapping ===
legendAx = axes('Position',[0.1, 0.93, 0.8, 0.05]); % [x y w h]
hold(legendAx,'on');

barHalf = 0.2; % half-length of each bar from center

for k = 1:numel(key_fracs)
    colk = COL(mod(k-1, size(COL,1)) + 1, :);
    % Draw shorter bar
    plot(legendAx, [k-barHalf, k+barHalf], [1 1], '-', ...
        'Color', colk, 'LineWidth', 4);
    % Time fraction text BELOW the bar
    text(k, 0.4, sprintf('%.2f', key_fracs(k)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontName','Montserrat', 'FontSize',12);
end

% Formatting for the legend axis
xlim(legendAx,[0.5 numel(key_fracs)+0.5]);
ylim(legendAx,[0 1.5]);
legendAx.Visible = 'off';

% $t/t_d$ label ABOVE all bars, centered
text(mean(xlim(legendAx)), 1.3, '$t/t_d$', ...
    'Interpreter','latex', ...
    'FontWeight','bold', 'FontSize',13, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'Parent', legendAx);



%% ============================== HELPERS ==================================
function [t, H, TH] = parseHydrusOut(fname)
% Parse HYDRUS column output (.out) to times, pressure head, and theta
% Returns:
%   t  [nt x 1]
%   H  [nt x nz]  (pressure head)
%   TH [nt x nz]  (moisture)

    fid = fopen(fname,'r'); assert(fid>0, 'Cannot open %s', fname);
    C = textscan(fid, '%s', 'Delimiter','\n', 'Whitespace','');
    fclose(fid); lines = C{1};

    heads_blocks = {}; thetas_blocks = {}; times = [];
    ch = []; ct = [];
    for i = 1:numel(lines)
        s = strtrim(lines{i});
        if startsWith(s, 'Time:')
            if ~isempty(ch)
                heads_blocks{end+1}  = ch; %#ok<*AGROW>
                thetas_blocks{end+1} = ct;
            end
            ch = []; ct = [];
            tv = sscanf(s, 'Time: %f');
            times(end+1) = tv;
        elseif ~isempty(regexp(s, '^\d+', 'once')) && ~contains(s,'Node')
            tok = sscanf(s, '%f');
            if numel(tok) >= 4
                ch(end+1) = tok(3);  % Column 3: pressure head
                ct(end+1) = tok(4);  % Column 4: theta
            end
        end
    end
    if ~isempty(ch)
        heads_blocks{end+1}  = ch;
        thetas_blocks{end+1} = ct;
    end
    % Ensure rectangular arrays [nt x nz]
    H  = padAndStack(heads_blocks);
    TH = padAndStack(thetas_blocks);
    t  = times(:);
end

function A = padAndStack(blocks)
% Convert a cell array of row vectors into a matrix [nBlocks x maxLen]
    n = numel(blocks);
    L = max(cellfun(@numel, blocks));
    A = nan(n, L);
    for i = 1:n
        v = blocks{i};
        A(i,1:numel(v)) = v;
    end
end

function r = RMSE(sim, obs)
% Root-mean-square error between two vectors (ignores NaNs)
    mask = ~(isnan(sim) | isnan(obs));
    sim = sim(mask); obs = obs(mask);
    r = sqrt(mean((sim - obs).^2));
end

function rgb = hex2rgb(hex)
% Convert '#RRGGBB' to [r g b] in 0..1
    if hex(1) == '#', hex = hex(2:end); end
    rgb = [hex2dec(hex(1:2)), hex2dec(hex(3:4)), hex2dec(hex(5:6))] / 255;
end

function formatAxes(ax)
% Uniform axis formatting for all subplots
    set(ax, 'Box','on', ...
        'TickDir','out', ...           % major ticks outside
        'TickLength',[0.025 0.02], ... % slightly longer major ticks
        'XMinorTick','on', ...
        'YMinorTick','on', ...
        'Layer','top');                % ticks above grid
    ax.MinorGridLineStyle = 'none';    % minor ticks only, no minor grid
    grid(ax, 'on');
end


function annotationInAxes(ax, txt)
    xlim_ = xlim(ax); ylim_ = ylim(ax);
    x = xlim_(1) + 0.97*(xlim_(2)-xlim_(1));
    y = ylim_(1) + 0.95*(ylim_(2)-ylim_(1));
    t = text(ax, x, y, txt, ...
        'HorizontalAlignment','right', 'VerticalAlignment','top', ...
        'BackgroundColor',[1 1 1], 'EdgeColor',[0 0 0], 'Margin',4, ...
        'FontWeight','bold', 'FontSize',10); % smaller font here
    t.FontName = 'Montserrat';
end

