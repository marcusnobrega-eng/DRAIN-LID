function [params, clog] = configure_clogging(params, opt)
% CONFIGURE_CLOGGING
% Build clogging state/config from user options and current params.
%
% Notes
% -----
% • This function does NOT change params.theta_s or params.Ks. It only
%   stores the "clean" reference copies (theta_s0, Ks0) and prepares
%   state variables used by the step-wise updater.
% • Call this ONCE during setup (after you assign params.theta_s / params.Ks).
%
% Required fields in params:
%   params.Nz  : number of vertical nodes
%   params.z   : node elevations [bottom..top] (length Nz)
%   params.theta_r, params.theta_s, params.Ks : vectors length Nz

    % --------- Defensive checks ----------
    assert(isfield(params,'Nz') && isscalar(params.Nz), ...
        'configure_clogging: params.Nz must exist and be a scalar.');
    assert(isfield(params,'z') && numel(params.z)==params.Nz, ...
        'configure_clogging: params.z must be a vector of length Nz.');
    mustHave = {'theta_r','theta_s','Ks'};
    for f = mustHave
        assert(isfield(params, f{1}) && numel(params.(f{1}))==params.Nz, ...
            'configure_clogging: params.%s must be a vector of length Nz.', f{1});
    end
    Nz = params.Nz;

    % --------- User knobs with safe defaults ----------
    % gamma:     exponential decay rate of porosity vs cumulative loading Ic [1/m]
    % mK:        porosity–permeability exponent (Ks ∝ (theta_s/theta_s0)^mK) [-]
    % phi_min:   irrecoverable porosity floor as a fraction of clean theta_s0 [-]
    % eta:       maintenance efficacy (fraction of recoverable gap restored) [-]
    clog.gamma        = getfielddef(opt, 'gamma',        2e-3);   % [1/m]
    clog.mK           = getfielddef(opt, 'mK',           6.0);    % [-]
    clog.phi_min_frac = getfielddef(opt, 'phi_min_frac', 0.10);   % [-]
    clog.eta          = getfielddef(opt, 'eta',          0.7);    % [-]

    % Maintenance policy: EITHER a fixed period (days) OR explicit list (seconds)
    clog.maintenance_period_days = getfielddef(opt, 'maintenance_period_days', []);
    clog.maintenance_dates_sec   = getfielddef(opt, 'maintenance_dates_sec',   []);

    % --------- Where clogging applies (mask) ----------
    % Option A: top_N_cells
    % Option B: depth_above_top_m (window measured downward from the top node)
    % Option C: explicit logical mask (vector length Nz)
    mask = false(1,Nz);

    if isfield(opt, 'top_N_cells') && ~isempty(opt.top_N_cells)
        n = max(1, min(Nz, round(opt.top_N_cells)));
        mask(end-n+1:end) = true;                % mark top n cells
    elseif isfield(opt, 'depth_above_top_m') && ~isempty(opt.depth_above_top_m)
        ztop = params.z(end);                    % top node elevation
        mask(params.z >= (ztop - opt.depth_above_top_m)) = true;
    elseif isfield(opt, 'mask') && numel(opt.mask) == Nz
        mask = logical(opt.mask(:));
    else
        % Sensible default: top 2 cells
        mask(max(1, Nz-1):Nz) = true;
    end
    clog.mask = mask;

    % --------- Immutable "clean" references ----------
    % These are never overwritten. Maintenance and updates always refer to them.
    clog.theta_s0 = params.theta_s;          % clean (reference) theta_s
    clog.Ks0      = params.Ks;               % clean (reference) Ks

    % φ_min per cell (respect theta_r for physical bounds)
    % φ_min = max(phi_min_frac * theta_s0, theta_r + eps)
    clog.phi_min  = max(clog.phi_min_frac * clog.theta_s0, params.theta_r + 1e-6);

    % --------- State variables ----------
    % Ic: cumulative downward loading since last maintenance [m]
    % next_maint_t: next scheduled maintenance time [s] since t=0
    clog.Ic           = zeros(1,Nz);
    clog.last_maint_t = 0;
    clog.next_maint_t = [];
    if ~isempty(clog.maintenance_period_days)
        clog.next_maint_t = clog.maintenance_period_days * 86400; % first maintenance
    elseif ~isempty(clog.maintenance_dates_sec)
        clog.next_maint_t = clog.maintenance_dates_sec(1);
    end

    % (Optional) make it easy for downstream code to know it's configured
    clog.enabled = true;
end

% ---------- tiny helper ----------
function v = getfielddef(s, name, default)
% GETFIELDDEF: return s.(name) if it exists and is non-empty; otherwise default.
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end
