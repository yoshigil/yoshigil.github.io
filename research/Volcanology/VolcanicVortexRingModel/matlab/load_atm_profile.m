function profile = load_atm_profile(filepath, varargin)
% LOAD_ATM_PROFILE  Load a user-supplied atmospheric sounding (e.g.
% radiosonde or reanalysis profile) for use with atmosphere_TP.m.
%
%   profile = load_atm_profile(filepath)
%   profile = load_atm_profile(filepath,'TempUnits','C','PressUnits','hPa')
%
% The file should be a delimited text table (CSV by default) with columns
% for height, temperature and pressure. Column names are auto-detected
% (case-insensitive match on 'height'/'alt'/'z', 'temp'/'t', 'press'/'p')
% unless given explicitly via the 'HeightCol','TempCol','PressCol'
% name-value pairs.
%
% Name-value options:
%   'TempUnits'  : 'K' (default) or 'C'
%   'PressUnits' : 'Pa' (default) or 'hPa'
%   'HeightCol','TempCol','PressCol' : explicit column names
%   'Delimiter'  : column delimiter (default ',')
%
% Returns a struct with fields z [m], T [K], P [Pa], sorted by ascending z.
%
% Johan Gilchrist

p = inputParser;
addParameter(p,'TempUnits','K');
addParameter(p,'PressUnits','Pa');
addParameter(p,'HeightCol','');
addParameter(p,'TempCol','');
addParameter(p,'PressCol','');
addParameter(p,'Delimiter',',');
parse(p,varargin{:});
opts = p.Results;

Tbl = readtable(filepath,'Delimiter',opts.Delimiter);
vn = lower(string(Tbl.Properties.VariableNames));

iz = find_col(vn, opts.HeightCol, ["height","alt","z"]);
it = find_col(vn, opts.TempCol,  ["temp","t"]);
ip = find_col(vn, opts.PressCol, ["press","p"]);

z = Tbl{:,iz};
T = Tbl{:,it};
P = Tbl{:,ip};

if strcmpi(opts.TempUnits,'C')
    T = T + 273.15;
end
if strcmpi(opts.PressUnits,'hPa')
    P = P*100;
end

[z,order] = sort(z);
profile.z = z;
profile.T = T(order);
profile.P = P(order);

end

function idx = find_col(varnames, explicitName, keywords)
if ~isempty(explicitName)
    idx = find(strcmpi(varnames, lower(explicitName)), 1);
    if isempty(idx)
        error('load_atm_profile:noColumn','Column ''%s'' not found.',explicitName);
    end
    return
end
for k = 1:numel(keywords)
    idx = find(contains(varnames, keywords(k)), 1);
    if ~isempty(idx); return; end
end
error('load_atm_profile:noColumn','Could not auto-detect a column matching %s.', ...
    strjoin(cellstr(keywords),', '));
end
