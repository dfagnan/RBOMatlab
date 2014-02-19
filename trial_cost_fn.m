%Determine cost of funding trials (stochastically) for phase comp.phase
%Input:
%@comp.phase is a cell of strings where each is 'DSC','PRE','P1','P2','P3','NDA', 'APP', 'SLD'
%@params is a struct containing the initial set of parameters set by the function init_params_fn
%Output:
%@cost is of type double and is the returned funding cost
function cost = trial_cost_fn(comp_phase, params, comp_types)

    cost = zeros(size(comp_phase));
    rowidxs    = phase2index_fn(comp_phase);
    elig_idx = ~( strcmp('DSC',comp_phase) | strcmp('NDA',comp_phase) | strcmp('APP',comp_phase) | strcmp('SLD',comp_phase) );
    if(~any(elig_idx))
       return;
    end
    
    mu = nan(size(comp_phase));
    sigma = nan(size(comp_phase));    
    mx = nan(size(comp_phase));    
    for i=1:length(unique(comp_types))
        mu(elig_idx) = params.assets(i).pricing_params(rowidxs(elig_idx), name2index_fn(params.assets(i), 'pricing_params', 'mu', false))';
        sigma(elig_idx) = params.assets(i).pricing_params(rowidxs(elig_idx), name2index_fn(params.assets(i), 'pricing_params', 'sigma', false));
        mx(elig_idx) = params.assets(i).pricing_params(rowidxs(elig_idx), name2index_fn(params.assets(i), 'pricing_params', 'max', false));
    end
    
    raw_cost  = lognrnd(mu',sigma');
    cost = min(raw_cost, mx')';
end
%COPYRIGHT 2012,2013
% This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
