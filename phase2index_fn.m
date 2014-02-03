%Convert character string representation of phase into an index
%Input:
%@phase is either a cell of strings or a single string
%Output:
%@idx is of type integer and maps the string phase to a corresponding index
function indxs = phase2index_fn(phase)
    
    names = {'DSC','PRE','P1','P2','P3','NDA', 'APP'};
    if(iscell(phase)) 
        indxs = zeros(size(phase));
        for i = 1:length(names)
            indxs = indxs + i*strcmp(names{i},phase);
        end
    else 
        indxs = find(strcmp(phase,names));
    end
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
