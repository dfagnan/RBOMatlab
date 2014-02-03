%Input:
%@p is of a vector of random draws from the uniform distribution across the open interval (0,1)
%@dist is a matrix where each row is a vector of type double and is row from the transition probability matrix params.assets.trans_prob
%@dist_col is a struct defined by params.assets.trans_prob_col
%Output:
%value  is of type String and maps the selected element in dist to its corresponding name from dist_col
function value = determine_current_state_fn(p, dist, dist_col)
    
    cdist = cumsum(dist,2);
    len   = size(cdist,2);
    cdist(:,len) = 1;
    i  = zeros(size(p));
    %Find the column that has cumulative sum of probability > p.
    for j = 1:len
        i = i +j*(cdist(:,j)'>p).*(i<=0);
    end
    %Convert the column index to the column name.
    value = cell(size(p));
    names = fieldnames(dist_col);
    for ctr = 1 : length(names)
        [value{dist_col.(names{ctr})==i}] = deal(names{ctr});
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
