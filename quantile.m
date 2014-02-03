%A simple quantile function
function qs = quantile(vec, percs)


    qs = nan(length(percs),size(vec,2));
    for ctr = 1 : size(vec,2)
        
        cvec = vec(:,ctr);
        cvec = cvec(~isnan(cvec));
        
        if ~isempty(cvec)
            p  = floor(size(cvec,1).*percs);
            p(p==0) = 1;

            [vals, pos] = sort(cvec);

            qs(:,ctr)   = vals(p);
        end
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
