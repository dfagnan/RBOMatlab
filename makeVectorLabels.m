function labelstruct = makeVectorLabels(lab, addctr)

    if ~exist('addctr','var')
        addctr = false;
    end

    labelstruct = struct;

    for ctr = 1 : length(lab)
        if ~addctr
            labelstruct.(lab{ctr}) = ctr;
        else
            nm  = sprintf('%s%d',lab{ctr},ctr);
            labelstruct.(nm) = ctr;
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
