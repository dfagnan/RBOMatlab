%This function assumes every matrix of name say 'mat' within a structure parentstruct
%has fields mat_col and mat_row which provides a mapping from name to row/col index
%Input:
%@parentstruct is of type struct
%@fieldname is of type string and parentstruct.(fieldname) should exist
%@rowcolname is of type string and is name of a row/column for the matrix parentstruct.(fieldname)
%@isrow is boolean and is used to decide whether we search for the index having name 'rowcolname' in the structure parentstruct.(fieldname_row) or parentstruct.(fieldname_col)
%Output:
%@idx is of type int and provides a mmaping from string (rowcolname) to integer (idx) 
function idx = name2index_fn(parentstruct, fieldname, rowcolname, isrow)

    if isrow
        idxfieldname = sprintf('%s_row',fieldname);
    else
        idxfieldname = sprintf('%s_col',fieldname);
    end
    
    idx = parentstruct.(idxfieldname).(rowcolname);

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
