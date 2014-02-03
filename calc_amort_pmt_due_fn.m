%Input:
%@per is an integer
%@a.start, and @a.end are integers representing period of start (end) of amortization
%@orig.par and @cur.par are doubles representing par at origination and outstanding
%Output:
%@pmt is a double represting the calculated principal payment due in period per
function pmt = calc_amort_pmt_due_fn(per, a_start, a_end, orig_par, cur_par)
    pmt       = 0;
    amort_per = a_end - a_start + 1;

    if (per>=a_start && per<=a_end)
        if (orig_par/amort_per < cur_par)
            pmt = orig_par/amort_per;
        else
            pmt = cur_par;
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
