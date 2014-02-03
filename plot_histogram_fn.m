%Input:
%@results is of type struct and is a simulation results summary object
%which is the resultant  output of the simulate_CFs_fn function 
function plot_histogram_fn(results,tit)
    
    SMALLFONTSZ = 10;
    BIGFONTSZ   = 12;
    
   	ROE_raw = results.ROE;
	NYEARS  = length(results.cash(1,:))/2;
    ROE_annualized = ((1+ROE_raw).^(1/NYEARS))-1;
    mean(ROE_annualized)
    hist(ROE_annualized,SMALLFONTSZ);
    ylabel('Frequency','FontSize',SMALLFONTSZ);
    title(tit,'FontSize',BIGFONTSZ);

end






%COPYRIGHT 2012
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

