%Input:
%@results is of type struct and is a simulation results summary object
%which is the resultant  output of the simulate_CFs_fn function 
function plot_standard_diagnostics_fn(results)
    
    SMALLFONTSZ = 10;
    BIGFONTSZ   = 12;

    subplot(3,2,1);
    plot_trajectories_fn(results.cash,true);
    ylabel('Period','FontSize',SMALLFONTSZ);
    title('Cash - Quantiles','FontSize',BIGFONTSZ);

    subplot(3,2,2);
   	ROE_raw = results.ROE;
	NYEARS  = length(results.cash(1,:))/2;
    ROE_annualized = ((1+ROE_raw).^(1/NYEARS))-1;
   
    hist(ROE_annualized,SMALLFONTSZ);
    ylabel('Frequency','FontSize',SMALLFONTSZ);
    title('ROE Annualized','FontSize',BIGFONTSZ);

    subplot(3,2,3);
    plot_trajectories_fn(results.A1_bals,true);
    ylabel('Period','FontSize',SMALLFONTSZ);
    title('A1 balance - Quantiles','FontSize',BIGFONTSZ);

    subplot(3,2,4);
    plot_trajectories_fn(results.A2_bals,true);
    ylabel('Period','FontSize',SMALLFONTSZ);
    title('A2 balance - Quantiles','FontSize',BIGFONTSZ);

    
    subplot(3,2,5);
    sale_phases = mean(results.sale_phases,1);
    withdrawal_phases = mean(results.withdrawal_phases,1);
   	col_names = fieldnames(results.withdrawal_phases_col);
   	row_names = {'Sales','WD'};
    exits     = [sale_phases; withdrawal_phases]';
    h         = bar(exits);
    positions = get(gca,'XTick');
    positions(positions>size(exits,1)) = 0;
    positions = setdiff(positions,0);
    set(gca,'XTick',positions);
    set(gca,'XTicklabel',col_names(positions));
    set(h(1),'FaceColor',[0 0.7 0.2]);
    set(h(2),'FaceColor',[1 0    0]);
    title('Mean # compounds sold (green) or withdrawn (red) in phase','FontSize',BIGFONTSZ);

   	
    subplot(3,2,6);
   	sale_times       = mean(results.sale_times,1);
    withdrawal_times = mean(results.withdrawal_times,1);
   	col_names        = fieldnames(results.withdrawal_times_col);
    exit_times       = [ sale_times; withdrawal_times]';
    h                = bar(exit_times);
    positions = get(gca,'XTick');
    positions(positions>size(exit_times,1)) = 0;
    positions = setdiff(positions,0);
    set(gca,'XTick',positions);
    set(gca,'XTicklabel',col_names(positions));
    set(h(1),'FaceColor',[0 0.7 0.2]);
    set(h(2),'FaceColor',[1 0    0]);
    title('Mean # compounds sold (green) or withdrawn (red) in period','FontSize',BIGFONTSZ);

   	%exit.times<-matrix(c(sale.times,withdrawal.times),2,length(col.names),dimnames=list(row.names,col.names),byrow=T)
   	%barplot(exit.times,space=c(0,0.5),col=c(3,2),beside=T,main="Mean # compounds sold (green) or withdrawn (red) in period")

end



function plot_trajectories_fn(paths, quant)
	if(quant)
        series = quantile(paths, [0 0.25 0.5 0.75 1]);
	else
		series = paths;
	end
    series = series';
    plot(series);
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

