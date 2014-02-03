%Input:
%@results is of type struct and is a simulation results summary object
%which is the resultant  output of the simulate_CFs_fn function 
%@plot is of type boolean and if true displays a graphical output of some
%summary statistics
function summarize_results_fn (results,plot)

    if ~exist('plot','var')
        plot = true;
    end

    HORIZON    = length(results.cash(1,:));
	NGUARPERS  = size(results.guarantee_draw,2);
    NYEARS     = HORIZON/2;
    NCOMPOUNDS = sum(results.params.simu.initial_compounds);
    NSIMUS     = results.params.simu.NSIMUS;
    sale_phases       = mean(results.sale_phases,1);
    withdrawal_phases = mean(results.withdrawal_phases,1);
    sale_times        = mean(results.sale_times,1);
    withdrawal_times  = mean(results.withdrawal_times,1);

    exits             = [sale_phases; withdrawal_phases];
    exits_row         = {'Sales:','WD:'};
    exits_col         = results.withdrawal_phases_col;
    exit_times        = [sale_times; withdrawal_times];
    exit_times_row    = exits_row;
    exit_times_col    = results.withdrawal_times_col;

    funds_needed = mean(results.compounds_to_fund,1);
    funds_done   = mean(results.compounds_funded,1);
    funding      = [funds_needed; funds_done];
    funding_row  = {'Funds needed:','Funds invested:'};
    funding_col  = exit_times_col;
    
    A1_gt_0      = results.A1_bals(:,results.params.bonds.amort_timing(1,end)+1)>0;
    A2_gt_0      = results.A2_bals(:,results.params.bonds.amort_timing(2,end)+1)>0;

    p1           = mean(A1_gt_0);
    p2           = mean(A2_gt_0);
    el1          = (mean(results.A1_bals(:,results.params.bonds.amort_timing(1,end)+1))./results.params.bonds.nominal(results.params.bonds.nominal_col.A1));
    el2          = (mean(results.A2_bals(:,results.params.bonds.amort_timing(2,end)+1))./results.params.bonds.nominal(results.params.bonds.nominal_col.A2));
    loss_stats     = [p1, p2; el1, el2];
    loss_stats_row = makeVectorLabels({'PD','EL'});
    loss_stats_col = makeVectorLabels({'A1','A2'});

    ROE_raw        = results.ROE;
    ROE_annualized = ((1+ROE_raw).^(1/NYEARS))-1;
    ROE_mean       = mean(ROE_raw);
    ROE_mean_a     = mean(ROE_annualized);
    ROE_std_a      = std(ROE_annualized)
    ROE            = [ROE_mean ROE_mean_a];
    ROE_row        = {'E(ROE)'};
    ROE_col        = {'TOT','ANN'};

    ROE_quant      = quantile(ROE_raw',[0 0.25 0.5 0.75 1]);
    ROE_quant_a    = quantile(ROE_annualized',[0 0.25 0.5 0.75 1]);
    ROE_q          = [ROE_quant ROE_quant_a]';
    ROE_q_row      = makeVectorLabels({'TOT','ANN'});
    ROE_q_col      = makeVectorLabels({'q0','q25','median','q75','q100'});
    
    p_EQ_wipeout   = mean(ROE_annualized==-1);
    p_EQ_loss      = mean(ROE_annualized<0);
    p_EQ_pos       = mean(ROE_annualized>=0);
    p_EQ_05        = mean(ROE_annualized>0.05);
    p_EQ_10        = mean(ROE_annualized>=0.1);
    p_EQ_15        = mean(ROE_annualized>0.15);
    p_EQ_25        = mean(ROE_annualized>0.25);
    
    EQ_probs       = [p_EQ_wipeout,p_EQ_loss,p_EQ_pos,p_EQ_05,p_EQ_10,p_EQ_15,p_EQ_25];
    EQ_probs_row   = {'p(x)'};
    EQ_probs_col   = {'EQ=0','EQ_loss','EQ_pos','EQ>0_05','EQ>0_10','EQ>0_15','EQ>0_25'};
        
	disc =  repmat((1+0.02/2).^(1:NGUARPERS),NSIMUS,1);
    cost_guar = sum(results.guarantee_draw./disc,2);
    
    GUAR_quant     = quantile(cost_guar,[0 0.01 0.02 0.05 0.10])';
    GUAR_quant_row = {'cost guarantee 2% dr'};
    GUAR_quant_col = makeVectorLabels({'q0','q01','q02','q05','q10'});
    
    mean_cost_guar = mean(cost_guar,1);
    prob_guar = sum(results.drew_on_guarantee)/NSIMUS;
    fprintf(1,'\n\nv==================================v\n');
    fprintf(1,'\n# Simu: %d Horizon: %4.1f years_\n',results.params.simu.NSIMUS,NYEARS);
    fprintf(1,'\n--------- ASSET  ANALYSIS ----------\n\n');

    fprintf(1,'Initial target portfoilio contains %d compounds:\n',NCOMPOUNDS);
    printvec(results.params.simu.initial_compounds,results.params.simu.initial_compounds_col);

    fprintf(1,'\nMean number of compounds initially purchased:\n');
    printvec(mean(results.compounds_initially_bought,1),results.compounds_initially_bought_col);

    fprintf(1,'\nMean number of compounds exiting in each state:\n');
    printmat(exits,exits_col, exits_row,3);

    fprintf(1,'\nMean number of compounds exiting in each period:\n');
    printmat(exit_times, exit_times_col, exit_times_row,3);

    fprintf(1,'\nMean number of compounds funded in each period:\n');
    printmat(funding, funding_col, funding_row,3)  ;
    
    fprintf(1,'\n--------- EQUITY ANALYSIS ----------\n\n');
    printmat(ROE, ROE_col, ROE_row,5);

    fprintf(1,'\n %15s \n','quant(ROE)');
    printmat(ROE_q, ROE_q_col, ROE_q_row, 5);

    fprintf(1,'\n');
    printmat(EQ_probs, EQ_probs_col, EQ_probs_row, 5);

    fprintf(1,'\n---------  BOND ANALYSIS  ----------\n\n');
    printmat(loss_stats, loss_stats_col, loss_stats_row, 5);
    
	   fprintf(1,'\n---------  GUARANTEE ANALYSIS  ----------\n\n');
    printmat(GUAR_quant,GUAR_quant_col,GUAR_quant_row,3);
    fprintf(1,'\nMean amount of guarantee drawn: %f\n',mean_cost_guar);
    fprintf(1,'\nProbability of guarantee drawn: %f\n',prob_guar);
    fprintf(1,'\n\n^==================================^\n\n');
    
    
    
    if(plot) 
        plot_standard_diagnostics_fn(results);
    end


end %function summarize_results_fn


function printvec(field, rowcol, prec)
    
    if ~exist('prec','var')
        prec = 2;
    end
    precstr = sprintf(' %%8.%df ',prec);

    if iscell(rowcol)
        names = rowcol;
    else
        names = fieldnames(rowcol);
    end
    fprintf(1,' %15s ','');
    fprintf(1,' %8s ',names{:});
    fprintf(1,'\n');
    fprintf(1,' %15s ','');
    fprintf(1,precstr,field);
    fprintf(1,'\n');

end

function printmat(field, cols, rows, prec)

    if ~exist('prec','var')
        prec = 2;
    end
    precstr = sprintf(' %%8.%df ',prec);

    if iscell(cols)
        names = cols;
    else
        names = fieldnames(cols);
    end
    fprintf(1,' %15s ','');
    fprintf(1,' %8s ',names{:});
    fprintf(1,'\n');


    if iscell(rows)
        names = rows;
    else
        names = fieldnames(rows);
    end
    
    for ctr = 1 : size(field,1)
        fprintf(1,' %15s ',names{ctr});
        fprintf(1,precstr,field(ctr,:));
        fprintf(1,'\n');
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
