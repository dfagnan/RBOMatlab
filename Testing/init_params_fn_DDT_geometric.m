%This function initializes all parameters used by the simulation code
%This version uses the orphan drug parameters discussed in 
%"Financing drug discovery for orphan diseases" Fagnan et al. 2013.
%@params is of type struct and has 4 fields simu (simulation parameters), 
%   assets (asset parameters), bonds (bond parameters), and ce (guarantee parameters)
function params = init_params_fn_DDT_geometric()

    % ASSET PARAMETERS
    assets = struct;
	assets.mode = 'General';
    assets.distribution = 'Geometric';
	%Correlation value for asset valuation.
    assets.rho = 0.2;
   
	%Present value upon succesful FDA approval.
	approvalValue = 817.643;
	
	%Costs for clinical trials by phase.
	costs = [5,5,8,43,0];
	
	%Discount rates for estimating valuation at phases prior to approval.
	disc = [0.3, 0.3, 0.3, 0.25, 0.15];
	
	%Duration of each phase in years.
	dur = [1, 1.6561, 2.0942, 2.1497, 0.8];
	
	%Probabilities of technical and regulatory success by phase.
	%%Note that the original paper by Fernandez et al. uses a more
	%%advanced generating matrix approach using data unavailable for
	%%the Orphan drug paper published in Drug Discovery Today.
	steady = [.69, .84, .53, .74, .96];
	assets.success = steady;
    assets.dur_mean = dur;
	%Target phase (eg. 4 = P3, 3 = P2), used for calculating budgeting rules.
   	target = 4;
	
	%Orphan drug utility functions build the capped log-normal distributions
	%Uses values for the standard deviations derived from previous work by Fernandez et al. 2012.
    assets.pricing_params = createTransMat(approvalValue,costs,disc,dur,steady,target);
	%%assets.trans_prob=createProbMat(steady,dur);
   
   %%Set the named variables to the corresponding column numbers.
   %DSC refers to a state of withdrawal or discontinued.
    assets.trans_prob_row.DSC = 1;
    assets.trans_prob_row.PRE = 2;
    assets.trans_prob_row.P1  = 3;
    assets.trans_prob_row.P2  = 4;
    assets.trans_prob_row.P3  = 5;
    assets.trans_prob_row.NDA = 6;
    assets.trans_prob_row.APP = 7;
    assets.trans_prob_col     = assets.trans_prob_row;
    
    assets.pricing_params_row     = assets.trans_prob_row;
    assets.pricing_params_col.vmu            = 1;
    assets.pricing_params_col.vsigma         = 2;
    assets.pricing_params_col.vmx            = 3;
    assets.pricing_params_col.UpFront        = 4;
    assets.pricing_params_col.Milestone      = 5;
    assets.pricing_params_col.mu             = 6;
    assets.pricing_params_col.sigma          = 7;
    assets.pricing_params_col.max            = 8;
    assets.pricing_params_col.FutureCostEst  = 9;
	
	%Sell compounds upon reaching this phase.
	%Note that it is possible for a compound to transition
	%past the target phase due to parallel trials (code
	%only checks for transition every six months).
    assets.sell_in_phase = {'P3'};
    assets.sale_time = ones(7,1).*1.5;
    assets.sale_time_col = assets.trans_prob_col;
	
	%Percent that assets are worth if they remain in the same stage
	%for the entirety of the simulation. 
	assets.ratio_unchanged = 1.00;
	
	%Share of valuation payouts that go to the megafund.
	%eg. 0.85 means 85% of the sale value goes to the megafund
	% with 15% paid out as royalties, lost to the megafund.
    assets.equity_stake    = 0.85;


    %BOND PARAMETERS
    bonds = struct;
	%Capital structure percentages e.g. senior, junior, equity
    bonds.capital_structure = [0.15 0.2 0.65]; %A1, A2, EQ

    bonds.capital_structure_col.A1  = 1;
    bonds.capital_structure_col.A2  = 2;
    bonds.capital_structure_col.EQ  = 3;

	%Number of periods for which to cover expected interest
	% and coupon payments with cash. (interest coverage)
    bonds.IC_pers = 2;

	% Interest coverage ratios
	% e.g. hold cash equal to 1.75x future senior tranche interest & coupon
	% and 3.5x the future junior tranche interest & coupon
    bonds.interest_coverage = [1.75, 2.75,  0];  %A1, A2, EQ
    bonds.interest_coverage_col = bonds.capital_structure_col;

	%Bond amortization schedule.
	%e.g. Senior tranche is paid out from Period 5 through 8
	% 	  Junior tranche is paid out from Period 9 through 12
	%Note that Coupon payments are paid semi-annually until completion of the schedule.
    bonds.amort_timing  = [...
    ... %start, stop
         5,      8; ... %A1
         9,     12; ... %A2
    ];
    bonds.amort_timing_row.A1 = 1;
    bonds.amort_timing_row.A2 = 2;
    bonds.amort_timing_col.start = 1;
    bonds.amort_timing_col.stop  = 2;
    
	%Annual coupon rates scaled for semi-annual payments.
    bonds.coupon = [0.05/2, 0.08/2, 0];   %A1, A2, EQ
    bonds.coupon_col = bonds.capital_structure_col;

	%Servicing rate or management fee for the fund.
    bonds.servicing_rate = 0.0025;
	
	%Accrual rate on missed service payments.
    bonds.svc_accrual_int_rate = 0.07/2;
	
	%Interest rate on unused cash on hand.
    bonds.cash_accrual_int_rate = 0.01/2;
    
    %CE PARAMETERS - Guarantees
    ce.ce_amt = 0;
    ce.ce_repay = false;
    ce.ce_premium = 0;
    
    
    % SIMU PARAMETERS
    simu = struct;
	%Number of simulations recommend 100K for quick estimates,
	% or 2 million for a complete profile.
    simu.NSIMUS = 100000;
    simu.TIMESTEPS =  max(reshape(bonds.amort_timing,1,[]))+1;
	%Total capital to be spread across tranches.
    simu.initial_cash  = 575;
    
    simu.initial_compounds = [0, 8, 8, 0, 0, 0, 0];
    simu.initial_compounds_col = assets.trans_prob_col;
    
    %Set the nominal amounts corresponding to the CS and total capital.
    bonds.nominal = simu.initial_cash*bonds.capital_structure;
    bonds.nominal_col = bonds.capital_structure_col;

    %create global parameter list
    params = struct;
    params.simu   = simu;
    params.assets = assets;
    params.bonds  = bonds;
    params.ce     = ce;
end
%COPYRIGHT 2013
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
