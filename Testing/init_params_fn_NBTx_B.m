%This function initializes all parameters used by the simulation code
%The function contains all default values used in the code
%Input:
%Output:
%@params is of type struct and has 3 fields simu (simulation parameters), assets (asset parameters), bonds (bond parameters)
function params = init_params_fn_NBTx_B()


    % ASSET PARAMETERS
    %
    assets = struct;

    assets.rho = 0.2;
    assets.mode = 'Markov';
	assets.trans_prob = [ ...


         ... %DSC    			PRE  P1       			P2 			        P3 			        NDA 		       APP 
			1.000000000000000,  0,   0,        			0,         			0.000000000000000, 	0,         			0; ...         				%DSC
			0.155000000000000,  0.5, 0.345,    			0,         			0.000000000000000,	0,         			0; ...         				%PRE





			0.053215560704185,  0,   0.807859886086702,	0.133437226118754, 	0.005307429700172,	0.000157443915296,	0.000022453474890; ...		%P1
			0.085093613975923,  0,   0,    				0.844734609586595,	0.066750087871169,	0.002866145785922, 	0.000555542780392; ...    	%P2
			0.062545675724766,  0,   0,        			0,         			0.848420924719526,	0.068234113621235, 	0.020799285934472; ...    	%P3
			0.021663504081486,  0,   0,        			0,         			0,         			0.566729918370271, 	0.411606577548242; ...    	%NDA
			0,			        0,   0,        			0,         			0,         			0,         			1.000000000000000; ...    	%APP
            ];
    assets.trans_prob_row.DSC = 1;
    assets.trans_prob_row.PRE = 2;
    assets.trans_prob_row.P1  = 3;
    assets.trans_prob_row.P2  = 4;
    assets.trans_prob_row.P3  = 5;
    assets.trans_prob_row.NDA = 6;
    assets.trans_prob_row.APP = 7;
    assets.trans_prob_col     = assets.trans_prob_row;
     
    assets.pricing_params  =  [ ... 
	...	%    vmu     			vsigma   vmx    UpFront   			Milestone    		mu       			sigma      max		FutureCostEst



             0,      			1,         0,  	0,		          	0,     				0,       			0,          0,      0;    ...     %DSC
             2.354733621857371, 0.939,   100,  	2.506648325603443, 	1.253324162801722,  1.573797446769595,  0.785050619780075,      20,    	134;  ...     %PRE
             2.976737648081240, 0.939,   250,  	7.519944976810328, 	3.759972488405164,  2.722501027059631,  0.732114106934504,    50,    	121;  ...     %P1
             3.996945974057142, 0.939,   500,   20.053186604827545,	10.026593302413772, 3.702027419856847,  0.795995957236353,     	120,   	85;   ...     %P2
             5.802748472629424, 0.939,  1000,  	75.199449768103278,	37.599724884051639, 5.065276957812039,	0.633123896262475,   	500,    0;    ...     %P3


             7.277762431143847, 0.939,  2500,  	0,          		0,    				0,       			0,          0,      0;    ...     %NDA
             7.239151725642049, 0.939,  5000, 	0,          		0,    				0,       			0,          0,      0;    ...     %APP
            ];
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

    assets.sell_in_phase = {'APP'};
    assets.sale_time = ones(7,1).*1.5;
    assets.sale_time_col = assets.trans_prob_col;
	assets.ratio_unchanged = 1;
    assets.equity_stake    = 0.85;


    %BOND PARAMETERS
    %
    bonds = struct;
    
    bonds.capital_structure = [0.4333333 0.0666667 0.5]; %A1, A2, EQ

    bonds.capital_structure_col.A1  = 1;
    bonds.capital_structure_col.A2  = 2;
    bonds.capital_structure_col.EQ  = 3;

    bonds.nominal = NaN;
    bonds.IC_pers = 2;

    bonds.interest_coverage = [1.75, 3.5,  0];  %A1, A2, EQ
    bonds.interest_coverage_col = bonds.capital_structure_col;

    bonds.amort_timing  = [...
    ... %start, stop


         5,      8; ... %A1
         9,     12; ... %A2
    ];
    bonds.amort_timing_row.A1 = 1;
    bonds.amort_timing_row.A2 = 2;
    bonds.amort_timing_col.start = 1;
    bonds.amort_timing_col.stop  = 2;
    
    bonds.coupon = [0.05/2, 0.08/2, 0];   %A1, A2, EQ
    bonds.coupon_col = bonds.capital_structure_col;

    bonds.servicing_rate = 0.0025;
    bonds.svc_accrual_int_rate = 0.07/2;
    bonds.cash_accrual_int_rate = 0.01/2;
    
    % SIMU PARAMETERS
    simu = struct;
    simu.NSIMUS = 1000;
    simu.TIMESTEPS =  max(reshape(bonds.amort_timing,1,[]))+1;
    simu.initial_cash  = 15000;
    
    simu.initial_compounds = [0, 0, 0, 100, 0, 0, 0];
    simu.initial_compounds_col = assets.trans_prob_col;

    %update any dependencies
    bonds.nominal = simu.initial_cash*bonds.capital_structure;
    bonds.nominal_col = bonds.capital_structure_col;

	%CE PARAMETERS - Guarantees
    ce.ce_amt = 0;
    ce.ce_repay = false;
    ce.ce_premium = 0;

    %create global parameter list
    params = struct;
    params.simu   = simu;
    params.assets = assets;
    params.bonds  = bonds;
	params.ce     = ce;

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
