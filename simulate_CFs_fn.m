%This code has been prepared to support the simulations presented in the paper "Can financial engineering cure cancer?" (David Fagnan, Jose-Maria Fernandez, Roger M. Stein, Andrew W.Lo; American Economic Association; 2013). We want to thank Lloyd Han, James Noraky, Allister Bernard and Ashutosh Singhal for their excellent coding support work. 
%Input:
%@show_progress is of type boolean and if true displays a status bar of the current simulation
%@params is of type struct and is parameter object which is the output of the init_params_fn function (optional)
%Output:
%@results is of type struct and provides details of each simulation run
function results = simulate_CFs_fn(show_progress, params)


function next_phase = check_for_transitions_ffn(compounds, invested, rho, params)
    next_phase = compounds;
    phase_idxs  = phase2index_fn(compounds);
    eligible_idx = ~strcmp('DSC',compounds) & ~strcmp('APP',compounds) & ~strcmp('SLD',compounds) & invested;
    
    p = rand(size(find(eligible_idx)));
    
    elig_phase = determine_current_state_fn(p, params.assets.trans_prob(phase_idxs(eligible_idx),:), params.assets.trans_prob_col);
    
    next_phase(eligible_idx)=elig_phase; 
end


function next_phase = check_for_transitions_ffn_distribution(compounds, invested, time_required, success, params)
    next_phase = compounds;
    phase_idxs  = phase2index_fn(compounds);
    eligible_idx = ~strcmp('DSC',compounds) & ~strcmp('APP',compounds) & ~strcmp('SLD',compounds) & invested ...
    & ((time_in_phase+1)>= time_required(sub2ind(size(time_required),[1:NCOMPOUNDS],max(phase_idxs,1))));
    
    success_idx = (rand(size(compounds)))<success(max(phase_idxs,1));
    phase_idxs(eligible_idx & success_idx) = phase_idxs(eligible_idx & success_idx) + 1;
    phase_idxs(eligible_idx & ~success_idx) = 1;
    names = {'DSC','PRE','P1','P2','P3','NDA', 'APP'};
    next_phase(eligible_idx)=names(phase_idxs(eligible_idx)); 
end

%Nested Function
%Changes guarantee_draw,drew_on_guarantee
function [cur_cash, cur_guarantee] = repay_guarantee_ffn(cur_per,cur_cash,cur_guarantee)
	%Experimental function, not currently used.
    ce_pmt_due = params.ce.ce_amt-cur_guarantee;
    ce_pmt = min(cur_cash,ce_pmt_due);
    cur_cash=cur_cash-ce_pmt;
    cur_guarantee = cur_guarantee+ce_pmt;
end

%Nested Function
%Changes guarantee_draw,drew_on_guarantee
function [cur_cash, cur_guarantee] = draw_from_guarantee_ffn(cur_per,cur_cash,cur_guarantee)
    draw_amt = cur_cash;
    guarantee_draw(s,cur_per) = guarantee_draw(s,cur_per) + draw_amt;
    cur_cash = 0;
    cur_guarantee = cur_guarantee + draw_amt;
    drew_on_guarantee(s) = true; 
end

%Nested Function
%Accesses params
function  svc_due_this_per = calc_svc_due_this_per_ffn(values)    
	%Use of ce.premium is experimental and untested.
    svc_due_this_per = params.bonds.servicing_rate.*sum(values)+params.ce.ce_amt*params.ce.ce_premium;
end

%Nested function
%Accesses params,bond_value
%Changes A1_in_default, A2_in_default
function default = check_for_default_ffn(cur_per,cur_cash,cur_guarantee)
    svc_due  = calc_svc_due_this_per_ffn(bond_value) + unpaid_svc*(1+params.bonds.svc_accrual_int_rate);
    A1_cash_required = AMORT_SCHED(1,cur_per)+bond_value(1)*bonds.coupon(1);
    A2_cash_required = AMORT_SCHED(2,cur_per)+bond_value(2)*bonds.coupon(2);
    if(cur_cash+cur_guarantee<svc_due)
        svc_in_default = true;
        if(bond_value(1)>0) 
            A1_in_default = true;
        end
        if(bond_value(2)>0)
            A2_in_default = true;
        end
    elseif(cur_cash+cur_guarantee<svc_due+A1_cash_required)
        A1_in_default = true;
        if(bond_value(2)>0)
            A2_in_default = true;
        end
    elseif((cur_cash+cur_guarantee<svc_due+A1_cash_required+A2_cash_required) && (bond_value(2)>0))
        A2_in_default = true;
    end
    default = svc_in_default || A1_in_default || A2_in_default;
    
end

%Nested function
%Changes interest_paid, drew_on_guarantee, guarantee_draw
%Accesses params
function [cur_cash,cur_guarantee] = pay_principal_ffn(cur_per,cur_cash,cur_guarantee)
   %PAY P&I
   for(b = 1:NBONDS) %pay principal
       P = AMORT_SCHED(b,cur_per);
       principal_paid(b,cur_per) = min(P,cur_cash+cur_guarantee);
       cur_cash        = cur_cash-principal_paid(b,cur_per);
       bond_value(b)       = bond_value(b)-principal_paid(b,cur_per);
       if(bond_value(b)<0.0001)
           bond_value(b)   = 0;
       end
       
	   if(cur_cash<0)
           [cur_cash cur_guarantee] = draw_from_guarantee_ffn(cur_per,cur_cash, cur_guarantee);
       end
   end
   
end

%Nested function
%Changes interest_paid, drew_on_guarantee, guarantee_draw
%Accesses params
function [cur_cash,cur_guarantee] = pay_interest_ffn(cur_per,cur_cash,cur_guarantee)
    for(b = 1:NBONDS) % pay int
        I = bonds.coupon(b)*bond_value(b);
		interest_paid(b,cur_per) = min(I,cur_cash+cur_guarantee); 
		cur_cash       = cur_cash-interest_paid(b,cur_per);
		if(cur_cash<0)
            [cur_cash cur_guarantee] = draw_from_guarantee_ffn(cur_per,cur_cash, cur_guarantee);
        end
    end			 
end

%Nested function
%Changes unpaid_svc
%Accesses params, in_default, bond_value, unpaid_svc, history_service, unpaid_svc
function [cur_cash , cur_guarantee] = pay_servicing_ffn(cur_per,cur_cash,cur_guarantee)
    svc_OK    = true;
    if (~in_default)
        svc_due_this_per  = calc_svc_due_this_per_ffn(bond_value);
        svc_past_due      = unpaid_svc*(1+params.bonds.svc_accrual_int_rate);
        svc_due           = svc_due_this_per + svc_past_due;
        svc_paid          = min(svc_due,cur_cash+cur_guarantee);

        %update main function variables
        history_service(s) =  history_service(s) + svc_paid;
        unpaid_svc         =  round((svc_due-svc_paid)*10^5)/10^5;
        cur_cash       	   =  cur_cash-svc_paid;
		if(cur_cash<0)
            [cur_cash , cur_guarantee] = draw_from_guarantee_ffn(cur_per,cur_cash,cur_guarantee);
        end
    end
 
end

%Nested function
%Accesses NPERS, unpaid_svc, bond_value, params, AMORT_SCHED
function [IC_ratio cash_due K ] = ic_test_result_ffn(b, i, cur_cash)
    if (i ~=NPERS)
        sched_svc = unpaid_svc;
        bv = bond_value(b);
        K  = cur_cash;
        %For bonds after the first, must subtract cash paid for previous
        %bonds.
        if (b>1)
            for(bb = 1:(b-1))
                btemp = bond_value(bb);
                N = min(params.bonds.IC_pers,NPERS-i);
                for(m = 1:N)
                    K = K-AMORT_SCHED(bb,i+m)-params.bonds.coupon(bb)*btemp-params.bonds.servicing_rate*btemp;
                    btemp = btemp-AMORT_SCHED(bb,i+m);
                end
            end
        end
        sched_P  =  0;
        sched_I  =  0;
       
        for(m = 1:min(params.bonds.IC_pers,NPERS-i))
            sched_P =  sched_P + AMORT_SCHED(b,i+m);
            sched_I =  sched_I + params.bonds.coupon(b)*bv;
            sched_svc = sched_svc+params.bonds.servicing_rate*bv;
            bv = bv-AMORT_SCHED(b,i+m);		
        end
        cash_due = sched_P+sched_I+sched_svc;
        IC_ratio = K/cash_due;
    else
        IC_ratio = params.bonds.interest_coverage(b);
        cash_due = 0;
        K = cur_cash;
    end

end

%Nested function
%Changes compounds,cash,sales,sale_times, sell_value or extra_cash_from_coverage,time_in_phase
function update_sale(price,remaining_indxs,cur_per,forCoverage)
    phase_idxs = phase2index_fn(compounds(remaining_indxs));
    sale_per  = min(cur_per+ceil(params.assets.sale_time(phase_idxs)), HORIZON);
    if(~forCoverage)
        sell_value(remaining_indxs) = price;
    end
    [unique_phase,~,temp_indx] = unique(phase_idxs);
    unique_counts = accumarray(temp_indx(:),1,[],@sum)'; 
	
	%update main function vars
    sales(s,unique_phase) = sales(s,unique_phase) + unique_counts;
    sale_times(s,cur_per) = sale_times(s,cur_per)+sum(unique_counts);
    [unique_per,~,~] = unique(sale_per);
    for loop = unique_per
        cash(loop) = cash(loop) + sum(price(loop==sale_per));
        if(forCoverage)
            extra_cash_from_coverage(s,loop) = extra_cash_from_coverage(s,loop)+sum(price(loop==sale_per)); 
        end
    end
    [compounds{remaining_indxs}] =  deal('SLD');
    time_in_phase(remaining_indxs) = 0;
end

%Nested Function
%Accesses NCOMPOUNDS, compounds, params, HORIZON, rho, z, equity_in_compound, invested, invest_cost, time_in_phase, sell_value, sales, sale_times, cash
function liquidate_portfolio_ffn(cur_per)
    remaining_indxs = ~strcmp('DSC',compounds) & ~strcmp('SLD',compounds);
    if(~any(remaining_indxs))
        return
    end

    phase_idxs = phase2index_fn( compounds(remaining_indxs));
    mx       = params.assets.pricing_params(phase_idxs, params.assets.pricing_params_col.('vmx'));
    mu       = params.assets.pricing_params(phase_idxs, params.assets.pricing_params_col.('vmu'));
    sigma    = params.assets.pricing_params(phase_idxs, params.assets.pricing_params_col.('vsigma'));
    price    = compound_sale_price_fn(mu, sigma, mx, rho, z).*equity_in_compound(remaining_indxs);
    
    % here is where we adjust for uninvested compounds
    rem_inv_cost = invest_cost(remaining_indxs);
    rem_not_invested = ~invested(remaining_indxs);
    price(rem_not_invested) = max(price(rem_not_invested)-rem_inv_cost(rem_not_invested),0);
    
    %adjust for compounds that have not transitioned at all until liquidation
    rem_not_transitioned = time_in_phase(remaining_indxs)==cur_per-1;
    price(rem_not_transitioned) = price(rem_not_transitioned)*params.assets.ratio_unchanged;

    % update bookkeeping AND cash AND compounds.
    update_sale(price,remaining_indxs,cur_per,false)
end

%Nested Function
%Accesses params, compounds, HORIZON, equity_in_compound, invested, invest_cost, s, cash, time_in_phase, extra_cash_from_coverage, sales, sale_times
function sell_compounds_to_cover_shortfall_ffn(shortfall,cur_per)
    %  We simply sell the most mature compounds first
    remaining_comp_idx = ~strcmp(compounds,'DSC') & ~strcmp(compounds,'SLD');
    if(~any(remaining_comp_idx))
        return
    end
    
    phase_idxs 	= phase2index_fn(compounds(remaining_comp_idx));
    mx       	= params.assets.pricing_params(phase_idxs, params.assets.pricing_params_col.('vmx'));
    mu       	= params.assets.pricing_params(phase_idxs, params.assets.pricing_params_col.('vmu'));
    sigma   	= params.assets.pricing_params(phase_idxs, params.assets.pricing_params_col.('vsigma'));
    prices   	= compound_sale_price_fn(mu, sigma, mx, rho, z).*equity_in_compound(remaining_comp_idx);

    % here is where we adjust for uninvested compounds
    rem_inv_cost 			= invest_cost(remaining_comp_idx);
    rem_not_invested		= ~invested(remaining_comp_idx);   
    prices(rem_not_invested)= max(prices(rem_not_invested)-rem_inv_cost(rem_not_invested),0);
	
    [~, sorted_idx] 		= sort(prices,'descend');
    cumSale 				= cumsum(prices(sorted_idx));
    numSold 				= find(cumSale>=shortfall,1);
    if(isempty(numSold))
        numSold = length(prices);
    end
    indxSold 				= sorted_idx(1:numSold);
    remaining_comp_loc 		= find(remaining_comp_idx);
    
	% update bookkeeping AND cash AND compounds.
    
    update_sale(prices(indxSold),remaining_comp_loc(indxSold),cur_per,true);     
end


%This needs to be a nested-function within simulate_CFs_fn
%Accesses 
%compounds, invested, params, sales, sale_times, cash, time_in_phase, HORIZON, equity_in_compound, invest_cost, sell_value, s, z
function prices = sell_compound_ffn(comp_idx, cur_per)
    remaining_comp_idx = comp_idx & ~strcmp('DSC',compounds) & ~strcmp('SLD',compounds);
    if(~any(remaining_comp_idx))
        prices = [];
        return
    end

    phase_idxs  = phase2index_fn( compounds(remaining_comp_idx));
    mx       	= params.assets.pricing_params(phase_idxs, params.assets.pricing_params_col.('vmx'));
    mu       	= params.assets.pricing_params(phase_idxs, params.assets.pricing_params_col.('vmu'));
    sigma    	= params.assets.pricing_params(phase_idxs, params.assets.pricing_params_col.('vsigma'));
    prices    	= compound_sale_price_fn(mu, sigma, mx, rho, z).*equity_in_compound(remaining_comp_idx);
    
    % here is where we adjust for uninvested compounds
    rem_inv_cost 	 = invest_cost(remaining_comp_idx);
    rem_not_invested = ~invested(remaining_comp_idx);   
    prices(rem_not_invested) = max(prices(rem_not_invested)-rem_inv_cost(rem_not_invested),0);
    
	% update bookkeeping AND cash AND compounds.
    update_sale(prices,remaining_comp_idx,cur_per,false);
end

%Nested function
%Changes withdrawals, withdrawal_times
function withdraw_compound_ffn(old_phase,cur_per)
    [uniq_tmp,~,unique_idx] = unique(old_phase);
    uniq_cnt = accumarray(unique_idx(:),1,[],@sum)'; 
    for lp = 1:length(uniq_tmp)
        withdrawals(s,withdrawals_col.(uniq_tmp{lp}))  = withdrawals(s,withdrawals_col.(uniq_tmp{lp})) + uniq_cnt(lp);
    end
    withdrawal_times(s, cur_per)    = withdrawal_times(s,cur_per) + sum(uniq_cnt);		
end

%Nested function
%Accesses compounds
%Changes  invested, invest_cost, time_in_phase
function transition_compound_ffn(idx)
    ndaidx 					= strcmp('NDA',compounds);
    invested(ndaidx & idx)  = true;	
    invested(~ndaidx & idx) = false;
    cindx             		= phase2index_fn(compounds(~ndaidx & idx));
    invest_cost(~ndaidx & idx)   = trial_cost_fn(compounds(~ndaidx & idx), params) + params.assets.pricing_params(cindx,params.assets.pricing_params_col.('Milestone'))';
    time_in_phase(idx) 		= 0;
end

    tStart = tic;
    if ~exist('params','var')
        params = init_params_fn();
    end
    if ~exist('show_progress','var')
        show_progress = true;
    end
    
    if show_progress
        prog_bar_click = round(params.simu.NSIMUS/10);
    end


    %setup global variables and constants
    simu                  = params.simu;
    assets                = params.assets;
    bonds                 = params.bonds;
	ce                    = params.ce;
   	NSIMUS                = simu.NSIMUS;
   	NPERS                 = simu.TIMESTEPS;
   	HORIZON               = NPERS+ceil(max(assets.sale_time));
   	NCOMPOUNDS            = sum(simu.initial_compounds(:)) ;
	ALL_POSSIBLE_STATES   = fieldnames(assets.trans_prob_col); %c("DSC","PRE","P1","P2","P3","NDA","APP")
	NSTATES               = size(simu.initial_compounds,2); %length(simu.initial_compounds); 
	NBONDS                = length(bonds.capital_structure)-1;
	STARTING_STATES       = {}; %repmat({'DSC'},1, simu.initial_compounds(1,1)) ;
	AMORT_SCHED           = zeros(NBONDS,NPERS);

    STARTING_PERIOD = []; 
	for t=1:size(simu.initial_compounds,1)
        for ctr = 1:NSTATES
            phase   =  ALL_POSSIBLE_STATES(ctr);
            tmpcell = repmat(phase, 1, simu.initial_compounds(t,ctr));
            STARTING_STATES = {STARTING_STATES{:} tmpcell{:} };
            STARTING_PERIOD = [STARTING_PERIOD repmat(t,1,length(tmpcell))];             
        end
    end
    
	temp_bond_value = bonds.nominal;
	for cctr = 1:NPERS
		for rctr = 1:NBONDS
			AMORT_SCHED(rctr,cctr) = calc_amort_pmt_due_fn(cctr, ...
            bonds.amort_timing(rctr,bonds.amort_timing_col.start), ...
            bonds.amort_timing(rctr,bonds.amort_timing_col.stop), ...
            bonds.nominal(rctr), temp_bond_value(rctr));
			temp_bond_value(rctr)  = temp_bond_value(rctr) - AMORT_SCHED(rctr,cctr);
		end
    end
    
    % ----------------------------------

	rho               	= assets.rho;
    if strcmp(params.assets.mode,'General')==1
       success             = [0, assets.success, 0];
       dur_mean            = assets.dur_mean;
    end
    
	bond_value        	= bonds.nominal;

	invested            = false(1,NCOMPOUNDS);
    acquired            = false(1,NCOMPOUNDS);
	sell_value          = zeros(1,NCOMPOUNDS);
	equity_in_compound  = ones(1,NCOMPOUNDS).*assets.equity_stake;
	compounds           = STARTING_STATES;
	time_in_phase       = zeros(1,NCOMPOUNDS);
	
    first_guarantee      = ce.ce_amt;
	first_cash          = simu.initial_cash;
	
	% state of securities  
	IC_ratio            = zeros(length(bonds.nominal), NPERS);
	
	%record keeping 
	sales     			= zeros(NSIMUS, NSTATES);
    sales_row 			= makeVectorLabels(repmat({'a'},1,NSIMUS),true);
    sales_col 			= makeVectorLabels(ALL_POSSIBLE_STATES);

	withdrawals       	= zeros(NSIMUS, NSTATES);
    withdrawals_row   	= sales_row;
    withdrawals_col   	= sales_col;

	sale_times      	= zeros(NSIMUS, HORIZON);
    sale_times_row   	= sales_row;
    sale_times_col   	= makeVectorLabels(repmat({'S'},1,HORIZON),true);

	withdrawal_times       		= zeros(NSIMUS, HORIZON);
    withdrawal_times_row   		= sales_row;
    withdrawal_times_col   		= sale_times_col;

	compounds_to_fund       	= zeros(NSIMUS, HORIZON);
    compounds_to_fund_row   	= sales_row;
    compounds_to_fund_col   	= sale_times_col;

	compounds_funded       		= zeros(NSIMUS, HORIZON);
    compounds_funded_row   		= sales_row;
    compounds_funded_col   		= sale_times_col;

	compounds_bought       		= zeros(NSIMUS, HORIZON, NSTATES);
    compounds_bought_row   		= sales_row;
    compounds_bought_col   		= sales_col;
    
    compounds_spent       		= zeros(NSIMUS, 1);

	final_compounds       		= zeros(1, NSTATES);
	guarantee_draw        = zeros(NSIMUS, HORIZON);
    drew_on_guarantee     = zeros(NSIMUS,1);
    
	extra_cash_from_coverage	= zeros(NSIMUS,HORIZON) ;
	cash_for_investment         = zeros(NSIMUS,NPERS);
    cash_for_buying             = zeros(NSIMUS,NPERS);
	interest_paid	            = zeros(length(bonds.nominal), NPERS);
	principal_paid              = zeros(length(bonds.nominal), NPERS);
	history_service		        = zeros(NSIMUS,1);
	cash_realized	            = zeros(NSIMUS,HORIZON);
	cash_begin_per              = zeros(NSIMUS,HORIZON);
    
	guarantee_begin_per          = zeros(NSIMUS,HORIZON);
	A1_realized_vals            = zeros(NSIMUS,NPERS);
	A2_realized_vals            = zeros(NSIMUS,NPERS);
	A1_paid                     = zeros(NSIMUS,NPERS);
	A2_paid                     = zeros(NSIMUS,NPERS);
	A1_defaults                 = true(1,NSIMUS);
	A2_defaults                 = true(1,NSIMUS);
    svc_defaults                = true(1,NSIMUS);
	return_on_equity            = zeros(1,NSIMUS);

	if (show_progress) 
        fprintf(1,'Simulation of %d paths: 0%% complete\n',NSIMUS);
    end

    for s = 1 : NSIMUS
		if(show_progress) 
			if (mod(s,prog_bar_click) == 0)
                fprintf(1,'Simulation %4.0f%% complete\n', (s/NSIMUS)*100);
            end
        end
		
		% reset variables for next path of simulation
		NCOMPOUNDS 			= sum(simu.initial_compounds(:));
		
		% state of compounds
		invested    		= false(1,NCOMPOUNDS); %true(1,NCOMPOUNDS);
        acquired            = false(1,NCOMPOUNDS);
		compounds   		= STARTING_STATES;
        %compounds           = compounds(randperm(length(compounds)));
		time_in_phase 		= zeros(1,NCOMPOUNDS);
        
        
		sell_value    		= zeros(1,NCOMPOUNDS);
		equity_in_compound 	= ones(1,NCOMPOUNDS).*assets.equity_stake;

		% state of securities
		cash			  	= zeros(1,HORIZON);
		guarantee           = zeros(1,HORIZON);
		
		IC_ratio        	= zeros(length(bonds.nominal)-1, NPERS);
		IC_shortfall   		= zeros(length(bonds.nominal)-1, NPERS);
		
		%record keeping 
		interest_paid	   	= zeros(length(bonds.nominal), NPERS);
		principal_paid     	= zeros(length(bonds.nominal), NPERS);
		
		bond_value          = bonds.nominal;
		cash(1)    = first_cash;
		cash_begin_per(s,1) = cash(1);
        
		guarantee(1) = first_guarantee;
        guarantee_begin_per(s,1) = guarantee(1);
        unpaid_svc  = 0;
        in_default  = false;
        invest_cost = trial_cost_fn(compounds,params);
        
		% temporary variables to make sure we spend only our target percentage
        money_spent =  zeros(1,NCOMPOUNDS);
		money_saved =  zeros(1, NCOMPOUNDS);
	    psindx      = phase2index_fn(assets.sell_in_phase)-1;
        cindxs     	=  phase2index_fn(compounds);
	    price_tmp 	=  invest_cost + assets.pricing_params(cindxs,assets.pricing_params_col.UpFront)';
        % If we have enough cash, buy it and update our initial cash
        %1/10*
        FutureCostEst  = (psindx>=1).*(assets.pricing_params(cindxs,assets.pricing_params_col.FutureCostEst)' - assets.pricing_params(psindx,assets.pricing_params_col.FutureCostEst));
            
        cash_left_to_spend = cash(1)-bonds.IC_pers*(bonds.nominal(1)*bonds.coupon(1)+bonds.nominal(2)*bonds.coupon(2)); 
          
        % number of compounds we actually buy per phase LAST PERIOD
    	num_compounds      =  zeros(1,NSTATES);
        for ctr = 1:NCOMPOUNDS
            if STARTING_PERIOD(ctr)==1
                if (cash_left_to_spend 	>= (price_tmp(ctr) + FutureCostEst(ctr)))
                    money_spent(ctr)     =  price_tmp(ctr);
                    money_saved(ctr)     =  FutureCostEst(ctr);
                    cash_left_to_spend   =  cash_left_to_spend - price_tmp(ctr) - FutureCostEst(ctr);
                    cash(1)              = cash(1) - price_tmp(ctr);
                    num_compounds(cindxs(ctr)) =  num_compounds(cindxs(ctr)) + 1;
                	invested(ctr)    = true;
                    acquired(ctr)    = true;
                    compounds_funded(s,1) = compounds_funded(s,1)+1;
                else 
                    compounds{ctr} = 'DSC';
                end
            end
        end
        compounds_bought(s,1,:) = num_compounds;				
            
        
		cash(2)             = cash(1);
		cash_begin_per(s,2) = cash(1);
        guarantee(2)             = guarantee(1);
        guarantee_begin_per(s,2) = guarantee(1);
        %----- simulate one trajectory ----
		%
		% I) For each time period
		%    1) Check for transitions
		%	 2) Apply waterfall rules
		%		A) if insufficient cash to pay this period's obligations transaction is in default
		%			i) liquidate portfolio 
		%			ii) use proceeds to pay servicing and then all P&I for A1, then P&I for A2, etc.
		%		B) if not in default
		%			i) check coverage test 
		%			ii) If coverage is not OK 
		%				a) sell enough assets to bring ratio back in line (coarse approximation)
		%			iii)Pay obligations 
		%				a) Pay servicing
		%				b) Pay P&I for A1, then for A2, etc.
		%				c) make investmensts in compounds as cash permits

		% II) If last time period
		%	1) liquidate portfolio
		%	2) Make final payements to bonds
		%	3) Residual cash goes to equity
		%

		A1_in_default = false;
		A2_in_default = false;
        svc_in_default = false;
        %Valuation for average market component.
		z = randn(1);

		%%TEMPORARY CODE FOR DURATION DISTRIBUTION
        if strcmp(params.assets.mode,'General')
            dur = [0, dur_mean, 0];
            v = dur/4;
            mu_dur = log(dur.^2./sqrt(v+dur.^2));
            sigma_dur = sqrt(log(1+v./dur.^2));
        
    		time_required = zeros(NCOMPOUNDS,NSTATES);
		
    		for i = 2:NSTATES-1
                if strcmp(assets.distribution,'Constant')
                    %% Constant Duration.
                    time_required(:,i) = round(2*dur(i)*ones(1,NCOMPOUNDS));
                elseif strcmp(assets.distribution,'LogNormal')
                    %% Log-Normal Duration.
                    time_required(:,i) = 2*lognrnd(mu_dur(i),sigma_dur(i),1,NCOMPOUNDS);
                elseif strcmp(assets.distribution,'Geometric')
                    %% Geometric Duration.
                    time_required(:,i) = 1+geornd(1-assets.trans_prob(i,i),1,NCOMPOUNDS);
                else
                    error('Unknown distribution.');
                end
            end
        end
		for (i = 2:NPERS) % during life of bonds

            current_cash  = cash(i);
			current_guarantee  = guarantee(i);
    
            % ----------- first determine current state of portfolio
            old_phase = compounds;
            
            elig_idxs = invested & ~(strcmp('DSC',compounds) | strcmp('SLD',compounds) | strcmp('APP',compounds));
            if strcmp(params.assets.mode,'General')==1
                compounds = check_for_transitions_ffn_distribution(compounds, invested, time_required, success, params);
            else
                compounds = check_for_transitions_ffn(compounds, invested, rho, params);
            end
            sell_idxs = elig_idxs & (strcmp(assets.sell_in_phase,compounds) | strcmp('APP',compounds)); 
            pr = sell_compound_ffn(sell_idxs,i);
           
    	    dsc_idxs = elig_idxs & strcmp('DSC',compounds);
            withdraw_compound_ffn(old_phase(dsc_idxs),i);
            time_in_phase(dsc_idxs) = 0;
            trans_idxs = elig_idxs & (~sell_idxs) & (~dsc_idxs) & (~strcmp(old_phase,compounds)); %transiton occured
           
            transition_compound_ffn(trans_idxs);
            stat_idxs = elig_idxs & ~sell_idxs & ~dsc_idxs & ~trans_idxs;
            time_in_phase(stat_idxs) = time_in_phase(stat_idxs)+1;
			
			% ----------- accounting balance
			A1_realized_vals(s,i) = bond_value(1);
			A2_realized_vals(s,i) = bond_value(2);
			
			%------------ next if not yet in default, check for default
			if(~in_default) 
              in_default = check_for_default_ffn(i,current_cash,current_guarantee);
			end %if ~in_default
			
			%------------ next if not in default pay P&I; if in default liquidate portfolio
			if (~in_default)
                %Pay service fees
                [current_cash, current_guarantee] = pay_servicing_ffn(i,current_cash,current_guarantee);
				%Pay P&I
                [current_cash, current_guarantee] = pay_interest_ffn(i,current_cash,current_guarantee);
                [current_cash, current_guarantee] = pay_principal_ffn(i,current_cash,current_guarantee);
				
				%------------ next check coverage
				for(b = 1:NBONDS) %calc IC coverage
					if(bond_value(b)>0)		
						cash_for_coverage=current_cash+sum(cash(i+(1:params.bonds.IC_pers)));
                        [IC_ratio(b,i) cost_due K] = ic_test_result_ffn(b,i,cash_for_coverage);  
						
						IC_shortfall(b,i) = (bonds.interest_coverage(b)*cost_due) - K;                     					
					else 
                        IC_shortfall(b) = 0;
                    end %if bond_value(b)>0
				end %for b = 1: NBONDS

                
				%------------ check equity svc coverage
				eq_svc_due = unpaid_svc;
				for (k = 1:ceil(max(assets.sale_time)))
					eq_svc_due = eq_svc_due + bond_value(length(bond_value))*bonds.servicing_rate;
				end
				eq_shortfall =  eq_svc_due - current_cash;
				
				%------------ next if IC too low, sell compounds to cover
				max_shortfall = max([IC_shortfall(:,i);eq_shortfall]); 
				if(max_shortfall > 0) 
					%Mutates cash,compounds.
					sell_compounds_to_cover_shortfall_ffn(max_shortfall,i);
					OK_to_fund = false;
				else 
					OK_to_fund = true;
				end
				
				%------------ next if IC is OK fund next phases of research for compounds 
				compounds_to_fund(s,i) = sum( ~strcmp(compounds,'DSC') & ~strcmp(compounds,'APP') & ~strcmp(compounds,'SLD') & acquired & ~invested);

				if(OK_to_fund)
					if(i<NPERS)
						if(i<NPERS-bonds.IC_pers) 
                            next_pers = (i+1):(i+bonds.IC_pers); 
							num       = bonds.IC_pers;
						elseif (i<(NPERS-bonds.IC_pers-1))
							next_pers = (i+1):(i+bonds.IC_pers-1);
							num       = bonds.IC_pers-1;
						else
							next_pers = i+1;
							num       = 1;
                        end
							
						cash_required_for_future = sum(sum(AMORT_SCHED(:,next_pers))) + num*( sum(bonds.coupon(1:NBONDS).*bond_value(1:NBONDS))+(bonds.servicing_rate*sum(bond_value))); 
                        
						cash_to_invest           = max(current_cash-cash_required_for_future,0);
                        
						if(cash_to_invest>0)
                            compound_nos           = 1:NCOMPOUNDS;
							compounds_to_fund_idx  = compound_nos(~strcmp(compounds,'DSC') & ~strcmp(compounds,'APP') & ~strcmp(compounds,'SLD') & acquired & ~invested); %index of compounds that need funding
							cash_spent             = 0;
							for k = compounds_to_fund_idx
								if(invest_cost(k) <= cash_to_invest) 
									invested(k)    = true;
									cash_to_invest = cash_to_invest-invest_cost(k);
									cash_spent     = cash_spent+invest_cost(k);
									compounds_funded(s,i) = compounds_funded(s,i)+1;
								end
                            end
                            
                            cash_for_investment(s,i) = cash_for_investment(s,i) + cash_spent;
                            current_cash = current_cash-cash_spent ;
                            cash_spent = 0;
                            cindxs     			=  phase2index_fn(compounds);
                            cindxs(cindxs<=0) 	=  1;
                            %%OLD PARAMETER FILES NEED THIS
                            %%%FutureCostEst  = (psindx>=1).*(assets.pricing_params(cindxs,assets.pricing_params_col.FutureCostEst)' - assets.pricing_params(psindx,assets.pricing_params_col.FutureCostEst));                        

                            %% NEW PARAMETER FILES NEED THIS
                            FutureCostEst = (cindxs>=1).*(assets.pricing_params(cindxs,assets.pricing_params_col.FutureCostEst)');                        
                            for ctr = 1:NCOMPOUNDS
                               if STARTING_PERIOD(ctr) < i & acquired(ctr)
                                   cash_to_invest = cash_to_invest - FutureCostEst(ctr);
                               end
                            end
                            
                            % number of compounds we actually buy per phase LAST PERIOD
                            num_compounds      =  zeros(1,NSTATES);
                            
                            for k = 1:NCOMPOUNDS
                                if STARTING_PERIOD(k)==i
                                    if (cash_to_invest 	>= (price_tmp(k) + FutureCostEst(k)))
                                        cash_to_invest   =  cash_to_invest - price_tmp(k) - FutureCostEst(k);
                                        cash_spent         = cash_spent + price_tmp(k);
                                        num_compounds(cindxs(k)) =  num_compounds(cindxs(k)) + 1;
                                        invested(k)    = true;
                                        acquired(k)    = true;
                                    else 
                                        compounds{k} = 'DSC';
                                    end
                                end
                            end
                            compounds_bought(s,i,:) = num_compounds;				
                            cash_for_buying(s,i) = cash_for_buying(s,i) + cash_spent;    

                            current_cash = current_cash-cash_spent ;
							
							
						end
                    end %if(i<NPERS)
                end %if(OK_to_fund)
			else % already in default

				liquidate_portfolio_ffn(i);
				[current_cash current_guarantee] = pay_servicing_ffn(i,current_cash, current_guarantee) ; % need to check to make sure missed payements are tracked

				%PAY P&I
				for(b = 1:NBONDS) % pay down bonds sequentially
                    
                    if(current_cash+current_guarantee>0)
                        I  = min(current_cash+current_guarantee,bonds.coupon(b)*bond_value(b));  
                        interest_paid(b,i)  = I;
                        current_cash        = current_cash-I;
                        
                        P                   = min(current_cash+current_guarantee,bond_value(b));
                        principal_paid(b,i) = P;
                        bond_value(b)       = bond_value(b)-P;
                        if(bond_value(b)<0.0001) 
                            bond_value(b) = 0;
                        end
                        current_cash = current_cash-P;
						if(current_cash<0)
                             [current_cash, current_guarantee] = draw_from_guarantee_ffn(i,current_cash, current_guarantee);
                        end
                    end
                end
				% Check to see if A1 defaults	
				if (~A1_in_default && A2_in_default)
					if (bond_value(1) > 0) 
						A1_in_default = true;
					end
				end
            end %else ~in_default

			cash(i) = current_cash;
			guarantee(i) = current_guarantee;
            
			if(i~=NPERS) 
				%cash left over at end of period acquires interest
				cash(i+1)             = cash(i+1)+(current_cash*(1+bonds.cash_accrual_int_rate));
				cash_begin_per(s,i+1) = cash(i+1);
                
				guarantee(i+1)        = guarantee(i);
                guarantee_begin_per(s,i+1) = guarantee(i+1);
            end
            
            
		end %for i = 2 : NPERS	%%END TIMESTEP ITERATION%%%%%

		final_compounds = final_compounds + [ sum(strcmp(compounds,'DSC')) sum(strcmp(compounds,'PRE')) sum(strcmp(compounds,'P1')) sum(strcmp(compounds,'P2')) sum(strcmp(compounds,'P3')) sum(strcmp(compounds,'NDA')) sum(strcmp(compounds,'APP'))]./NSIMUS;

		liquidate_portfolio_ffn(NPERS);
		%Experimental, untested repayment parameter.
        if(params.ce.ce_repay)
           [current_cash current_guarantee] = repay_guarantee_ffn(NPERS,current_cash, current_guarantee);
        end
        for(i = NPERS:(HORIZON-1)) 
			current_cash = cash(i);
			[current_cash current_guarantee] = pay_servicing_ffn(i,current_cash, current_guarantee);
			%Experimental, untested repayment parameter.
			if(params.ce.ce_repay)
                [current_cash current_guarantee] = repay_guarantee_ffn(i,current_cash, current_guarantee);
            end
            cash(i)=current_cash;
			guarantee(i)=current_guarantee;
			cash(i+1)    = cash(i)+cash(i+1); % last sales; need to add check for NPERS>=HORIZON
			guarantee(i+1) = guarantee(i);
            
			cash_begin_per(s,i+1) = cash(i+1);
        end
		cash_realized(s,:)    = cash;
        compounds_spent(s)  = sum(money_spent);
		A1_paid(s,:)   = principal_paid(1,:)+interest_paid(1,:);
		A2_paid(s,:)   = principal_paid(2,:)+interest_paid(2,:);
		A1_defaults(s) = A1_in_default;
		A2_defaults(s) = A2_in_default;
        svc_defaults(s)= svc_in_default;
       
        %keyboard
		return_on_equity(s) = max((cash(end) - bonds.nominal(bonds.nominal_col.EQ))/(bonds.nominal(bonds.nominal_col.EQ)),-1);
 
    end %for s = 1:NSIMUS 	###END SIM COUNT LOOP#####
    

    results = struct;
    
    results.runtime = toc(tStart);
    results.cash_begin_per      = cash_begin_per;
    results.cash                = cash_realized;
    results.cash_for_investment = cash_for_investment;
    results.cash_for_buying     = cash_for_buying;
    results.A1_bals           = A1_realized_vals;
    results.A2_bals           = A2_realized_vals;
    results.amort_sched       = AMORT_SCHED;
    results.A1_payment        = A1_paid;
    results.A2_payment        = A2_paid;
    results.IC_ratio          = IC_ratio;
    results.ROE               = return_on_equity;
    results.sale_phases           = sales;
    results.sale_phases_row       = sales_row;
    results.sale_phases_col       = sales_col;
    results.withdrawal_phases     = withdrawals;
    results.withdrawal_phases_row = withdrawals_row;
    results.withdrawal_phases_col = withdrawals_col;
    results.sale_times            = sale_times;
    results.sale_times_row        = sale_times_row;
    results.sale_times_col        = sale_times_col;
    results.withdrawal_times      = withdrawal_times;
    results.withdrawal_times_row  = withdrawal_times_row;
    results.withdrawal_times_col  = withdrawal_times_col;
    results.extra_cash_for_coverage    = extra_cash_from_coverage;
    results.servicing                  = history_service;
    results.compounds_to_fund              = compounds_to_fund;
    results.compounds_to_fund_row          = compounds_to_fund_row;
    results.compounds_to_fund_col          = compounds_to_fund_col;
    results.compounds_funded               = compounds_funded;
    results.compounds_funded_row           = compounds_funded_row;
    results.compounds_funded_col           = compounds_funded_col;
    results.compounds_bought     = compounds_bought;
    results.compounds_bought_row = compounds_bought_row;
    results.compounds_bought_col = compounds_bought_col;
    results.compounds_spent     = compounds_spent;
    results.params            = params;
    results.final_sale_values = sell_value;
    results.compounds_left    = final_compounds;
    results.A1_defaults       = A1_defaults;
    results.A2_defaults       = A2_defaults;
    results.svc_defaults      = svc_defaults;
	results.guarantee_begin_per = guarantee_begin_per;
    results.guarantee_draw      = guarantee_draw;
    results.drew_on_guarantee   = drew_on_guarantee;
end

%Usage


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




