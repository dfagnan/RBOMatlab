%% Function for converting log-normal distribution to a capped log-normal distributions
%% Inputs : vmu - mu parameter for (uncapped) log-normal distribution
%%			vsigma - sigma parameter for (uncapped) log-normal distribution
%% 			vmx - cap on log-normal distribution
%%			value - desired expected value (mean) used for fitting
%% Outputs: output mu parameter for capped log-normal distributions.
function vmu=convertLNtoCappedLN(vmu,vsigma,vmx,value)
	
	function err=errorCapped(vmu,vsigma,vmx,value)
		P = normcdf(log(vmx),vmu(1),vsigma);
		err=1/2*(value-vmx*(1-P)-temp(vmu(1),vsigma,vmx)*P)^2;
	end

	function g=temp(vmu,vsigma,vmx)
		vmx_scaled = (log(vmx)-vmu)/vsigma;
		g = exp(vmu+vsigma^2/2)*normcdf(vmx_scaled-vsigma)/normcdf(vmx_scaled);
	end
	if vmx == 0
		vmu = 0;
	elseif isnan(vmu) 
		error('Log-normal mu parameter is NaN, possibly due to zero valuation or cost.')
    else
		options = optimset('LargeScale','off','Display','off');
        vmu=fminunc(@(x) errorCapped(x,vsigma,vmx,value),vmu,options);
	end
		
end

