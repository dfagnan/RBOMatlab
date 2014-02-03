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
        else
			options = optimset('LargeScale','off','Display','off');
            vmu=fminunc(@(x) errorCapped(x,vsigma,vmx,value),vmu,options);
		end
		
end

