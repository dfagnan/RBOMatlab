function mat=createTransMat(appValue,costs,disc,dur,steady,target)
		n = length(disc);
		value = zeros(n+1,1);
		value(n+1) = appValue;		
		stdvalue=zeros(n+1,1);
		stdvalue(n+1) = appValue*(2240/1870);
		stdratio = [5.2/4.8, 15.2/12.8, 23.5/22.1, 86.3/60.6, 0];
		stdcost = zeros(n+1,1);
		for i = 1:n
			stdcost(i)=costs(i)*stdratio(i);
		end	

		disc_total = ones(n+1,1);
		futureCost = zeros(n+1,1);
		upfront    = zeros(n+1,1);
%		upfront    = (2,7.5,40,200,0,0)
		milestone  = [1, 3.75, 10, 37.6,0,0]';
	%	milestone  = zeros(n+1)
		
		
        cmx = [20, 20, 100, 500, 0, 0]'; % ORG cmx = [20, 20, 100, 500, 0, 0]';
		vmx = [100, 250, 500, 1000, 5000, 5000]'; % ORIG [100, 250, 500, 1000, 2500, 5000]';
        
		vmu = zeros(n+1,1);
		vsigma = zeros(n+1,1);
		csigma = zeros(n+1,1);
		cmu = zeros(n+1,1);
		
		vsigma(n+1)=sqrt(log(1+stdvalue(n+1)^2/value(n+1)^2));
		vmu(n+1) = log(value(n+1))-0.5*vsigma(n+1)^2;
	
		vsigma(n+1)=0.939;
		
		vmu(n+1)=convertLNtoCappedLN(vmu(n+1),vsigma(n+1),vmx(n+1),value(n+1));
	
	
		for i = n:-1:1
			disc_total(i)=1/(1+disc(i))^dur(i)*steady(i);
			if i<n
				futureCost(i)=(futureCost(i+1)+costs(i+1))*steady(i);
			end
		
			value(i)=max(appValue*prod(disc_total(i:n))-futureCost(i),0);
			upfront(i)=max(value(i)-costs(i),0);
%			milestone(i)=upfront(i)/2
		
			stdvalue(i)=value(i)*2240/1870;
			if value(i) < 1E-6 
				if vmx(i) >1E-6
					warning('Zero valuation and non-zero maximum not compatible, assuming identically zero.')
				end
				vmu(i)=0;
				vmx(i)=0;
				csigma(i)=1;
			else
				vsigma(i)=sqrt(log(1+stdvalue(i)^2/value(i)^2));
				vmu(i)= log(value(i))-0.5*vsigma(i)^2; 
				vmu(i)=convertLNtoCappedLN(vmu(i),vsigma(i),vmx(i),value(i));
			end
			vsigma(i)=0.939;
			
			if costs(i) < 1E-6
				if cmx(i) >1E-6
					warning('Zero cost and non-zero maximum not compatible, assuming identically zero.')
				end
				cmu(i)=0;
				cmx(i)=0;
				csigma(i)=1;
			else
				csigma(i)=sqrt(log(1+stdcost(i)^2/costs(i)^2));
				cmu(i)=log(costs(i))-0.5*csigma(i)^2;
				cmu(i)=convertLNtoCappedLN(cmu(i),csigma(i),cmx(i),costs(i));
			end
		end
		
		budget = zeros(n+1,1);
		for i = (target-2):-1:1
			budget(i) = steady(i)*(costs(i+1)+budget(i+1)+milestone(i+1));
		end
		%budget += upfront

		mat=[[1, 0, 0, 0, 0, 0, 0, 0, 0]; [vmu,vsigma,vmx,upfront,milestone,cmu,csigma,cmx,budget]];
end
