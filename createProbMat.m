function prob=createProbMat(steady,dur)
	n=length(steady);
	prob = zeros(n+2,n+2);
	prob(1,1)=1;
	prob(n+2,n+2)=1;
	for i = 1:n
		prob(i+1,i+1)=(dur(i)*(12/6)-1)/(dur(i)*(12/6));
		prob(i+1,i+2)=(1-prob(i+1,i+1))*steady(i);
		prob(i+1,1)=1-sum(prob(i+1,:));
	end	
end
