%%Input 	: steady - vector (length n) of success probabilities (steady-state)
%%			: dur - vector (length n) of expected time to transition 
%%Optional input
%% 		: timestep_per_unit - number of time steps per time unit 
%%Output: transition matrix (n+2 x n+2 ) for all states, adding row/column 
%% 		  for withdrawn (row 1) and approved (row n+2) states.

function prob=createProbMat(steady,dur,varargin)
	p = inputParser;
	p.addRequired('steady',@(x) all(x>=0) & all(x<=1));
	p.addRequired('dur',@(x) all(x>0));
	
	%%Default to 2 timesteps.
	p.addOptional('timestep_per_unit',2,@isnumeric);
	p.parse(steady,dur,varargin{:});
	n=length(steady);
	prob = zeros(n+2,n+2);
	
	%%withdraw state is row one
	prob(1,1)=1;
	%%approved state is row n+2
	prob(n+2,n+2)=1;
	for i = 1:n
		prob(i+1,i+1)=(dur(i)*p.Results.timestep_per_unit-1)/(dur(i)*p.Results.timestep_per_unit);
		prob(i+1,i+2)=(1-prob(i+1,i+1))*steady(i);
		prob(i+1,1)=1-sum(prob(i+1,:));
	end	
end
