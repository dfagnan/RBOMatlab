%Input:
%@mu is of type double and specifies the mean of a normal distribution
%@sigma is of type double and specified the standard deviation of the normal distribution
%@rho is of type double and acts as a control parameter which if 0 returns fraws from the corresponding lognormal distribution else it returns draws from a biased log normal 
%@z1 is of type double and specifies the amount of bias if rho is non-zero
%@mx is of type double and specifies the maximum value beyong which draws from the lognorma distribution are rejected
%Output:
%@value is of type double and contains the resulting draw of the log-normal distribution
function value = compound_sale_price_fn(mu, sigma, mx, rho, z1)
    if(length(mu)~=length(sigma) || length(mu) ~= length(mx))
        error('Mu and Sigma different lengths.');
    end
    if (rho==0)
        value = min(lognrnd(mu, sigma), mx);
    else
        ej    = randn(length(mu),1);
        Z     = (sqrt(rho)*z1)+(sqrt(1-rho)*ej);
        X     = Z.*sigma + mu;
        value = min(exp(X),mx)';
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
