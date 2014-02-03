%Generate random draws from a lognormal distribution
function r = lognrnd(mu,sigma, varargin)

    if nargin < 2
        error(message('stats:lognrnd:TooFewInputs'));
    end

    [err, sizeOut] = statsizechk(2,mu,sigma,varargin{:});
    if err > 0
        error(message('stats:lognrnd:InputSizeMismatch'));
    end

    % Return NaN for elements corresponding to illegal parameter values.
    sigma(sigma < 0) = NaN;

    r = exp(randn(sizeOut) .* sigma + mu);


end


function [err, commonSize, numElements] = statsizechk(nparams,varargin)

    try
        tmp = 0;
        for argnum = 1:nparams
            tmp = tmp + varargin{argnum};
        end
        if nargin > nparams+1
            tmp = tmp + zeros(varargin{nparams+1:end});
        end
        err = 0;
        commonSize = size(tmp);
        numElements = numel(tmp);

    catch
        err = 1;
        commonSize = [];
        numElements = 0;
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
