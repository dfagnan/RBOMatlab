COPYRIGHT 2015

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
=========================================================================================================================

INSTRUCTIONS: 

This folder contains MATLAB code for the following publications:
Fernandez, Jose-Maria, Roger M. Stein, and Andrew W. Lo. "Commercializing biomedical research through securitization techniques." Nature biotechnology 30.10 (2012): 964-975.
Fagnan, David E., et al. "Can financial engineering cure cancer?." The American Economic Review 103.3 (2013): 406-411.
Fagnan, David E., et al. "Financing drug discovery for orphan diseases." Drug discovery today 19.5 (2014): 533-538.
Fagnan, David E., et al. "Financing translation: Analysis of the NCATS rare-diseases portfolio." Science translational medicine 7.276 (2015): 276ps3-276ps3.

Each is provided as an init params file in the "Testing" folder.

For convenience the script therein "run_tests" runs each of these to validate the results of all the papers, 
which is useful for testing the consistency of the framework after making modifications to the simulation assumptions.

By using these examples you can create your own simulations using a new init_params file.

Brief instructions are provided for each language below, but see the toolbox for an explanation of the key functions.

For running the MATLAB code:
1. Add the matlab directory to your path at the matlab prompt using the function addpath.
2. Run the following command at the matlab prompt to produce a simulation results object
>> results=simulate_CFs_fn();
3. You can view some statistics and plots of the resultant simulation results object using
>> summarize_results_fn(results);
4. Feel free to change the parameters of the simulation by editing the file
init_params_fn.m


Please contact David Fagnan dfagnan@mit.edu for technical support for any issues with using or modifying the software
(or even if you want to know more about trying Julia).
=========================================================================================================================

David Fagnan, Austin Gromatzky, Jose-Maria Fernandez, Roger M. Stein and Andrew W. Lo

In Cambridge, Massachusetts on the 11th of July 2015

