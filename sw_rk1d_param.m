% Sat 31 Aug 12:49:37 CEST 2024
% Karl Kastner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
function param = sw_rk1d_param(R,cva,psl,T,L,dx)
	% model parameters
	param         = struct();

	% spatial extent of computational domain (determines spectral resolution)
	param.L  = L;

	% number of grid cells
	param.nx = param.L(1)/dx;

	% simulated time interval
	param.T  = [0, T];

	% advection in x-direction
	param.pmu.vx = {0, 0, 10};
	% diffusion in x-direction
	param.pmu.ex = {0.1, 0.1, 50};
	% precipitation
	param.pmu.R = R;
	param.opt.R0 = param.pmu.R;
 
	% mean infiltration coefficient
	param.pmu.a = 0.2;
	    
	% spatial heterogeneity of parameter a
	param.pss.a  = param.pmu.a*cva;
        if (cva ~= 0)
		param.psdist.a = 'geometric-ornstein-uhlenbeck';
		param.psl.a    = psl;
       	end

	% boundary condition
	param.boundary_condition    = {'circular','circular'};
	    
	% data type during compution
	param.opt.compute_class = @single;
	% reload values of intermediate time steps
	% TODO this should not be a param, but argument to run
	param.opt.loadfinal = false;
	% data type of output file
	param.opt.output.class  = @half;
	% time step for storing output
	param.opt.output.dt     = T/100;
	param.opt.rng           = 0;
	% activate time step adaptation
	param.opt.time_integration.adapt_step_length=1;
	% initial time step
	param.opt.time_integration.dt_initial   = 1/400;
	param.opt.time_integration.dt_max       = inf;
	param.opt.time_integration.dt_min       = 1/400;
	param.opt.time_integration.abstol       = 1e-3;
	param.opt.time_integration.reltol       = 1e-2;
	param.opt.time_integration.dt_max_scale_up   = sqrt(2);
	param.opt.time_integration.dt_min_scale_down = 0;
	param.opt.time_integration.dt_safety_factor = 0.9;
	% integrate in time with trapezoidal scheme
	param.opt.time_integration.q = 0.5;
	% integrate in frequency space
        %param.opt.time_integration.scheme = 'step_advect_diffuse_spectral';
	%param.opt.time_integration.scheme = 'step_advect_diffuse_implicit_q_fft';
	param.opt.time_integration.scheme = 'step_advect_diffuse_implicit_q';
	
	% random initial condition
	param.initial_condition = 'obj.random_state()';

	% output directory
        param.opt.output.path_str = sprintf('mat/rk-vx-%g-ex-%g-Lc-%g-Lx-%g-T-%g/',param.pmu.vx{3},param.pmu.ex{3},psl,L,T);
end

