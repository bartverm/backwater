classdef Backwater < handle
% Computes backwater curves
%
%   Backwater properties (Constant):
%   g - acceleration of gravity, constant 9.81 m^2/s
%
%   Backwater properties:
%   Q - discharge (m^3/s)
%   b - width (m)
%   So - bed slope (-)
%   Chez - Chèzy coefficient (m^.5/s)
%   a0 - depth boundary condition for x=x0 (m)
%   x0 - location of boundary condition (m)
%   x_end - extent of computation (m)
%   target_fract - fraction for ending computation to eq. (-)
%   color_water - set the color of water for plotting (RGB)
%   color_bed - set the color of the bed for plotting (RGB)
%   bed_offset - additional thickness of bed for plotting (m)
%
%   Backwater properties (Computed, Read only):
%   a_equilibrium - equilibrium depth (m)
%   u_equilibrium - equilibrium velocity (m/s)
%   a_critical - critical depth (m)
%   a_target - depth at which computation stops (m)
%   x_target - target distance for computation (m)
%   Sc - critical slope for given friction (-)
%   slope_type - type of slope (backwater classification)
%   curve_type - type of curve (backwater classification)
%   bed_level - bed level elevation (m)
%   is_supercritical - whether flow is supercritical (true|false)
%   is_equilibrium - whether flow is in equilibrium (true|false)
%   x_curve - the x-coordinate of the solved backwater (m)
%   a_curve - depth profile of the backwater (m)
%   dadx_curve - x-gradient of depth (-)
%   u_curve - velocity profile along the backwater (m/s)
%   dudx_curve - x-gradient of velocity (-)
%   qs_curve - sediment transport rate along backwater (m^2/s)
%   dqsdx_curve - gradient of sediment transport along backwater (m/s)
%   dzbdt_curve - erosion/deposition rate of the bed along backwater (m/s)
%   a_morf_equilibrium - depth at morphological equilibrium (m)
%   So_morf_equilibrium - bed-slope at morphological equilibrium (-)
%
%   Backwater methods:
%   solve - solve the backwater curve numerically
%   bresse - analytical backwater solver
%   plot - plot the backwater curve
%   plot_velocity - plot the velocity profile
%   plot_vel_gradient - plot the velocity gradient
%   plot_qs_gradient - plot sediment transport gradient
%   plot_zb_gradient - plot erosion/deposition rate
%   plot_gradient - generic gradient plotting function
%   plot_curve_interfaces - plot dashed lines at transition of curves
%   plot_initial_ersed - plot backwater and gradient for init eros/sediment
%   plot_morf_equilibrium - plot the morphological equilibrium
%
%   Backwater static methods:
%   belanger_static - Equation of Belanger (for use with ode solver)
%   bresse_int - Solution of Bresse's integral
%   bresse_static - Bresse's solution to equation of Belanger
%   bresse_zero_static - Bresse's solution modified for zero bed slope
%
%   Copyright 2020, Bart Vermeulen

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

    properties(Constant)
        % Backwater/g - gravity
        %
        %   Acceleration of gravity. Constant property equal to 9.81 m^2/s
        %
        % see also: Backwater
        g(1,1) double = 9.81; % accelaration of gravity (L T^-2)
    end

    properties(SetObservable)
        %% set input variables
        % Backwater/Q - Discharge
        %
        %   Water discharge (m^3/s). Requirements: double, scalar, finite, 
        %   real and positive. Default=1000 m^3/s.
        %
        % see also: Backwater
        Q (1,1) double {mustBeFinite, mustBeReal, mustBePositive} = 1000
        
        % Backwater/b - width
        %
        %   Width of the channel (m). Requirements: double, scalar, finite,
        %   positive and real. Default=100 m.
        %
        % see also: Backwater
        b (1,1) double {mustBePositive, mustBeFinite, mustBeReal} = 100
        
        % Backwater/b - bed slope
        %
        %   Bed slope of the river (-). Defined positive for downsloping
        %   river, i.e. So=-dz/dx. Requirements: double, scalar, finite,
        %   and real. Default=1e-4.
        %
        % see also: Backwater
        So (1,1) double {mustBeFinite, mustBeReal} = 1e-4; % bed slope (-)
        
        % Backwater/Chez - Chezy coefficient
        %
        %   Chezy coefficient for rough (m^.5/s). Defined positive for 
        %   downsloping river, i.e. So=-dz/dx. Requirements: double, 
        %   scalar, finite, real and positive. Default = 50 m^.5/s
        %
        % see also: Backwater        
        Chez (1,1) double {mustBeFinite, mustBePositive, mustBeReal} = 50; 

        % Backwater/a0 - Depth at boundary condition
        %
        %   Depth to impose on boundary x=x0 (m). Requirements: double, 
        %   scalar, finite, real and positive. Default = 4 m.
        %
        % see also: Backwater        
        a0 (1,1) double {mustBeFinite, mustBePositive, mustBeReal} = 4;
        
        % Backwater/x0 - x coordinate of boundary condition
        %
        %   x coordinate of boundary condition (m). Requirements: double, 
        %   scalar, finite, real. Default = 0 m.
        %
        % see also: Backwater 
        x0 (1,1) double {mustBeFinite, mustBeReal} = 0;
        
        % Backwater/x_end - ending x coordinate for backwater
        %
        %   End coordinate of computation (m). Requirements: double, 
        %   scalar, finite, real. Default = -40000 m. An estimate for the
        %   ending coordinate can be obtained from the x_target property.
        %
        % see also: Backwater, Backwater/x_target, Backwater/solve
        x_end (1,1) double {mustBeFinite, mustBeReal} = -40000;
    end
    properties
        m_sed_transp (1,1) double = 2.3e-4 
        n_sed_transp (1,1) double = 5
        porosity (1,1) double = 0.6
        qs_morf_eq (1,1) double = 1.3e-4
        
        % Backwater/target_fract - fraction to estimate end of computation
        %
        %   Fraction of the difference between the initial depth and the
        %   equilibrium depth to estimate x_target when approaching
        %   equilibrium depth (-). Requirements: scalar, finite, positive,
        %   real and less than 1. Default: exp(3).
        %
        %   see also: Backwater, Backwater/x_target, Backwater/a_target
        target_fract(1,1) double {mustBeFinite, mustBePositive, mustBeLessThan(target_fract,1), mustBeReal} = exp(-3);
        
        % Backwater/color_water - color of the water used in plotsobj
        %
        %   Defines the color of water as RGB values. Requirements: double,
        %   3 element row vector, finite, >=0, <=1, real.
        %   Default:[204, 204, 253]/255
        %
        %   see also: Backwater, Backwater/color_bed, Backwater/bed_offset
        %             Backwater/plot
        color_water(1,3) double {mustBeFinite, mustBeNonnegative, mustBeLessThanOrEqual(color_water, 1), mustBeReal} = [204, 204, 253]/255

        % Backwater/color_bed - color of the bed used in plots
        %
        %   Defines the color of water as RGB values. Requirements: double,
        %   3 element row vector, finite, >=0, <=1, real.
        %   Default:[204, 204, 253]/255
        %
        %   See also: Backwater, Backwater/color_water Backwater/bed_offset
        %             Backwater/plot
        color_bed(1,3) double {mustBeFinite, mustBeNonnegative, mustBeLessThanOrEqual(color_bed, 1), mustBeReal} = [204, 178, 178]/255
        
        zb0(1,1) double {mustBeFinite, mustBeReal}=0;
        
        % Backwater/bed_offset - additional thickness of the bed in plots
        %
        %   Sets an additional thickness of the bed in plots (m).
        %   Requirements: double, scalar, finite, >=0, real. Default: 0 m.
        %
        %   See also: Backwater, Backwater/color_water Backwater/color_bed
        %             Backwater/plot
        bed_offset(1,1) double {mustBeFinite, mustBeNonnegative, mustBeReal} = 0
    end
    properties(Dependent, SetAccess=private)
        % Backwater/a_equilibrium - computes the equilibrium depth (m)
        %
        %   Compute the equilibrium depth (m). The value can be complex and
        %   infinite
        %
        %   See also: Backwater, Backwater/a_critical
        a_equilibrium

        % Backwater/e_equilibrium - computes the equilibrium velocity (m/s)
        %
        %   Compute the equilibrium velocity (m/s). 
        %
        %   See also: Backwater, Backwater/a_equilibrium
        u_equilibrium
        
        % Backwater/a_critical - computes the critical depth (m)
        %
        %   Compute the critical depth (m). The value is always positive
        %   and marks the transition between subcritical and supercritical
        %   flow.
        %
        %   See also: Backwater, Backwater/a_equilibrium
        a_critical
        
        % Backwater/a_target - computes a target depth for computations (m)
        %
        %   Compute a logical depth at which the computation should end.
        %   This equals the critical depth, when the backwater curve tends
        %   to the critical depth, and it equals a fraction of the initial
        %   offset between boundary condition depth and equilibrium depth
        %   when the curve thends to approach the equilibrium depth. This
        %   is used to analytically evaluate x_target with analytical
        %   solutions of Bresse. For subcritical adverse and horizontal
        %   backwater curves, that do not approach any value, but are ever
        %   increasing away from the boundary condition the target depth
        %   will equal the boundary depth plus three times the difference
        %   between the boundary depth and the critical depth
        %
        %   See also: Backwater, Backwater/x_target, Backwater/x_end
        a_target
        
        % Backwater/x_target - returns a logical value for x_end
        % 
        %   From the target depth a_target the distance at which this
        %   target depth is reached is computed usin analytical solutions
        %   of bresse. This value can be used as a first guess for x_end.
        %
        %   See also: Backwater, Backwater/x_target, Backwater/x_end,
        %             Backwater/bresse
        x_target
        
        % Backwater/Sc - compute the critical slope (-)
        %
        %   The critical slope defines the boundary between a mild slope
        %   and the critical slope. It is a function of gravity and Chezy
        %   only. When the slope is smaller than this value, the
        %   equilibrium depth will be larger than the critical depth, i.e.
        %   when the river is in equilibrium it will be in subcritical
        %   flow. When the slope is larger than this value the equilibrium
        %   depth will be smaller than the critical depth, which means that
        %   the river in equilibrium is in supercritical flow.
        %
        %   See also; Backwater, Backwater/So, Backwater/a_critical,
        %             Backwater/a_equilibrium.
        Sc
        
        % Backwater/slope_type - slope type for backwater classification
        %
        %   This returns the type of slope based on the value of So and Sc
        %   (which only depends on the Chezy coefficient):
        %   'A' when So < 0
        %   'H' when So = 0
        %   'M' when 0 > So > Sc
        %   'C' when So = Sc
        %   'S' when So > Sc
        %
        %   See also: Backwater, Backwater/So, Backwater/Sc,
        %             Backwater/curve_type
        slope_type
        
        % Backwater/curve_type - curve type for backwater classification
        %
        %   This returns the type of curve base on the boundary condition
        %   and the equilibrium and critical depth:
        %
        %   1 when a0 > ae and a0 > ac
        %   2 when a0 is between ae and ac
        %   3 when a0 < ae and a0 < ac
        %
        %   See also: Backwater, Backwater/a_critical, 
        %       Backwater/a_equilibrium, Backwater/slope_type
        curve_type
        
        % Backwater/bed_level - starting and ending bed-level
        %
        %   starting and ending bed-level (m)
        %
        %   See also: Backwater, So
        bed_level
        
        % Backwater/is_supercritical - whether flow is supercritical
        %
        %   returns boolean indicating whether flow is supercritical
        %
        %   See also: Backwater, is_equilibrium
        is_supercritical
        
        % Backwater/is_equilibrium - whether flow is in equilibrium
        %
        %   returns boolean indicating whether flow is in equilibrium, or
        %   uniform
        %
        %   See also: Backwater, is_supercritical
        is_equilibrium
        
        % Backwater/x_curve - x-coordinate of backwater solution
        %
        %   returns the x-coordinate of the backwater solution (m). All
        %   properties ending with '_curve' are given at the locations in
        %   the x_curve property. This property is the output of the solve
        %   method, and is cached. When changing a property it is
        %   recomputed. Using this instead of directly calling the solve
        %   routine saves unnecessary recomputations of the backwater curve
        %
        %   See also: Backwater, solve, a_curve, u_curve
        x_curve
        
        % Backwater/a_curve - depth of backwater solution
        %
        %   returns the depth along the backwater (m) at the locations
        %   given in the x_curve property. This property is the output of 
        %   the solve method, and is cached. When changing a property it is
        %   recomputed. Using this instead of directly calling the solve
        %   routine saves unnecessary recomputations of the backwater curve
        %
        %   See also: Backwater, solve, x_curve, u_curve
        a_curve
        
        % Backwater/dadx_curve - x-gradient of depth
        %
        %   returns the gradient in x-direction of the depth (-) at the
        %   locations given in the x_curve property.
        %
        %   See also: Backwater, x_curve, u_curve
        dadx_curve
        
        % Backwater/u_curve - velocity along the backwater
        %
        %   return the profile of velocity along the backwater (m/s) at the
        %   locations given in the x_curve property.
        %
        %   See also: Backwater, x_curve, a_curve
        u_curve
        
        % Backwater/dudx_curve - x-gradient of velocity
        %
        %   returns the gradient in x-direction of the velocity (1/s) at
        %   the locations given in the x_curve property.
        %
        %   See also: Backwater, x_curve, u_curve
        dudx_curve
        
        % Backwater/qs_curve - sediment transport rates along the backwater
        %
        %   returns the sediment transport rates along the backwater curve
        %   (m^2/s) at the locations given in x-curve.
        %
        %   See also: Backwater, x_curve, m_sed_transport, n_sed_transport
        qs_curve
        
        % Backwater/dqsdx_curve - x-gradient of sediment transport rate
        %
        %   returns the gradient in x-direction of the sediment transport
        %   rate (m/s) at the locations given in the x_curve property.
        %
        %   See also: Backwater, x_curve, qs_curve
        dqsdx_curve
        
        % Backwater/dzbdt_curve - erosion/sedimentation rates
        %
        %   returns the erosion and sedimentation rates (m/s) at the
        %   locations given in the x_curve property. Computed using Exner's
        %   equation.
        %
        %   See also: Backwater, x_curve, porosity
        dzbdt_curve
        
        % Backwater/a_morf_equilibrium - depth at  morphological equilibr.
        %
        %   returns the depth (m) at the morphological equilibrium.
        %
        %   See also: Backwater, qs_morf_eq
        a_morf_equilibrium
        
        % Backwater/So_morf_equilibrium - slope at morphological equilibr.
        %
        %   returns the bed-slope (-) at morphological equilibrium.
        %
        %   See also: Backwater, qs_morf_eq
        So_morf_equilibrium

    end
    properties(Access=protected)
        % whether the cached solution is still up to date
        solution_is_up_to_date=false
        
        % cached x solution
        x_cache=[]
        
        % cached a solution
        a_cache=[]
    end
    methods
        %%% Set and get methods %%%
        function set.Q(obj,val)
            obj.Q=val;
            obj.invalidate_solution();
        end
        function set.b(obj,val)
            obj.b=val;
            obj.invalidate_solution();
        end
        function set.So(obj,val)
            obj.So=val;
            obj.invalidate_solution();
        end
        function set.Chez(obj,val)
            obj.Chez=val;
            obj.invalidate_solution();
        end
        function set.x0(obj,val)
            obj.x0=val;
            obj.invalidate_solution();
        end
        function set.a0(obj,val)
            obj.a0=val;
            obj.invalidate_solution();
        end
       
        function val=get.a_equilibrium(obj)
            val=(obj.Q^2/obj.Chez^2/obj.So/obj.b^2)^(1/3);
        end
        function val=get.a_morf_equilibrium(obj)
            val=(obj.m_sed_transp/obj.qs_morf_eq).^(1/obj.n_sed_transp)*obj.Q./obj.b;
        end
        function val=get.So_morf_equilibrium(obj)
            val=obj.b/obj.Chez^2/obj.Q*(obj.qs_morf_eq/obj.m_sed_transp)^(3./obj.n_sed_transp);
        end
        function val=get.u_equilibrium(obj)
            val=obj.Q/obj.b/obj.a_equilibrium;
        end
        function val=get.a_critical(obj)
            val=(obj.Q^2/obj.g/obj.b^2)^(1/3); 
        end
        function at=get.a_target(obj)
            if obj.a0<obj.a_critical && obj.So<obj.Sc || obj.a0 > obj.a_critical && obj.So>obj.Sc
                at=obj.a_critical; % curve approaches critical depth               
            else
                if obj.So<eps
                    at=obj.a0+(obj.a0-obj.a_critical)*3; % What length scale for subcritical H or A slopes?
                else
                    at=obj.a_equilibrium+obj.target_fract*(obj.a0-obj.a_equilibrium);
                end
            end
        end
        function val=get.is_supercritical(obj)
            if obj.a0-obj.a_critical<eps
                val=true;
            else
                val=false;
            end
        end
        function val=get.is_equilibrium(obj)
            if abs(obj.a0-obj.a_equilibrium)<eps
                val=true;
            else
                val=false;
            end
        end
        function set.x_end(obj,val)
            if obj.is_supercritical
                if val < obj.x0
                    warning('End of curve must be downstream of bc in supercritical flow, negating end')
                    obj.x_end=-val;
                else
                    obj.x_end=val;
                end
            else
                if val > obj.x0
                    warning('End of curve must be upstream of bc in subcritical flow, negating end')
                    obj.x_end=-val;
                else
                    obj.x_end=val;
                end
            end
            obj.invalidate_solution();
        end
        function val=get.x_target(obj)
            if obj.is_equilibrium
                val=obj.a_critical/obj.Sc;
            else
                val=obj.x0+obj.bresse(obj.a_target);
            end
        end
        function val=get.Sc(obj)
            val=obj.g/obj.Chez^2;
        end
        function val=get.slope_type(obj)
            if obj.is_equilibrium
                val='';
                return
            end
            if obj.So<=-eps
                val='A';
            elseif abs(obj.So)< eps
                val='H';
            elseif obj.So<obj.Sc-eps
                val='M';
            elseif abs(obj.So-obj.Sc)<eps
                val='C';
            else
                val='S';
            end
        end
        function val=get.curve_type(obj)
            if obj.is_equilibrium
                val='';
                return
            end
            if obj.So>eps
                if obj.a0>obj.a_critical && obj.a0 > obj.a_equilibrium
                    val='1';
                elseif obj.a0<obj.a_critical && obj.a0<obj.a_equilibrium
                    val='3';
                else 
                    val='2';
                end
            else
                if obj.a0 >= obj.a_critical
                    val='2';
                elseif obj.a0 < obj.a_critical
                    val='3';
                end
            end                   
        end
        
        function val=get.bed_level(obj)
            val=Backwater.bed_level_static([obj.x0, obj.x_end],obj.x0,obj.zb0,obj.So);
        end
        
        function val=get.x_curve(obj)
            if ~obj.solution_is_up_to_date
                obj.set_solution_cache();
            end
            val=obj.x_cache;
        end
        function val=get.a_curve(obj)
            if ~obj.solution_is_up_to_date
                obj.set_solution_cache();
            end
            val=obj.a_cache;
        end   
        function val=get.u_curve(obj)
            val=obj.Q/obj.b./obj.a_curve;
        end
        function val=get.qs_curve(obj)
            val=obj.m_sed_transp*obj.u_curve.^obj.n_sed_transp;
        end
        function val=get.dadx_curve(obj)
            val=Backwater.belanger_static(obj.a_curve,obj.So,obj.Q./obj.b,obj.Chez,obj.g);
        end
        function val=get.dudx_curve(obj)
            val=-obj.Q./obj.b./obj.a_curve.^2.*obj.dadx_curve;
        end
        function val=get.dqsdx_curve(obj)
            val=-obj.n_sed_transp*obj.qs_curve./obj.a_curve.*obj.dadx_curve;
        end
        function val=get.dzbdt_curve(obj)
            val=-1/(1-obj.porosity).*obj.dqsdx_curve;
        end
        %%% Ordinary methods %%%
        function [x,a,c]=solve(obj)
        % Backwater/solve - solves the backwater curve
        %
        %   [x,a]=obj.solve() solve the backwater curve for the Backwater
        %   object obj. This solves the belanger equation defined in
        %   belanger_static. It will return a backwater curve between x0
        %   and x_end. The function returns the x coordinates where the
        %   depths a are determined. Makes use of the matlab ode45 solver.
        %
        %   See also: Backwater, Backwater/x0, Backwater/x_end, ode45,
        %             Backwater/belanger_static
            if isscalar(obj)
                if obj.x0==obj.x_end
                    x=[obj.x0; obj.x_end];
                    a=obj.a0*[1; 1];
                else
                    if obj.a_target==obj.a_critical && ...
                            ((obj.is_supercritical && obj.x_end>obj.x_target) ||...
                            (~obj.is_supercritical && obj.x_end<obj.x_target))
                        warning('Trying to solve backwater past critical depth, changing end of reach...')
                        obj.x_end=obj.x_target;
                    end
                    [x,a]=ode45(@(x,a) Backwater.belanger_static(a, obj.So, obj.Q./obj.b, obj.Chez, obj.g),[obj.x0, obj.x_end],obj.a0);
                end
                x=x';
                a=a';
                c=ones(size(x));
            else
                % Solve the pieces
                x=[];
                a=[];
                c=[];
                for co=1:numel(obj)
                    [tx,ta,tc]=solve(obj(co));
                    c=[c tc*co];
                    x=[x tx];
                    a=[a ta];
                end
            end
        end
        
        function dx=bresse(obj,a)
        % Backwater/bresse - solves the backwater curve analytically
        %
        %   dx=obj.bresse(a) computes at what distance from the boundary
        %   condition the depth a is reached. Makes use of bresse_static
        %   for the computations except when So = 0, in which case it uses
        %   bresse_zero_static 
        %
        %   See also: Backwater, Backwater/bresse_static, 
        %             Backwater/bresse_zero_static, Backwater/x_target,
        %             Backwater/a_target
            if abs(obj.So)<eps
                eta_0=obj.a0/obj.a_critical;
                eta_1=a/obj.a_critical;
                dx=Backwater.bresse_zero_static(obj.a_critical,obj.Chez,obj.g,eta_1,eta_0);
            else
                eta_0=obj.a0/obj.a_equilibrium;
                eta_1=a/obj.a_equilibrium;
                dx=real(Backwater.bresse_static(obj.So,obj.a_equilibrium,obj.Chez, obj.g,eta_1,eta_0));
            end
        end

        function plot(obj)
        % Backwater/plot - plots the backwater curve
        %
        %   plot(obj) plots the backwater curve, with all kinds of
        %   annotations indicating the equilibrium and critical depth and
        %   the location of the imposed boundary condition
        %
        %   See also: Backwater
        hold_status=get(gca,'NextPlot');
        if isscalar(obj)
            z_b=Backwater.bed_level_static(obj.x_curve,obj.x0,obj.zb0,obj.So);
            z_w=z_b+obj.a_curve;
            z_b=[obj.bed_level];
            z_c=z_b+obj.a_critical;
            z_e=z_b+obj.a_equilibrium;
            hold on
            patch([obj.x0 obj.x_end([1 1]) obj.x0([1 1])], [min(z_b)*[1 1]-obj.bed_offset  z_b(2)  z_b(1)  z_b(1)-obj.bed_offset],obj.color_bed, 'linestyle','none');
            patch([obj.x0 obj.x_end flip(obj.x_curve) obj.x0],[z_b flip(z_w) z_b(1)], obj.color_water, 'linestyle','none')
            if obj.is_supercritical
                h_align='left';
            else
                h_align='right';
            end
            if all(isreal(z_e)) && all(isfinite(z_e))
                plot(obj.x_curve([1 end]),z_e,'b-.');
                text(obj.x_curve(1),z_e(1),'a_e','color','b','verticalalignment','middle','HorizontalAlignment',h_align);
            end
            plot(obj.x_curve([1 end]),z_c,'r--');
            text(obj.x_curve(1),z_c(1),'a_c','color','r','verticalalignment','middle','HorizontalAlignment',h_align);
            if ~obj.is_equilibrium
                plot(obj.x0,obj.zb0+obj.a0,'b.','markersize',15)
                text(obj.x0,obj.zb0+obj.a0,'a_0','verticalalignment','middle','HorizontalAlignment',h_align,'color','b');
                idx=max(1,floor(numel(obj.a_curve)/2));
                text(obj.x_curve(idx), z_w(idx),[obj.slope_type,obj.curve_type],'VerticalAlignment','bottom','HorizontalAlignment',h_align)
            end
        else
            minzb=min([obj.bed_level]-reshape(repmat([obj.bed_offset],2,1),1,[]));
            for co=1:numel(obj)
                plot(obj(co))
                z_b=min(obj(co).bed_level)-obj(co).bed_offset;
                patch([obj(co).x0 obj(co).x_end([1 1]) obj(co).x0([1 1])], [minzb*[1 1] z_b([1 1]) minzb],obj(co).color_bed, 'linestyle','none');
            end
            obj.plot_curve_interfaces();
        end
        set(gca,'tickdir','out')
        xlabel('x (m)')
        ylabel('z (m)')
        axis tight
        set(gca,'NextPlot',hold_status);
        end
        
        function plot_velocity(obj)
        % Backwater/plot_velocity - plots the computed velocity profile
        %
        %   This function makes a plot of the computed velocity (m/s)
        %   against the x-coordinate (m). 
        %
        %   See also: Backwater, plot
            hold_status=get(gca,'NextPlot');
            if isscalar(obj)
                plot(obj.x_curve,obj.u_curve,'b-')
                hold on
                plot(obj.x_curve([1 end]), obj.u_equilibrium([1 1]),'k--')
                if obj.is_supercritical
                    h_align='left';
                else
                    h_align='right';
                end
                text(obj.x_curve(1),obj.u_equilibrium,'u_e','horizontalalignment',h_align,'verticalalignment','top')
                set(gca,'tickdir','out')
                xlabel('x (m)')
                ylabel('u (m/s)')
                axis tight
            else
                for co=1:numel(obj)
                    obj(co).plot_velocity()
                end
                obj.plot_curve_interfaces()
                set(gca,'NextPlot',hold_status)
            end
        end
        function plot_vel_gradient(obj)
            obj.plot_gradient('u_curve', 'dudx_curve', 'dudx (1/s)')
        end
        function plot_qs_gradient(obj)
            obj.plot_gradient('qs_curve', 'dqsdx_curve', 'dqsdx (m/s)')
        end
        function plot_zb_gradient(obj)
            obj.plot_gradient('qs_curve', 'dzbdt_curve', 'dzbdt (m/s)', true)
        end
        
        function plot_curve_interfaces(obj)
            xbnd=deal(nan(2,numel(obj)));
            for co=1:numel(obj)
                xbnd(:,co)=[obj(co).x_end;  obj(co).x0];
            end
            xbnd=unique(xbnd(:));
            xbnd([1 end])=[];
            hold on
            for cb=1:numel(xbnd)
                plot(xbnd(cb)*[1 1],ylim(),'k-.');
            end
        end
        function plot_initial_ersed(obj)
            a(1)=subplot(4,1,1);
            plot(obj)
            a(end+1)=subplot(4,1,2);
            plot_velocity(obj)
            a(end+1)=subplot(4,1,3);
            plot_qs_gradient(obj)
            a(end+1)=subplot(4,1,4);
            plot_zb_gradient(obj)
            linkaxes(a,'x');
        end
        function plot_gradient(obj,curve_name, gradient_curve_name, curve_label, invert_flag)
            if nargin<5
                invert_flag=false;
            end
            hold_status=get(gca,'NextPlot');
            if isscalar(obj)
                plot(obj.x_curve,obj.(gradient_curve_name),'b-','linewidth',1.5)
                hold on
                plot(obj.x_curve([1 end]), [0 0],'k--')
                set(gca,'tickdir','out')
                xlabel('x (m)')
                ylabel(curve_label)
                axis tight
            else
                % draw interfaces between curves
                for co=1:numel(obj)
                    obj(co).plot_gradient(curve_name,gradient_curve_name, curve_label)
                end
                obj.plot_curve_interfaces();
                
                % find infinite gradients
                x=[obj.x_curve];
                [x,idx]=sort(x);
                u=[obj.(curve_name)];
                u=u(idx);
                dudx=[obj.(gradient_curve_name)];
                dudx=dudx(idx);
                fbnd=find(diff(x)==0);
                arrow_dir=(-1).^(double(invert_flag))*sign(u(fbnd+1)-u(fbnd));
                f_cont=arrow_dir==0;
                fbnd(f_cont)=[];
                arrow_dir(f_cont)=[];
                arrow_idx=(arrow_dir+1)/2+1;
                quiver_x=x(fbnd);
                yl=ylim;
                for cbnd=1:numel(arrow_idx)
                    y_start=dudx(fbnd(cbnd));
                    y_cent=yl(arrow_idx(cbnd));
                    arrow_ysiz=0.05*(y_cent-y_start);
                    y_lr=y_cent-arrow_ysiz;
                    x_cent=quiver_x(cbnd);
                    plot(x_cent*[1 1], [y_start y_cent],'b','linewidth',1.5)
                    sgn=num2str(arrow_dir(cbnd),'%+d');
                    sgn(end)=[];
                    text(x_cent, y_lr,[sgn,'\infty'])
                end
                set(gca,'NextPlot',hold_status)
           end            
        end
        function plot_morf_equilibrium(obj)
            hold_status=get(gca,'nextplot');
            x0s=[obj.x0];
            if all([obj.is_supercritical])
                order='ascend';
            elseif all(~[obj.is_supercritical])
                order='descend';
            else
                error('Cannot plot morphological equilibrium when reaches with sub- and supercritical flow are mixed')
            end
            [x0s, idx]=sort(x0s,order);
            Sos=[obj.So_morf_equilibrium];
            Sos=Sos(idx);
            aes=[obj.a_morf_equilibrium];
            aes=aes(idx);
            xends=[obj.x_end];
            xends=xends(idx);
            zwend=cumsum(Backwater.bed_level_static(xends,x0s,0,Sos))+obj(idx(1)).zb0+aes(1);
            zw=[obj(idx(1)).zb0+aes(1) zwend(1:end-1); zwend];
            zw=zw(:);
            xs=[x0s;xends];
            xs=xs(:);
            zb=zw-reshape(repmat(aes,[2 1]),[],1);
            patch([xs; xs([end 1 1])], [zw; 0; 0; zw(1)],obj(1).color_water,'linestyle','none')
            hold on
            patch([xs; xs([end 1 1])], [zb; 0; 0; zb(1)],obj(1).color_bed,'linestyle','none')
            obj.plot_curve_interfaces()
            zb_orig=vertcat(obj(idx).bed_level)';
            zb_orig=zb_orig(:);
            h=plot(xs, zb_orig,'k--');
            xlabel('x (m)')
            ylabel('z (m)')
            legend(h,'original bed level')
            set(gca,'tickdir','out')
            axis tight
            set(gca,'nextplot',hold_status)
        end
    end
    methods(Access=protected)
        function set_solution_cache(obj)
            if isscalar(obj)
                [obj.x_cache,obj.a_cache]=obj.solve();
                obj.solution_is_up_to_date=true;
            else
                for co=1:numel(obj)
                    obj(co).set_solution_cache;
                end
            end
        end
        function invalidate_solution(obj)
            obj.solution_is_up_to_date=false;
        end
    end
    methods (Static)
        function dadx=belanger_static(a, So, q, C, g)
        % Backwater/belanger_static - definition of equation of belanger
        %
        %   dadx=belanger_static(a, So, q, C, g) computes the slope of
        %   the depth, given the current depth a, the bed slope So, the
        %   specific discharge q, the chezy coefficient and the
        %   acceleration of gravity g. This function can be used with the
        %   matlab builtin ode solvers.
        %
        %   See also: Backwater, Backwater/solve
            dadx=(So-q.^2/C.^2./a.^3)./(1-q.^2./g./a.^3);
        end
        function phi=bresse_int(eta)
        % Backwater/bresse_int - solution to Bresse's integral
        %
        %   phi=bresse_int(eta) this function returns the analytical
        %   solution to the integral of Bresse, needed for Bresse's
        %   solution of the equation of Belanger.
        %
        %   See also: Backwater, Backwater/bresse, Backwater/bresse_static
            phi=1/6*log((eta.^2+eta+1)./(eta-1).^2)+1/sqrt(3)*(atan((2*eta+1)/sqrt(3))-atan(1/sqrt(3)));
        end
        function dx=bresse_static(So,ae, C, g, eta1, eta0)
        % Backwater/bresse_static - Bresse's solution to Belanger's eq.
        %
        %   dx=bresse_static(So,ae, C, g, eta1, eta0) Compute the distance
        %   from the x location where the normalized depth eta0 is imposed,
        %   to the location where the normalized depth eta1 is reached. The
        %   depth is normalized with the equilibrium depth. So is the bed
        %   slope, ae the equilibrium depth, C the coefficient of Chezy, g
        %   the acceleration of gravity. Does not work for horizontal bed
        %   slope. In that case use bresse_zero_static.
        %
        %   See also: Backwater, Backwater/bresse,
        %             Backwater/bresse_zero_static
            Ladapt=ae/So;
            gamma=(1-C^2*So/g);
            dx=Ladapt*(eta1-eta0-gamma*(Backwater.bresse_int(eta1)-Backwater.bresse_int(eta0)));            
        end
        function dx=bresse_zero_static(ac, C, g, eta1, eta0)
        % Backwater/bresse_zero_static - Bresse's solution for flat bed
        %
        %   dx=bresse_static(ac, C, g, eta1, eta0) Compute the distance
        %   from the x location where the normalized depth eta0 is imposed,
        %   to the location where the normalized depth eta1 is reached. The
        %   depth is normalized with the critical depth. ac is the 
        %   critical depth, C the coefficient of Chezy, g  the 
        %   acceleration of gravity. Only works for horizontal slopes.
        %
        %   See also: Backwater, Backwater/bresse, Backwater/bresse_static
            Sc=g/C^2;
            dx=ac/Sc*(eta1-eta0+.25*(eta0.^4-eta1.^4));
        end
        function z=bed_level_static(x, x0, zb0, So)
            z=zb0-(x-x0).*So;
        end
    end
end
