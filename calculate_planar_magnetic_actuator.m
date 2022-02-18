clear;
%close;

%% Options
use_render = true;
use_calculate_resistance = true;
use_resistance_visualization = true;

%% Constants
xsize = 45; %board size x, mm
ysize = 100; %board size y, mm
N = 10; %winding count
dp = pi/1000; %angle step
clr = 1; %clearance between lines, mm
start_angle = 0*pi/4; %angle of starting point
clearance_angle_for_check = pi/10; %angle in which we check clearance

%% var
angles = 0:dp:2*pi-dp; %angle grid
len = length(angles); %number of points in grid (single winding_
image = zeros(N+1, len); %main results are stored in here (radiuses of windings for each angle)
angle1 = atan(ysize/xsize);
angle2 = pi-angle1;
angle3 = pi+angle1;
angle4 = 2*pi-angle1;
offset_angle_index = ceil(start_angle/dp)+1;
clearance_points_for_check = ceil(clearance_angle_for_check / dp);

sinp = sin(dp);
cosp = cos(dp);

sinp2 = sin(dp/2);
tgp2 = tan(dp/2);
cosp2 = cos(dp/2);

%% bounding box
for i=1:len
    %uncomment these 2 lines to have round bounding
    %image(1,i) = xsize/2;
    %continue;
    if angles(i) < angle1
        image(1,i) = xsize/2 / cos(angles(i));
        continue;
    end
    if angles(i) < angle2
        image(1,i) = ysize/2 / sin(angles(i));
        continue;
    end
    if angles(i) < angle3
        image(1,i) = -xsize/2 / cos(angles(i));
        continue;
    end
    if angles(i) < angle4
        image(1,i) = -ysize/2 / sin(angles(i));
        continue;
    end
    image(1,i) = xsize/2 / cos(angles(i));
end
image(1,:) = [image(1,offset_angle_index:end) image(1,1:offset_angle_index-1)];
angles = angles([offset_angle_index:end, 1:offset_angle_index-1]);

%% coil
r_prev = sqrt(xsize^2+ysize^2)*0.47;
r_next = sqrt(xsize^2+ysize^2)*0.47;
for i=2:N+1
    for j=1:len
		if i==2 && j==len
			i=i;
		end
		%rs - radiuses of points to check clearance
		%ps - indexes of angles of these points
		%minus side
		if j<clearance_points_for_check
			if i>2 
				rs = [image(i-2, end-clearance_points_for_check+j:end) image(i-1,1:j-1)];
				ps = [len-clearance_points_for_check+j:len, 1:j-1];
			else %before begin, need cycling repeat of bounding box
				rs = [image(i-1, end-clearance_points_for_check+j:end) image(i-1,1:j-1)];
				ps = [len-clearance_points_for_check+j:len, 1:j-1];
			end
		else	
			rs = image(i-1, 1+j-clearance_points_for_check:j-1);
			ps = 1+j-clearance_points_for_check:j-1;
		end
		%plus side
		if j+clearance_points_for_check>len
			rs = [rs image(i-1, j:end) image(i,1:clearance_points_for_check-(len-j+1))];
			ps = [ps, j:len, 1:clearance_points_for_check-(len-j+1)];
		else	
			rs = [rs image(i-1,j:j+clearance_points_for_check)];
			ps = [ps j:j+clearance_points_for_check];
		end
		
		%cosine theorem: x^2-2*ri*x*cos(p) + ri^2-clr^2 = 0;
		xs = zeros(length(rs),1);
		for k=1:length(rs)
			p = abs(angles(ps(k)) - angles(j));
			c = cos(p);
			D = (rs(k)*c)^2+clr^2-rs(k)^2;
			if D>=0
				roo = roots([1 -2*rs(k)*cos(p) rs(k)^2-clr^2]);
				xs(k) = min(roo);
                if xs(k)<0
                    xs(k) = max(roo);
                    %warning('negative root!!!');
                end
                if xs(k)<0
                    xs(k) = inf;
                    %warning('both negative roots!!!');
                end
			else
				xs(k) = inf;
            end
		end
		m = min(xs);
		if m==inf
			m = image(i-1,j)-clr;
		end
		image(i,j) = m;
    end
end

%% Magnetic moment
%ri - vector of current direction, p - phi, dl = r*dp
%m = 1/2 intagral [r,j] dV = 1/2*I integral [r,ri] r*dp
m1 = 0; %without I
m2 = 0;
ms1 = zeros(1,len);
ms2 = zeros(1,len);
for i=2:N+1
    for j=1:len
        if i==2 && j==1
            continue;
        end
        if j==1
            rim1 = image(i-1,end);
        else
            rim1 = image(i,j-1);
        end
        ri  = image(i,j);
        dr = abs(ri-rim1);
        
        % variant1
        %{
        AC = rim1*sinp;
        AB2 = ri^2+rim1^2-2*ri*rim1*cosp;
        cosg = (dr^2+AB2-AC^2)/(2*dr*sqrt(AB2));
        sina = sqrt(1-cosg^2);
        dm1 = 0.5e-6*dp*ri^2*sina; %0.5e-6 = 1/2 * 0.001 * 0.001. 1/2 from eqn, 0.001 - because r in mm
        m1 = m1 + dm1;
        ms1(j) = m1;
        %}
        
        %variant2
        % {
        AO = (ri-dr/2)*cosp2;
        AC = AO*tgp2;
        AB2 = ri^2+AO^2-2*ri*AO*cosp2;
        sina = (AC^2+AB2-(dr/2)^2) / (2*AC*sqrt(AB2));
        dm2 = 0.5e-6*dp*ri^2*sina; %0.5e-6 = 1/2 * 0.001 * 0.001. 
        %1/2 from eqn, 0.001 - because r in mm
        m2 = m2 +dm2;
        ms2(j) = m2;
        %}
    end
end

%% Resistance calculation
if use_calculate_resistance
% place code here

end

%% Visualizing
if use_render
    figure(1);
    clf;
    hold on;
    set(gcf, 'Color', 'w');
    [x,y] = pol2cart(angles,image(1,:));
    plot(x,y,'r');
    xs = [];
    ys = [];
    for i=2:N+1
        [x,y] = pol2cart(angles,image(i,:));
        %plot(x,y,'.-b');
        xs = [xs x];
        ys = [ys y];
    end
    plot(xs,ys,'.-b');
    xlim([-xsize/2*1.1 xsize/2*1.1]);
    ylim([-ysize/2*1.1 ysize/2*1.1]);
    %axis equal
end

%% Visualize resistance
if use_resistance_visualization
    
end

%% test ideal square
%{
ms_formula_angle = zeros(1,len);
ms_ideal = zeros(1,len);
for j=1:len
    as = pi/4:pi/4:2*pi;
    for i = 0:length(as)-1
        if angles(j) <= as(i+1)
            ms_ideal(j) = i/8*xsize^2 + 1/8*xsize^2*tan(angles(j)-pi/4*i);
            break
        end
    end
end
ms_ideal = 1e-6*ms_ideal;
m_ideal = ms_ideal(end);
m_formula_angle = 0;
for j=1:len
    if angles(j) < angle1
        ms_formula_angle(j) = abs(sin(pi/2-angles(j)));
        continue;
    end
    if angles(j) < angle2
        ms_formula_angle(j) = abs(cos(pi/2-angles(j)));
        continue;
    end
    if angles(j) < angle3
        ms_formula_angle(j) = abs(sin(pi/2-angles(j)));
        continue;
    end
    if angles(j) < angle4
        ms_formula_angle(j) = abs(cos(pi/2-angles(j)));
        continue;
    end
    ms_formula_angle(j) = abs(sin(pi/2-angles(j)));
end
for j=1:len
    ri = image(2,j);
    m_formula_angle_old = m_formula_angle;
    %m_formula_angle = m_formula_angle + 0.5e-6*dp*ri^2*ms_formula_angle(j);
    m_formula_angle = m_formula_angle + 0.5e-6*dp*ri^2;%*cos(angles(j));
    %m_formula_angle_old = 0.25e-6*xsize*ri*sin(angles(j)); %no error prodused
    ms_formula_angle(j) = m_formula_angle_old;
end
plot(angles*180/pi, ms1, 'r', angles*180/pi, ms2, 'b', ...
     angles*180/pi, ms_formula_angle, '*g', angles*180/pi, ms_ideal, 'k');
%}
