function plt_arrow(haxs, scale_fator, x, y, u, v)
%% Plot vector field on axes haxs as defined by locations (x, y) and vector
% components (u, v). scale_fator stretches or shrinks the length of the
% vectors.


%%

%
arrowsplt.scale_fator = scale_fator;


%%

%
arrowsplt.x = x;
arrowsplt.y = y;
arrowsplt.u = u;
arrowsplt.v = v;


%%

%
arrowsplt.uvmag = sqrt(arrowsplt.u.^2 + arrowsplt.v.^2);
arrowsplt.angrot = atan2(arrowsplt.v, arrowsplt.u);

%
arrowsplt.Narrows = length(arrowsplt.x);

%
arrowsplt.color = [1, 0, 1];


%% Define base arrow (that will be scaled and translated for plotting)

% %
% dx_arhead = 40;
% dy_arhead = 20;
% dy_arshaft = 10;
%
dx_arhead = 20;
dy_arhead = 10;
dy_arshaft = 7;
%
scale_head = 1;

%
x_shaft_arrow = [0, 1, NaN, NaN, NaN, 1, 0, 0];
y_shaft_arrow = [dy_arshaft, dy_arshaft, NaN, NaN, NaN];
%
y_shaft_arrow = [y_shaft_arrow, -fliplr(y_shaft_arrow(1:2)), dy_arshaft];


%
dx_head_arrow = [0, dx_arhead, 0];
dy_head_arrow = [(dy_arshaft+dy_arhead), 0, -(dy_arshaft+dy_arhead)];
%
dx_head_arrow = dx_head_arrow * scale_head;
dy_head_arrow = dy_head_arrow * scale_head;


%% Calculate coordinates of scaled and translated arrows

%
arrowsplt.x_array = NaN(8, arrowsplt.Narrows);
arrowsplt.y_array = arrowsplt.x_array;

%
for i = 1:arrowsplt.Narrows
    %
    x_aux = arrowsplt.uvmag(i) .* arrowsplt.scale_fator .* x_shaft_arrow;
    y_aux = y_shaft_arrow;
    %
    x_aux(3:5) = x_aux(2) + dx_head_arrow;
    y_aux(3:5) = dy_head_arrow;

    %
    ang_aux = arrowsplt.angrot(i);
    rot_aux = [cos(ang_aux), -sin(ang_aux); sin(ang_aux), cos(ang_aux)];

    %
    xyrot_aux = rot_aux * [x_aux(:).'; y_aux(:).'];

    %
    arrowsplt.x_array(:, i) = xyrot_aux(1, :);
    arrowsplt.y_array(:, i) = xyrot_aux(2, :);

    %
    arrowsplt.x_array(:, i) = arrowsplt.x_array(:, i) + arrowsplt.x(i);
    arrowsplt.y_array(:, i) = arrowsplt.y_array(:, i) + arrowsplt.y(i);

end



%% Plot vector field

%
fill(haxs, arrowsplt.x_array, arrowsplt.y_array, arrowsplt.color)


