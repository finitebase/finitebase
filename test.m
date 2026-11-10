%===================blade design tools===================

function C = bezier_cubic(t, P)
    B0 = (1-t).^3;
    B1 = 3*(1-t).^2 .* t;
    B2 = 3*(1-t) .* t.^2;
    B3 = t.^3;
    C = B0 * P(1,:) + B1 * P(2,:) + B2 * P(3,:) + B3 * P(4,:);
end

function dC = bezier_derivative(t, P)
    dB0 = -3*(1-t).^2;
    dB1 = 3*(1-t).^2 - 6*(1-t).*t;
    dB2 = 6*(1-t).*t - 3*t.^2;
    dB3 = 3*t.^2;
    dC = dB0 * P(1,:) + dB1 * P(2,:) + dB2 * P(3,:) + dB3 * P(4,:);
end

function T = thickness_nonzero_edges(relative_t, relative_t_maxium_thickness, relative_maxium_thickness, t_le, t_te)
    % THICKNESS_NONZERO_EDGES Compute thickness distribution with specified max thickness position
    n = 1.5;  % Shape parameter
    
    % Linear component connecting LE and TE thicknesses
    linear_component = t_le + (t_te - t_le) * relative_t;
    
    % Create shape function that peaks at t_maxium_thickness
    if relative_t_maxium_thickness > 0 && relative_t_maxium_thickness < 1
        trans = relative_t.^(log(0.5) / log(relative_t_maxium_thickness));
        %trans = ((0.5 - t_maxium_thickness)*t.^2 + (t_maxium_thickness^2 - 0.5).*t) / (t_maxium_thickness*(t_maxium_thickness-1));
        shape_function = sin(pi * trans).^n;
    else
        % Fallback to symmetric if position is invalid
        shape_function = sin(pi * relative_t).^n;
    end
    
    % Amplitude to achieve desired max thickness at t_maxium_thickness
    %linear_at_max = interpn(t, linear_component, t_maxium_thickness .* ones(size(linear_component)));
    linear_at_max = t_le + (t_te - t_le) * relative_t_maxium_thickness;
    A = relative_maxium_thickness - linear_at_max;
    
    % Final thickness distribution
    T = linear_component + A .* shape_function;
end

function point = compute_upper_point(t, P, relative_t_maxium_thickness, relative_maxium_thickness, t_le, t_te)
    
    camber = bezier_cubic(1 - t, P);
    deriv = bezier_derivative(1 - t, P);
    normal = [-deriv(2), deriv(1)];
    normal = normal / norm(normal);
    T = thickness_nonzero_edges(1 - t, relative_t_maxium_thickness, relative_maxium_thickness, t_le, t_te);
    point = camber + 0.5 * T * normal * (norm(bezier_cubic(1, P) - bezier_cubic(0, P)));
end

function point = compute_lower_point(t, P, relative_t_maxium_thickness, relative_maxium_thickness, t_le, t_te)
    camber = bezier_cubic(t, P);
    deriv = bezier_derivative(t, P);
    normal = [-deriv(2), deriv(1)];
    normal = normal / norm(normal);
    T = thickness_nonzero_edges(t, relative_t_maxium_thickness, relative_maxium_thickness, t_le, t_te);
    point = camber - 0.5 * T * normal * (norm(bezier_cubic(1, P) - bezier_cubic(0, P)));
end

function lower_func = get_lower_surface_func(P, relative_maxium_thickness, relative_t_maxium_thickness, t_le, t_te)
% GET_LOWER_SURFACE_FUNC Return parameterized function for lower surface
%   Input:  t ∈ [0,1], t=0 at trailing edge, t=1 at leading edge
%   Output: 1×2 coordinate [x, y]

    lower_func = @(t) compute_lower_point(t, P, relative_t_maxium_thickness, relative_maxium_thickness, t_le, t_te);
end

function trailing_func = get_trailing_edge_func(P, relative_maxium_thickness, relative_t_maxium_thickness, t_le, t_te)
% GET_TRAILING_EDGE_FUNC Return parameterized function for trailing edge arc
%   Input:  t ∈ [0,1], t=0 at upper surface TE, t=1 at lower surface TE
%   Output: 1×2 coordinate [x, y]

    % 获取后缘上下端点
    te_upper = compute_upper_point(0, P, relative_t_maxium_thickness, relative_maxium_thickness, t_le, t_te);
    te_lower = compute_lower_point(1, P, relative_t_maxium_thickness, relative_maxium_thickness, t_le, t_te);
    
    % 计算后缘切线方向
    eps_val = 1e-6;
    te_upper_near = compute_upper_point(eps_val, P, relative_t_maxium_thickness, relative_maxium_thickness, t_le, t_te);
    tangent_te = (te_upper - te_upper_near) / eps_val;
    
    % 计算后缘圆弧的圆心和半径
    A = [2*(te_upper(1)-te_lower(1)), 2*(te_upper(2)-te_lower(2)); 
         tangent_te(1), tangent_te(2)];
    b = [te_upper(1)^2 + te_upper(2)^2 - te_lower(1)^2 - te_lower(2)^2;
         te_upper(1)*tangent_te(1) + te_upper(2)*tangent_te(2)];
    te_center = (A \ b)';
    te_radius = norm(te_upper - te_center);
    
    % 计算角度范围
    angle_upper = atan2(te_upper(2) - te_center(2), te_upper(1) - te_center(1));
    angle_lower = atan2(te_lower(2) - te_center(2), te_lower(1) - te_center(1));
    if angle_lower > angle_upper
        angle_lower = angle_lower - 2*pi;
    end
    
    trailing_func = @(t) te_center + te_radius * ...
        [cos(angle_lower + t*(angle_upper - angle_lower)), ...
         sin(angle_lower + t*(angle_upper - angle_lower))];
end

function upper_func = get_upper_surface_func(P, relative_maxium_thickness, relative_t_maxium_thickness, t_le, t_te)
% GET_UPPER_SURFACE_FUNC Return parameterized function for upper surface
%   Input:  t ∈ [0,1], t=0 at leading edge, t=1 at trailing edge
%   Output: 1×2 coordinate [x, y]

    upper_func = @(t) compute_upper_point(t, P, relative_t_maxium_thickness, relative_maxium_thickness, t_le, t_te);
end

function leading_func = get_leading_edge_func(P, relative_maxium_thickness, relative_t_maxium_thickness, t_le, t_te)
% GET_LEADING_EDGE_FUNC Return parameterized function for leading edge arc
%   Input:  t ∈ [0,1], t=0 at lower surface LE, t=1 at upper surface LE
%   Output: 1×2 coordinate [x, y]

    % 获取前缘上下端点
    le_upper = compute_upper_point(1, P, relative_t_maxium_thickness, relative_maxium_thickness, t_le, t_te);
    le_lower = compute_lower_point(0, P, relative_t_maxium_thickness, relative_maxium_thickness, t_le, t_te);
    
    % 计算前缘切线方向
    eps_val = 1e-6;
    le_upper_near = compute_upper_point(1 - eps_val, P, relative_t_maxium_thickness, relative_maxium_thickness, t_le, t_te);
    tangent_le = (le_upper_near - le_upper) / eps_val;
    
    % 计算前缘圆弧的圆心和半径
    A = [2*(le_upper(1)-le_lower(1)), 2*(le_upper(2)-le_lower(2)); 
         tangent_le(1), tangent_le(2)];
    b = [le_upper(1)^2 + le_upper(2)^2 - le_lower(1)^2 - le_lower(2)^2;
         le_upper(1)*tangent_le(1) + le_upper(2)*tangent_le(2)];
    le_center = (A \ b)';
    le_radius = norm(le_upper - le_center);
    
    % 计算角度范围
    angle_upper = atan2(le_upper(2) - le_center(2), le_upper(1) - le_center(1));
    angle_lower = atan2(le_lower(2) - le_center(2), le_lower(1) - le_center(1));
    if angle_lower < angle_upper
        angle_lower = angle_lower + 2*pi;
    end
    
    leading_func = @(t) le_center + le_radius * ...
        [cos(angle_upper + t*(angle_lower - angle_upper)), ...
         sin(angle_upper + t*(angle_lower - angle_upper))];
end

function L = function_limit(f, t, initial_dt)
    dt  = initial_dt;
    value = f(t + dt);
    dt = dt / 2;
    new_value = f(t + dt);
    osc = norm(new_value - value);
    dt = dt / 2;
    new_value = f(t + dt);
    osc_ = 2 * osc;
    %osc = norm(new_value - val);
    
    while(abs(osc_ / osc - 1) < 0.2 && osc > 0)
        dt = dt / 2;
        new_value = f(t + dt);
        osc_ = osc;
        osc = norm(new_value - value);
    end
    L = new_value;
end

function tangent_vector = curve_tangent_vector(f,t)
    dft = @(x) (f(t + x) - f(t)) / x;
    
    if t == 0
        tangent_vector = function_limit(dft, 0,1e-4);
    elseif t == 1
        tangent_vector = function_limit(dft, 0, -1e-4);
    else
        tangent_vector = 0.5 * function_limit(dft, 0, min(0.5 * abs(1 - t), 1e-4)) + 0.5 * function_limit(dft, 0, -min(0.5 * abs(t), 1e-4));
    end
    if norm(tangent_vector) > 0
        tangent_vector = tangent_vector / norm(tangent_vector);
    end
    
end

%===================mesh tools===================

function ij = zigzag(m, n, k)
% ZIGZAG Maps 1D index k to 2D coordinates (i,j) using zigzag scan order.
%   ij = ZIGZAG(m, n, k) returns an N-by-2 matrix of [i, j] coordinates
%   corresponding to the linear indices in k for an m-by-n matrix
%   traversed in zigzag order.
%
%   Inputs:
%       m  - number of rows (positive integer)
%       n  - number of columns (positive integer)
%       k  - vector or array of linear indices (values in [1, m*n])
%            Internally reshaped to a column vector.
%
%   Output:
%       ij - N-by-2 matrix where each row is [i, j] for the corresponding
%            element of k(:).  Size matches the number of elements in k.
%
%   Example:
%       ij = zigzag(3, 4, [1, 6, 12]);
%       % Returns coordinates for indices 1, 6, and 12 in zigzag order.

    p = min(m, n);          % shorter dimension
    q = max(m, n);          % longer dimension
    mn = m * n;             % total number of elements
    S1 = p*(p+1)/2;         % cumulative sum at the end of the growing triangle
    S2 = S1 + (q-p)*p;      % cumulative sum at the end of the constant middle section

    % --- Force input to column vector ---
    % Reshape k into a column vector regardless of original orientation.
    % This ensures consistent behavior for row vectors, column vectors,
    % and higher-dimensional arrays.
    k = k(:);
    N = length(k);          % number of input indices

    % Preallocate output arrays as column vectors
    d    = zeros(N, 1);     % diagonal index: d = i + j - 1  (1 to m+n-1)
    pos  = zeros(N, 1);     % position along the d-th diagonal

    % --- Region 1: growing triangle (d <= p) ---
    % First p diagonals where diagonal length increases: 1, 2, ..., p
    % Cumulative sum before diagonal d: (d-1)*d/2
    mask1 = (k <= S1);
    if any(mask1)
        k1 = k(mask1);
        d(mask1)   = ceil((sqrt(8*k1 + 1) - 1) / 2);
        prev       = (d(mask1) - 1) .* d(mask1) / 2;
        pos(mask1) = k1 - prev;
    end

    % --- Region 2: constant middle section (p < d <= q) ---
    % Middle diagonals where length is constant and equal to p
    % Only exists if q > p (non-square matrices)
    mask2 = (k > S1) & (k <= S2);
    if any(mask2)
        k2 = k(mask2);
        d(mask2)   = p + ceil((k2 - S1) / p);
        prev       = S1 + (d(mask2) - p - 1) * p;
        pos(mask2) = k2 - prev;
    end

    % --- Region 3: shrinking triangle (d > q) ---
    % Last p-1 diagonals where length decreases: p-1, p-2, ..., 1
    % Work backwards: reversed index k_rev runs through shrinking triangle
    mask3 = (k > S2);
    if any(mask3)
        k3 = k(mask3);
        k_rev     = mn - k3 + 1;              % reversed index counting from end
        t         = ceil((sqrt(8*k_rev + 1) - 1) / 2);  % tail position
        d(mask3)  = m + n - t;                % diagonal from the end
        prev      = mn - t .* (t + 1) / 2;    % cumulative sum at start of diagonal
        pos(mask3)= k3 - prev;
    end

    % --- Convert diagonal index and position to 2D coordinates ---
    i = zeros(N, 1);     % row indices (column vector)
    j = zeros(N, 1);     % column indices (column vector)

    % Odd diagonals: traverse top-to-bottom
    % Starting row: max(1, d-n+1)  — ensures we stay within bounds
    odd = (mod(d, 2) == 1);
    if any(odd)
        i_min    = max(1, d(odd) - n + 1);
        i(odd)   = i_min + pos(odd) - 1;
        j(odd)   = d(odd) - i(odd) + 1;
    end

    % Even diagonals: traverse bottom-to-top
    % Starting row: min(m, d)  — ensures we stay within bounds
    even = (mod(d, 2) == 0);
    if any(even)
        i_max    = min(m, d(even));
        i(even)  = i_max - pos(even) + 1;
        j(even)  = d(even) - i(even) + 1;
    end

    % Return as N-by-2 matrix
    ij = [i, j];
end

function k = dezigzag(m, n, i, j)
    % DEZIGZAG Maps 2D coordinates (i,j) to 1D index k using zigzag scan order.
    %   k = DEZIGZAG(m, n, i, j) returns the linear index for coordinates
    %   (i,j) in an m-by-n matrix traversed in zigzag order.
    %
    %   Inputs:
        %       m  - number of rows (positive integer)
    %       n  - number of columns (positive integer)
        %       i  - vector or array of row indices (values in [1, m])
    %            Internally reshaped to a column vector.
    %       j  - vector or array of column indices (values in [1, n])
    %            Internally reshaped to a column vector.  Must have the
    %            same number of elements as i.
    %
        %   Output:
    %       k  - column vector of linear indices in [1, m*n], same size as i(:).
    %
    %   Example:
    %       k = dezigzag(3, 4, [1, 2, 3], [1, 2, 1]);
    %       % Returns zigzag linear indices for coordinates (1,1), (2,2), (3,1).


    oob = (i < 1) | (i > m) | (j < 1) | (j > n);
   

    p = min(m, n);          % shorter dimension
    q = max(m, n);          % longer dimension
    mn = m * n;             % total number of elements

    % --- Force inputs to column vectors ---
    % Reshape i and j into column vectors regardless of original orientation.
    % This guarantees element-wise correspondence and consistent output shape.
    i = i(:);
    j = j(:);

    % Diagonal index: d = i + j - 1, ranges from 1 to m+n-1
    d = i + j - 1;

    % Number of input coordinates
    N = length(d);

    % Preallocate cumulative sum of previous diagonals
    prev = zeros(N, 1);     % S(d-1): total elements before diagonal d

    % --- Compute cumulative sum of elements in previous diagonals ---

    % Region 1: growing triangle (d <= p)
    % Diagonals 1 through p: lengths are 1, 2, ..., p
    % Cumulative sum before diagonal d: (d-1)*d / 2
    mask1 = (d <= p);
    if any(mask1)
        d1 = d(mask1);
        prev(mask1) = (d1 - 1) .* d1 / 2;
    end

    % Region 2: constant middle section (p < d <= q)
    % Diagonals p+1 through q: all have length p
    % Only exists for non-square matrices (q > p)
    mask2 = (d > p) & (d <= q);
    if any(mask2)
        d2 = d(mask2);
        prev(mask2) = p*(p+1)/2 + (d2 - p - 1) * p;
    end

    % Region 3: shrinking triangle (d > q)
    % Final p-1 diagonals: lengths are p-1, p-2, ..., 1
    % Use symmetry: remaining diagonals form a triangle of size
    % t = m + n - d, which contains t*(t+1)/2 elements from the end
    mask3 = (d > q);
    if any(mask3)
        d3 = d(mask3);
        t = m + n - d3;                      % tail length from the end
        prev(mask3) = mn - t .* (t + 1) / 2; % subtract tail from total
    end

    % --- Compute position along the current diagonal ---
    pos = zeros(N, 1);     % r: 1-based position within diagonal d

    % Odd diagonals: traversed top-to-bottom
    % Position determined by how far i is from the starting row
    odd = (mod(d, 2) == 1);
    if any(odd)
        i_min    = max(1, d(odd) - n + 1);   % starting row of diagonal
        pos(odd) = i(odd) - i_min + 1;        % position = row offset + 1
    end

    % Even diagonals: traversed bottom-to-top
    % Position determined by how far i is from the ending row
    even = (mod(d, 2) == 0);
    if any(even)
        i_max     = min(m, d(even));          % ending row of diagonal
        pos(even) = i_max - i(even) + 1;       % position = offset from bottom + 1
    end

    % Final 1D index = cumulative sum before diagonal + position in diagonal
    k = prev + pos;
    k(oob) = -1;
end

function mesh = block_mesh(C1, C2, C3, C4, n_x, n_y, x_control, y_control)
    % this function generate a "block mesh", which is a different type from the 
    % mesh generated by the latter functions, the later one is generated by "stich"
    % them together.
    %
    %
    %              C4
    %        *--<----------*            -> y
    %        |             |
    %        |             Λ
    %    C1  |             |   C3
    %        |             |
    %        V             |
    %        |             |
    %        *---------->--*
    %              C2
    %
    %        |
    %        v
    %        x
    %
    %
    %
    %
    % The last 2 argument are optional ones, either put none of them, or
    % both of them, in format of [relative_first_offset, growth_rate, direction]
    % the type should be respectively real number, real number greater then
    % 1, and {-1,1}



    x_list = transpose(linspace(0,1,n_x));
    y_list = transpose(linspace(0,1,n_y));
    
    
    

    if nargin > 6
        
        if ~isempty(x_control)
            relative_first_offset_x = x_control(1);
            growth_rate_x = x_control(2);
            direction_x = x_control(3);



            x_list = zeros([n_x,1]);
            x_list(n_x) = 1;
        
            current_position = 0;
            last_stet_length = relative_first_offset_x;
            growth_rate = growth_rate_x;


            count = 1;
            while current_position + growth_rate * last_stet_length < 1 && count < n_x && (1 - current_position) / (n_x - count) >= last_stet_length
                current_position = current_position + last_stet_length;
                last_stet_length = growth_rate * last_stet_length;
                x_list(count + 1) = current_position;
                count = count + 1;
            end
            
            average_length = (1 - current_position) / (n_x - count);
            for i = count + 1 : n_x - 1
                x_list(i) = (i - count) * average_length + x_list(count);
            end
            if direction_x == -1
                x_list = flip(ones(n_x,1) - x_list);
            end
        end


    

        if ~isempty(y_control)
            relative_first_offset_y = y_control(1);
            growth_rate_y = y_control(2);
            direction_y = y_control(3);


            y_list = zeros([n_y,1]);
            y_list(n_y) = 1;
        
            current_position = 0;
            last_stet_length = relative_first_offset_y;
            growth_rate = growth_rate_y;


            count = 1;
            while current_position + growth_rate * last_stet_length < 1 && count < n_y && (1 - current_position) / (n_y - count) >= last_stet_length
                current_position = current_position + last_stet_length;
                last_stet_length = growth_rate * last_stet_length;
                y_list(count + 1) = current_position;
                count = count + 1;
            end
            
            
            average_length = (1 - current_position) / (n_y - count);
            for i = count + 1 : n_y - 1
                y_list(i) = (i - count) * average_length + y_list(count);
            end
            if direction_y == -1
                y_list = flip(ones(n_y,1) - y_list);
            end
        end

    




    end




    nodes = zeros([n_x * n_y,2]);
    cell_node_pointer = zeros([n_x * n_y + 1 - n_x - n_y, 4]);
    cell_neighbour_pointer = zeros([n_x * n_y + 1 - n_x - n_y, 4]);

    for i = 1:n_x
        for j = 1:n_y
            % this part is to interpolate the curve approperly, basically
            % using linear interpolation to get curve shape, and use a
            % orthogonal transformation to adjust the extreme points.
            
            % this algorthim is not invariant by rotation, if you turn the
            % block around (change the name of curve, for example  C1, C2, C3, C4 -> C2, C3, C4, C1)
            % the mesh will be totally different, not only rotated.


            k = dezigzag(n_x, n_y, i, j);


            x = x_list(i);
            y = y_list(j);
            %if nargin > 6 && ~isempty(x_control) && ~isempty(y_control)
                %x = x_list(i);
                %y = y_list(j);
            %else
                %x = (i - 1) / (n_x - 1);
                %y = (j - 1) / (n_y - 1);
            %end
            



            nodes(k,:) = (1 - x) * C4(1 - y) + x * C2(y);
            C = (1 - x) * C4(1 - y) + x * C2(y);
            % known values, by row vector
            A = (1 - x) * C4(1 - 0) + x * C2(0);
            B = (1 - x) * C4(1 - 1) + x * C2(1);
            A_prime = C1(x);
            B_prime = C3(1 - x);
            

            % 1. compute the variance vector
            u = B(1) - A(1);
            v = B(2) - A(2);
            u_prime = B_prime(1) - A_prime(1);
            v_prime = B_prime(2) - A_prime(2);

            % 2. compute the rotatiion matrix, avoid using arcsin or arccos
            % to have better precision
            denom = u^2 + v^2;
            c = (u*u_prime + v*v_prime) / denom;  % cos(θ)
            s = (u*v_prime - v*u_prime) / denom;  % sin(θ)

            % 3. construct the matrix
            R = [c, -s; s, c];

            % 4. translation
            A_rot = (A * transpose(R));
            T = A_prime - A_rot;
    
            % 5. transformation of C
            C_prime = (C * transpose(R)) + T;
            nodes(k,:) = C_prime;
        end
    end
    
    for i = 1:n_x - 1
        for j = 1:n_y - 1
            k = dezigzag(n_x - 1, n_y - 1, i, j);
            cell_node_pointer(k,:) = [dezigzag(n_x, n_y, i, j), dezigzag(n_x, n_y, i + 1, j), dezigzag(n_x, n_y, i + 1, j + 1), dezigzag(n_x, n_y, i, j + 1)];


            %A = dezigzag(n_13, n_24, i + 1, j);
            %B = dezigzag(n_13, n_24, i, j + 1);
            %C = dezigzag(n_13, n_24, i - 1, j);
            %D = dezigzag(n_13, n_24, i, j - 1);
            cell_neighbour_pointer(k,:) = [dezigzag(n_x - 1, n_y - 1, i, j - 1), dezigzag(n_x - 1, n_y - 1, i + 1, j), dezigzag(n_x - 1, n_y - 1, i, j + 1), dezigzag(n_x - 1, n_y - 1, i - 1, j)];
        end
    end


    % return the mesh struct
    mesh = struct();
    mesh.nodes = nodes;
    mesh.cell_node_pointer = cell_node_pointer;
    mesh.cell_neighbour_pointer = cell_neighbour_pointer;
    mesh.m = n_x - 1;
    mesh.n = n_y - 1;

end

function mesh = stich_block_mesh(block_mesh_list, stich_list, boundry_list)
    % this function stich up those block meshes, and give a boundry list
    %
    % input:
    %   block_mesh_list : block_mesh struct list, which shoulld be generated by the former function
    %   stich_list      : N×4，for each row [bA, eA, bB, eB] they are block index of side one, its edge index,
    %                     block index of the other side one, its edge index, respectively
    %                     edge index means: 1=C1(A-B), 2=C2(B-C), 3=C3(D-C), 4=C4(A-D)
    %   boundary_list   : M×3，each row [block_index, edge_index, boundary_index], it is a discription of boundry edge
    %                     （optional, put an [] if not necessary）
    %
    % edge agreement: 1=C1(A-B), 2=C2(B-C), 3=C3(D-C), 4=C4(A-D)
    %
    % output mesh:
    %
    %   mesh.nodes                  : number_of_cells x 2 table, node
    %                                   position
    %
    %   mesh.cell_node_pointer      : number_of_cells x 4 table, refer to
    %                                   node index
    %
    %   mesh.cell_neighbour_pointer : number_of_cells x 4 table, refer to
    %                                   the index of neighbourhood cells, in same order as the node
    %
    %   mesh.boundary_faces         : K×3，for each row [node_a, node_b, owner_cell_global]
    %                                 （golobal node index + owner cell index）
    %   mesh.boundary_index         : K×1，owner boundry group's boundary_index for each boundry faces 
    %
    %   mesh.boundary_groups        : cell list，sorted or classification by boundary_index for output into OpenFOAM afterwards
    %                                   notice that the index in here should be
    %                                   from mesh.boundary_index, the rest part of
    %                                   this cell are empty(undefined). In those
    %                                   defined part, it is a table of face
    %                                   information, in each row it represent a
    %                                   face, as 
    %                                   [node_a, node_b, owner_cell_global], 
    %                                   or say, it is a regroup
    %                                   of mesh.boundary_faces
    %
    %
    %                               
    %
    %
    %
    %






    function ids = edge_local_ids(nx, ny, edge)
    % edge 约定: 1=左(j=1), 2=下(i=nx), 3=右(j=ny), 4=上(i=1)
    % 标准方向: 沿 i 递增 (左/右) 或沿 j 递增 (下/上)
    switch edge
        case 1   % 左: j=1, i=1..nx
            ids = dezigzag(nx, ny, (1:nx)',         ones(nx, 1));
        case 2   % 下: i=nx, j=1..ny
            ids = dezigzag(nx, ny, nx*ones(ny,1),  (1:ny)');
        case 3   % 右: j=ny, i=1..nx
            ids = dezigzag(nx, ny, (1:nx)',         ny*ones(nx, 1));
        case 4   % 上: i=1, j=1..ny
            ids = dezigzag(nx, ny, ones(ny,1),     (1:ny)');
        otherwise
            error('Unknown edge id: %d', edge);
        end
    end
    
    function [cell_ids, nbr_idx] = edge_cell_ids(bm, edge, c_off)
        % edge 约定: 1=左, 2=下, 3=右, 4=上
        % cell_neighbour_pointer 列号 == edge 编号
        mc = bm.m;   nc = bm.n;
        switch edge
            case 1   % 左: j=1
                cell_ids = dezigzag(mc, nc, (1:mc)',         ones(mc,1));
            case 2   % 下: i=mc
                cell_ids = dezigzag(mc, nc, mc*ones(nc,1),  (1:nc)');
            case 3   % 右: j=nc
                cell_ids = dezigzag(mc, nc, (1:mc)',         nc*ones(mc,1));
            case 4   % 上: i=1
                cell_ids = dezigzag(mc, nc, ones(nc,1),     (1:nc)');
        end
        cell_ids = cell_ids + c_off;
        nbr_idx  = edge;
    end

    function p = compute_local_coord(bm, edge, which)
        % edge 约定: 1=左(j=1), 2=下(i=nx), 3=右(j=ny), 4=上(i=1)
        % 标准方向起点 = (1, 1) 或 (1, *) 或 (*, 1) 中靠近原点那个
        nx = bm.m + 1;   ny = bm.n + 1;
        switch edge
            case 1   % 左: i=1..nx, j=1
                if strcmp(which, 'first'),  k_ = dezigzag(nx, ny, 1,  1);
                else,                       k_ = dezigzag(nx, ny, nx, 1); end
            case 2   % 下: i=nx, j=1..ny
                if strcmp(which, 'first'),  k_ = dezigzag(nx, ny, nx, 1);
                else,                       k_ = dezigzag(nx, ny, nx, ny); end
            case 3   % 右: i=1..nx, j=ny
                if strcmp(which, 'first'),  k_ = dezigzag(nx, ny, 1,  ny);
                else,                       k_ = dezigzag(nx, ny, nx, ny); end
            case 4   % 上: i=1, j=1..ny
                if strcmp(which, 'first'),  k_ = dezigzag(nx, ny, 1,  1);
                else,                       k_ = dezigzag(nx, ny, 1,  ny); end
        end
        p = bm.nodes(k_, :);
    end







    if nargin < 3
        boundry_list = zeros(0, 3);
    end

    block_n = length(block_mesh_list);

    % ---------- 1. stich 邻接表 ----------
    stich_of = cell(block_n, 4);
    for s = 1:size(stich_list, 1)
        bA = stich_list(s, 1);  eA = stich_list(s, 2);
        bB = stich_list(s, 3);  eB = stich_list(s, 4);
        stich_of{bA, eA} = [bB, eB, NaN];
        stich_of{bB, eB} = [bA, eA, NaN];
    end

    % ---------- 2. 预估容量 ----------
    max_nodes = 0;
    total_cells = 0;
    cell_offset = zeros(block_n + 1, 1);
    for b = 1:block_n
        bm = block_mesh_list(b);
        nx = bm.m + 1;   ny = bm.n + 1;
        max_nodes   = max_nodes + nx * ny;
        total_cells = total_cells + bm.m * bm.n;
        cell_offset(b+1) = total_cells;
    end

    nodes = zeros(max_nodes, 2);
    node_count = 0;
    cell_node_pointer      = zeros(total_cells, 4);
    cell_neighbour_pointer = zeros(total_cells, 4);
    edge_node_table = cell(block_n, 4);
    added = false(block_n, 1);

    % ---------- 3. 依次加入每个 block ----------
    for b = 1:block_n
        bm = block_mesh_list(b);
        nx = bm.m + 1;   ny = bm.n + 1;
        nn = nx * ny;
        mc = bm.m;       nc = bm.n;

        local_to_global = zeros(nn, 1);

        % --- 3a. 复用已加入 block 缝合边的节点 ---
        for e = 1:4
            info = stich_of{b, e};
            if isempty(info), continue; end
            other_b = info(1);
            if ~added(other_b), continue; end
            other_e = info(2);

            local_ids  = edge_local_ids(nx, ny, e);
            other_glob = edge_node_table{other_b, other_e};

            p_my1 = compute_local_coord(bm, e, 'first');
            p_my2 = compute_local_coord(bm, e, 'last');
            p_ot1 = nodes(other_glob(1),   :);
            p_ot2 = nodes(other_glob(end), :);
            if (norm(p_my1 - p_ot2) + norm(p_my2 - p_ot1)) < ...
               (norm(p_my1 - p_ot1) + norm(p_my2 - p_ot2))
                other_glob = flipud(other_glob);
            end

            if (size(other_glob) ~= size(local_ids))
                error("segment number does not match");
            end
            local_to_global(local_ids) = other_glob;
            stich_of{b, e}(3) = 0;
        end

        % --- 3b. 新增剩余节点 ---
        new_mask = (local_to_global == 0);
        n_new = sum(new_mask);
        new_ids = node_count + (1:n_new)';
        local_to_global(new_mask) = new_ids;
        nodes(new_ids, :) = bm.nodes(new_mask, :);
        node_count = node_count + n_new;

        % --- 3c. 填 edge_node_table ---
        for e = 1:4
            local_ids = edge_local_ids(nx, ny, e);
            edge_node_table{b, e} = local_to_global(local_ids);
        end

        % --- 3d. 转换到全局编号 ---
        c_off = cell_offset(b);
        nc_block = mc * nc;

        cnp = bm.cell_node_pointer;
        cnp_global = zeros(size(cnp));
        valid = cnp > 0;
        cnp_global(valid) = local_to_global(cnp(valid));
        cell_node_pointer(c_off + (1:nc_block), :) = cnp_global;

        cnb = bm.cell_neighbour_pointer;
        cnb_global = -ones(size(cnb));
        valid = cnb > 0;
        cnb_global(valid) = cnb(valid) + c_off;
        cell_neighbour_pointer(c_off + (1:nc_block), :) = cnb_global;

        added(b) = true;
    end

    % ---------- 4. 跨 block 邻居更新 ----------
    for s = 1:size(stich_list, 1)
        bA = stich_list(s, 1);  eA = stich_list(s, 2);
        bB = stich_list(s, 3);  eB = stich_list(s, 4);

        [cellsA, nbrA_idx] = edge_cell_ids(block_mesh_list(bA), eA, cell_offset(bA));
        [cellsB, nbrB_idx] = edge_cell_ids(block_mesh_list(bB), eB, cell_offset(bB));

        nodesA = edge_node_table{bA, eA};
        nodesB = edge_node_table{bB, eB};
        if nodesA(1) == nodesB(end) && nodesA(end) == nodesB(1)
            cellsB = flipud(cellsB);
        end

        cell_neighbour_pointer(cellsA, nbrA_idx) = cellsB;
        cell_neighbour_pointer(cellsB, nbrB_idx) = cellsA;
    end

    % ---------- 5. 处理边界标记 ----------
    % 对 boundary_list 中每条 [block, edge, boundary_index]，
    % 把该边上每一个 cell-edge 转换为一个边界面（两节点 + owner cell）。
    boundary_faces = zeros(0, 3);
    boundary_index = zeros(0, 1);
    boundary_type = [];

    for k_ = 1:size(boundry_list, 1)
        b   = boundry_list(k_, 1);
        e   = boundry_list(k_, 2);
        bindex = boundry_list(k_, 3);
        

        bm = block_mesh_list(b);
        [cell_ids, ~] = edge_cell_ids(bm, e, cell_offset(b));
        edge_nodes = edge_node_table{b, e};   % length = mc+1 或 nc+1

        % 沿边的每段 (i -> i+1) 是一个边界面，owner 就是对应 cell
        for i = 1:length(cell_ids)
            face = [edge_nodes(i), edge_nodes(i+1), cell_ids(i)];
            boundary_faces(end+1, :) = face; %#ok<AGROW>
            boundary_index(end+1, 1) = bindex; %#ok<AGROW>
            
        end

        % 同时把这些 cell 在该 edge 方向的邻居指针明确设为 -btype
        % （保留 -1 表示未标记边界；负数编码 boundary_index 便于求解器识别）
        %cell_neighbour_pointer(cell_ids, e) = -bindex;
    end

    % 按 boundary_index 分组（方便 OpenFOAM patch 输出）
    boundary_groups = struct();
    if ~isempty(boundary_index)
        uniq = unique(boundary_index);
        boundary_groups.index   = uniq;
        boundary_groups.faces = cell(numel(uniq), 1);
        
        for q = 1:numel(uniq)
            boundary_groups.faces{q} = boundary_faces(boundary_index == uniq(q), :);
        end
    else
        boundary_groups.index   = [];
        boundary_groups.faces = {};
    end

    % ---------- 6. 截断并输出 ----------
    nodes = nodes(1:node_count, :);

    mesh = struct();
    mesh.nodes                  = nodes;
    mesh.cell_node_pointer      = cell_node_pointer;
    mesh.cell_neighbour_pointer = cell_neighbour_pointer;
    mesh.boundary_faces         = boundary_faces;   % [n1, n2, owner_cell]
    mesh.boundary_index         = boundary_index;   % boundary_index of those faces, corresponding to the list boundary_faces
    mesh.boundary_groups        = boundary_groups;  % 按类型分组
end

function f = segment(A, B)
    A = reshape(A, 1, []);
    B = reshape(B, 1, []);
    f = @(t) (1 - t(:)) .* A + t(:) .* B;
end

function block = cluster_block_cell_initialize(A,B,C,D,C1,C2,C3,C4)
    block = struct();
    block.A = A;
    block.B = B;
    block.C = C;
    block.D = D;
    block.C1 = C1;
    block.C2 = C2;
    block.C3 = C3;
    block.C4 = C4;
end

function parallel_class = compute_block_cluster_parallel_class(block_cluster)
    % compute the equivalence class generated by thoes blocks, which is the
    % parallel edge of quadrilateeral
    % for each block, which is a quadrilateeral ABCD，they are stored as
    % index in vertice_list (as [u,v])
    % we defined those parallel edges, orientanted：
    %   (A->B) // (D->C)      —— which is C1 C3 are of oppsite direction
    %   (A->D) // (B->C)      —— which is C4 C2 are of oppsite direction
    % if 2 blocks share one orientanted edge (u,v), we say it is
    % equivalence to itself
    % we only regist a edge once, either itself or its opposite, to avoid
    % messed up
    %
    % output:
    %   parallel_class : cell list，each element is a Nx2 matrix，
    %                    each row is a orientanted edge [u, v]。
    %                    each Nx2 matrix is a equivalence class

    block_list   = block_cluster.block_list;
    n_blocks     = length(block_list);

    % ---------- 1. 收集所有出现过的有向边 ----------
    % 每个块贡献 4 条有向边: A->B, B->C, C->D, D->A
    % 用 map: key = "u,v" 字符串, value = 该有向边在 edge_list 中的索引
    edge_map  = containers.Map('KeyType', 'char', 'ValueType', 'int32');
    edge_list = zeros(0, 2);   % 第 i 行 = [u, v]
    n_edges   = 0;

    % 辅助函数：登记一条无向边，返回它在 edge_list 中的索引
    function [idx,orientation] = register(u, v)
        orientation = 1;
        key = sprintf('%d,%d', u, v);
        key_ = sprintf('%d,%d', v, u);
        if isKey(edge_map, key)
            idx = edge_map(key);
        elseif isKey(edge_map, key_)
            idx = edge_map(key_);
            orientation = -1;
        else
            n_edges = n_edges + 1;
            edge_list(n_edges, :) = [u, v];
            edge_map(key) = int32(n_edges);
            idx = n_edges;
            orientation = 1;
        end
    end

    % 记录每个块的 4 条有向边在 edge_list 中的索引，用于后面建立平行关系
    block_edges = zeros(n_blocks, 4,2);  % 列对应 [AB, BC, CD, DA]
    for b = 1:n_blocks
        blk = block_list(b);
        A = blk.A;  B = blk.B;  C = blk.C;  D = blk.D;
        block_edges(b, 1, :) = register(A, B);   % AB
        block_edges(b, 2, :) = register(B, C);   % BC
        block_edges(b, 3, :) = register(C, D);   % CD
        block_edges(b, 4, :) = register(D, A);   % DA
    end

    % ---------- 2. 并查集初始化 ----------
    parent = (1:n_edges)';

    function r = find_root(x)
        % 路径压缩 find
        r = x;
        orientation = 1;
        while abs(parent(r)) ~= r
            if parent(r) < 0
                orientation = -orientation;
            end
            r = abs(parent(r));
        end
        r = orientation * r;
    end

    function union(x, y, orientation)
        ry = find_root(y);
        rx = find_root(x);
        if rx ~= ry && rx ~= -ry
            parent(abs(x)) = x / abs(rx) * orientation * ry;
        end


        for update = 1:n_edges
            if update == x
                continue
            elseif parent(update) == rx
                parent(update) = orientation * ry;
            elseif parent(update) == -rx
                parent(update) = -orientation * ry;
            end
        end
    end

    % ---------- 3. 在每个块内合并平行有向边 ----------
    % 块内平行关系（同向）:
    %   AB  ~  DC  (= reverse of CD)
    %   AD  (= reverse of DA)  ~  BC
    % 所以我们需要 CD 的反向边和 DA 的反向边也作为可被合并的对象。
    for b = 1:n_blocks

        blk = block_list(b);
        A = blk.A;  B = blk.B;  C = blk.C;  D = blk.D;
        [e1, orientation1] = register(A, B);                 % A -> B
        [e2, orientation2] = register(D, C);                    % D -> C （CD 的反向，单独登记）
        union(e1, e2, orientation2 * orientation1);

        [e1, orientation1] = register(B, C);                 % A -> B
        [e2, orientation2] = register(A, D);                    % D -> C （CD 的反向，单独登记）
        union(e1, e2, orientation2 * orientation1);
    end

    % ---------- 4. 按根节点分组 ----------
    roots = zeros(n_edges, 1);
    for k = 1:n_edges
        roots(k) = abs(parent(k));
        if parent(k) < 0
            edge_list(k,:) = flip(edge_list(k,:));
        end
    end
    [unique_roots, ~, group_id] = unique(roots);

    parallel_class = cell(numel(unique_roots), 1);
    for g = 1:numel(unique_roots)
        parallel_class{g} = edge_list(group_id == g, :);
    end
end

function stich_list = generate_stich_list(cluster)
    % GENERATE_STICH_LIST 根据 cluster 中各 block 的拓扑生成缝合表
    %
    % 输入:
    %   cluster.block_list(i).A/B/C/D : 顶点在 vertice_list 中的 index
    %
    % 输出:
    %   stich_list : N×4，每行 [bA, eA, bB, eB]
    %                edge 编号: 1=C1(A-B), 2=C2(B-C), 3=C3(D-C), 4=C4(A-D)
    %
    % 规则: 两个 block 若共享同一条无向边 {u,v}，则在 stich_list 中登记一次。

    block_list = cluster.block_list;
    n_blocks   = length(block_list);

    % 每个 block 的 4 条边（顶点 index 对），按 edge 编号顺序:
    %   edge 1 (C1): A-B
    %   edge 2 (C2): B-C
    %   edge 3 (C3): D-C
    %   edge 4 (C4): A-D
    edge_def = [1 2;   % AB
                2 3;   % BC
                4 3;   % DC
                1 4];  % AD
    vert_names = {'A','B','C','D'};

    % 收集所有 (block_idx, edge_idx, key) 三元组
    % key 用 "min,max" 字符串保证无向边唯一
    edge_map = containers.Map('KeyType','char','ValueType','any');

    for b = 1:n_blocks
        blk = block_list(b);
        verts = [blk.A, blk.B, blk.C, blk.D];
        for e = 1:4
            u = verts(edge_def(e,1));
            v = verts(edge_def(e,2));
            key = sprintf('%d,%d', min(u,v), max(u,v));
            if isKey(edge_map, key)
                edge_map(key) = [edge_map(key); b, e];
            else
                edge_map(key) = [b, e];
            end
        end
    end

    % 对每条出现 ≥2 次的无向边，两两配对生成 stich
    stich_list = zeros(0, 4);
    keys = edge_map.keys;
    for k = 1:numel(keys)
        entries = edge_map(keys{k});   % 每行 [block, edge]
        if size(entries,1) < 2
            continue;
        end
        % 一般几何中一条边最多被两个 block 共享，但也兼容多于 2 的情况
        for i = 1:size(entries,1)
            for j = i+1:size(entries,1) 
                stich_list(end+1, :) = ...
                    [entries(i,1), entries(i,2), entries(j,1), entries(j,2)]; %#ok<AGROW>
            end
        end
    end
end

function cluster = block_cluster_initialize(block_list, vertice_list)
    cluster = struct();
    cluster.block_list = block_list;
    cluster.vertice_list = vertice_list;
    cluster.parallel_class = [];
    cluster.division_class = [];
    cluster.stich_list = [];
    cluster.parallel_class = compute_block_cluster_parallel_class(cluster);
    cluster.division_class = ones(size(cluster.parallel_class,1),1) * 4;
    cluster.stich_list = generate_stich_list(cluster);
end

function mesh = generate_block_cluster_mesh(cluster, boundary_list, edge_control_option)
    % the edge_control_option is an optional input, in the form of N * 5
    % matrix, each row is of form [vertex1, vertex2, relative_first_offset, growth_rate, direction]
    mesh = struct();
    mesh.nodes = [];
    mesh.cell_node_pointer = [];
    mesh.cell_neighbour_pointer = [];
    mesh.m = 2;
    mesh.n = 2;
    block_mesh_list = [mesh];

    

    block_edge_control = -ones(length(cluster.block_list),6);
    if nargin > 2 && ~isempty(edge_control_option)
        [N,~] = size(edge_control_option);
            
        for setting_loop = 1:N

            control = edge_control_option(setting_loop, :);
            edge = control(:,1:2);
            options = control(:,3:5);
            parallel_class_index = 1;

            for parallel_class_search = 1:size(cluster.parallel_class,1)
                class = cluster.parallel_class{parallel_class_search};
                if ismember(edge, class, "rows")
                    parallel_class_index = parallel_class_search;
                    break;
                elseif ismember(flip(edge), class, "rows")
                    options(3) = -options(3);
                    parallel_class_index = parallel_class_search;
                    break;
                end
            end
            
            


            









            for block_edge_search = 1: length(cluster.block_list)
                % this loop is to fill the block_edge_control, or say, fill
                % it with a direct control vector for the block_mesh
                % function






                if ismember([cluster.block_list(block_edge_search).A, cluster.block_list(block_edge_search).B] , cluster.parallel_class{parallel_class_index}, 'rows')
                    %options(3) = -options(3);
                    block_edge_control(block_edge_search,:) = [options, [-1,-1,-1]];
                    %options(3) = -options(3);
                    continue;
                end
                if ismember([cluster.block_list(block_edge_search).B, cluster.block_list(block_edge_search).C] , cluster.parallel_class{parallel_class_index}, 'rows')
                    %options(3) = -options(3);
                    block_edge_control(block_edge_search,:) = [[-1,-1,-1], options];
                    %options(3) = -options(3);
                    continue;
                end

                if ismember([cluster.block_list(block_edge_search).C, cluster.block_list(block_edge_search).D] , cluster.parallel_class{parallel_class_index}, 'rows')
                    options(3) = -options(3);
                    block_edge_control(block_edge_search,:) = [options, [-1,-1,-1]];
                    options(3) = -options(3);
                    continue;
                end
                if ismember([cluster.block_list(block_edge_search).D, cluster.block_list(block_edge_search).A] , cluster.parallel_class{parallel_class_index}, 'rows')
                    options(3) = -options(3);
                    block_edge_control(block_edge_search,:) = [[-1,-1,-1], options];
                    options(3) = -options(3);
                    continue;
                end

                







                if ismember([cluster.block_list(block_edge_search).B, cluster.block_list(block_edge_search).A] , cluster.parallel_class{parallel_class_index}, 'rows')
                    options(3) = -options(3);
                    block_edge_control(block_edge_search,:) = [options, [-1,-1,-1]];
                    options(3) = -options(3);
                    continue;
                end
                if ismember([cluster.block_list(block_edge_search).C, cluster.block_list(block_edge_search).B] , cluster.parallel_class{parallel_class_index}, 'rows')
                    options(3) = -options(3);
                    block_edge_control(block_edge_search,:) = [[-1,-1,-1], options];
                    options(3) = -options(3);
                    continue;
                end

                if ismember([cluster.block_list(block_edge_search).D, cluster.block_list(block_edge_search).C] , cluster.parallel_class{parallel_class_index}, 'rows')
                    %options(3) = -options(3);
                    block_edge_control(block_edge_search,:) = [options, [-1,-1,-1]];
                    %options(3) = -options(3);
                    continue;
                end
                if ismember([cluster.block_list(block_edge_search).A, cluster.block_list(block_edge_search).D] , cluster.parallel_class{parallel_class_index}, 'rows')
                    %options(3) = -options(3);
                    block_edge_control(block_edge_search,:) = [[-1,-1,-1], options];
                    %options(3) = -options(3);
                    continue;
                end
            end
        end


        
    end





    for i = 1:length(cluster.block_list)
        nx = 2;
        ny = 2;
        for j = 1:length(cluster.division_class)
            if ismember([cluster.block_list(i).A, cluster.block_list(i).B] , cluster.parallel_class{j}, 'rows')
                nx = cluster.division_class(j);
            end
            if ismember([cluster.block_list(i).B, cluster.block_list(i).A] , cluster.parallel_class{j}, 'rows')
                nx = cluster.division_class(j);
            end

            if ismember([cluster.block_list(i).B, cluster.block_list(i).C] , cluster.parallel_class{j}, 'rows')
                ny = cluster.division_class(j);
            end
            if ismember([cluster.block_list(i).C, cluster.block_list(i).B] , cluster.parallel_class{j}, 'rows')
                ny = cluster.division_class(j);
            end
        end








        if (isempty(cluster.block_list(i).C1))
            C1 = segment(cluster.vertice_list(cluster.block_list(i).A, :), cluster.vertice_list(cluster.block_list(i).B, :));
        else
            C1 = cluster.block_list(i).C1;
        end

        if (isempty(cluster.block_list(i).C2))
            C2 = segment(cluster.vertice_list(cluster.block_list(i).B, :), cluster.vertice_list(cluster.block_list(i).C, :));
        else
            C2 = cluster.block_list(i).C2;
        end
        

        if (isempty(cluster.block_list(i).C3))
            C3 = segment(cluster.vertice_list(cluster.block_list(i).C, :), cluster.vertice_list(cluster.block_list(i).D, :));
        else
            C3 = cluster.block_list(i).C3;
        end

        if (isempty(cluster.block_list(i).C4))
            C4 = segment(cluster.vertice_list(cluster.block_list(i).D, :), cluster.vertice_list(cluster.block_list(i).A, :));
        else
            C4 = cluster.block_list(i).C4;
        end




        if nargin <= 2 || (block_edge_control(i, 1) == -1 && block_edge_control(i,4) == -1)
            block_mesh_list(i) = block_mesh(C1,C2,C3,C4,nx,ny);
        else
            if block_edge_control(i, 1) ~= -1 && block_edge_control(i, 4) ~= -1
                block_mesh_list(i) = block_mesh(C1,C2,C3,C4,nx,ny, block_edge_control(i, 1:3),block_edge_control(i, 4:6));
            elseif block_edge_control(i, 1) ~= -1
                block_mesh_list(i) = block_mesh(C1,C2,C3,C4,nx,ny, block_edge_control(i, 1:3),[]);
            elseif block_edge_control(i, 4) ~= -1
                block_mesh_list(i) = block_mesh(C1,C2,C3,C4,nx,ny, [],block_edge_control(i, 4:6));
            end
        end
    end
    
    mesh = stich_block_mesh(block_mesh_list, cluster.stich_list, boundary_list);
    

end

function plot_block_cluster(cluster)
    % PLOT_BLOCK_CLUSTER 在当前坐标轴中绘制 block_cluster 的顶点和块
    figure;
    vertice_list = cluster.vertice_list;
    block_list   = cluster.block_list;
    n_blocks     = length(block_list);
    n_vertices   = size(vertice_list, 1);

    edge_color = [0.8 0.4 0.8];
    line_width = 0.8;

    held = ishold;
    hold on; axis equal;

    for b = 1:n_blocks
        blk = block_list(b);
        verts_idx = [blk.A, blk.B, blk.C, blk.D];
        verts_xy  = vertice_list(verts_idx, :);

        % 闭合四边形 A-B-C-D-A
        poly = [verts_xy; verts_xy(1, :)];
        plot(poly(:,1), poly(:,2), '-', ...
             'Color', edge_color, 'LineWidth', line_width);

        % 块号标注
        cx = mean(verts_xy(:,1));
        cy = mean(verts_xy(:,2));
        text(cx, cy, sprintf('B%d', b), ...
             'FontSize', 8, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center', ...
             'Color', [0.8 0.4 0.8]);
    end

    % 顶点
    scatter(vertice_list(:,1), vertice_list(:,2), 50, ...
            'filled', 'MarkerFaceColor', [0.85 0.2 0.2], ...
            'MarkerEdgeColor', 'k');

    % 顶点编号
    xr = max(vertice_list(:,1)) - min(vertice_list(:,1));
    yr = max(vertice_list(:,2)) - min(vertice_list(:,2));
    off = 0.002 * max(xr, yr);
    if off == 0, off = 0.02; end
    for v = 1:n_vertices
        text(vertice_list(v,1) + off, vertice_list(v,2) + off, ...
             sprintf('V%d', v), ...
             'FontSize', 8, 'FontWeight', 'bold', ...
             'Color', [0.6 0.1 0.1]);
    end

    if ~held
        hold off;
    end
end

function check = plot_mesh(mesh)
    % input mesh:
    %
    %   mesh.nodes                  : number_of_cells x 2 table, node
    %                                   position
    %
    %   mesh.cell_node_pointer      : number_of_cells x 4 table, refer to
    %                                   node index
    %
    %   mesh.cell_neighbour_pointer : number_of_cells x 4 table, refer to
    %                                   the index of neighbourhood cells, in same order as the node
    %
    %   mesh.boundary_faces         : K×3，for each row [node_a, node_b, owner_cell_global]
    %                                 （golobal node index + owner cell index）
    %   mesh.boundary_index         : K×1，owner boundry group's boundary_index for each boundry faces 
    %
    %   mesh.boundary_groups        : cell list，sorted or classification by boundary_index for output into OpenFOAM afterwards
    %                                   notice that the index in here should be
    %                                   from mesh.boundary_index, the rest part of
    %                                   this cell are empty(undefined). In those
    %                                   defined part, it is a table of face
    %                                   information, in each row it represent a
    %                                   face, as 
    %                                   [node_a, node_b, owner_cell_global], 
    %                                   or say, it is a regroup
    %                                   of mesh.boundary_faces
    %
    %
    %
    %
    %
    %
    %
    % output:
    %   
    %   check if the boundary faces are assigned, and if the neighbourhood list is symmetry.   
    %
    
    
    check = true;

    % -------- 6. Validate results --------
    fprintf('\n========== Mesh statistics ==========\n');
    fprintf('Total nodes: %d\n', size(mesh.nodes, 1));
    fprintf('Total cells: %d\n', size(mesh.cell_node_pointer, 1));
    fprintf('Boundary faces (from boundary_faces): %d\n', size(mesh.boundary_faces, 1));

    % Check cell_node_pointer range
    cnp = mesh.cell_node_pointer;
    if all(cnp(:) >= 1) && all(cnp(:) <= size(mesh.nodes, 1))
        fprintf('[PASS] cell_node_pointer indices are within valid range\n');
    else
        fprintf('[FAIL] cell_node_pointer contains invalid indices\n');
        check = false;
    end

    % Check neighbour symmetry
    cnb = mesh.cell_neighbour_pointer;
    sym_ok = true;
    for c = 1:size(cnb, 1)
        for col = 1:4
            nb = cnb(c, col);
            if nb > 0
                if ~ismember(c, cnb(nb,:))
                    sym_ok = false;
                end
            elseif nb == -1
                % Boundary face, skip
            else
                sym_ok = false;
            end
        end
    end

    assign_ok = true;
    for c = 1:size(cnb, 1)
        for col = 1:4
            nb = cnb(c, col);
            if nb > 0
                % Internal edge, skip
            elseif nb == -1
                cell_node = mesh.cell_node_pointer(c, :);
                if col <= 3
                    n1 = col;
                    n2 = col + 1;
                else
                    n1 = 4;
                    n2 = 1;
                end
                n1 = cell_node(n1);
                n2 = cell_node(n2);
    
                
                assign_ok = assign_ok && (ismember([n1,n2],mesh.boundary_faces(:,[1,2]), "rows") || ismember([n2,n1],mesh.boundary_faces(:,[1,2]), "rows"));
            else
                assign_ok = false;
            end
        end
    end

    if sym_ok
        fprintf('[PASS] Neighbour pointers are mutually symmetric\n');
    else
        fprintf('[FAIL] Neighbour pointers are asymmetric or contain invalid values\n');
        check = false;
    end

    if assign_ok
        fprintf('[PASS] All boundary faces are assigned\n');
    else
        fprintf('[FAIL] Some boundary faces are not assigned\n');
        check = false;
    end

    % -------- 7. Visualization --------
    figure('Name', 'Block Cluster Mesh', 'Position', [100, 100, 800, 600]);
    hold on; axis equal;
    title(sprintf('Cluster mesh: %d nodes, %d cells', size(mesh.nodes,1), size(cnp,1)));

    % Draw cell edges (all)
    for c = 1:size(cnp, 1)
        verts = mesh.nodes(cnp(c, :), :);
        verts = [verts; verts(1,:)];
        plot(verts(:,1), verts(:,2), 'b-', 'LineWidth', 0.5);
    end

    % Draw nodes
    plot(mesh.nodes(:,1), mesh.nodes(:,2), 'r.', 'MarkerSize', 6);

    % -------- NEW: Highlight all one‑sided edges (boundary edges) --------
    % Build a unique set of boundary edges from cell_neighbour_pointer
    edge_set = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    for c = 1:size(cnp, 1)
        for e = 1:4
            if cnb(c, e) < 0   % boundary edge
                n1 = cnp(c, e);
                n2 = cnp(c, mod(e,4) + 1);   % next node in the quad
                if n1 > n2
                    [n1, n2] = deal(n2, n1);
                end
                key = sprintf('%d,%d', n1, n2);
                edge_set(key) = true;
            end
        end
    end
    all_boundary_edges = edge_set.keys();
    n_boundary_edges = numel(all_boundary_edges);
    fprintf('Total one‑sided edges (including boundary_faces): %d\n', n_boundary_edges);

    % Identify which of these are already in mesh.boundary_faces
    bf_set = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    if ~isempty(mesh.boundary_faces)
        for i = 1:size(mesh.boundary_faces, 1)
            n1 = mesh.boundary_faces(i,1);
            n2 = mesh.boundary_faces(i,2);
            if n1 > n2, [n1,n2] = deal(n2,n1); end
            bf_set(sprintf('%d,%d',n1,n2)) = true;
        end
    end

    % Plot the extra boundary edges (those not in boundary_faces) in dashed black
    extra_edges = {};
    for k = 1:numel(all_boundary_edges)
        key = all_boundary_edges{k};
        if ~isKey(bf_set, key)
            extra_edges{end+1} = key; %#ok<AGROW>
        end
    end

    if ~isempty(extra_edges)
        fprintf('Additional one‑sided edges NOT in boundary_faces: %d\n', numel(extra_edges));
        for k = 1:numel(extra_edges)
            parts = strsplit(extra_edges{k}, ',');
            n1 = str2double(parts{1});
            n2 = str2double(parts{2});
            pts = mesh.nodes([n1, n2], :);
            plot(pts(:,1), pts(:,2), '--k', 'LineWidth', 2.5);
        end
    else
        fprintf('All one‑sided edges are already present in boundary_faces.\n');
    end

    % Highlight boundary_faces with colors (existing code)
    if ~isempty(mesh.boundary_faces)
        colors = lines(max(mesh.boundary_index));
        for i = 1:size(mesh.boundary_faces, 1)
            face = mesh.boundary_faces(i, :);
            n1 = face(1); n2 = face(2);
            pts = mesh.nodes([n1, n2], :);
            bidx = mesh.boundary_index(i);
            plot(pts(:,1), pts(:,2), 'Color', colors(bidx, :), 'LineWidth', 0.6);
        end
    end

    xlabel('x'); ylabel('y');
    hold off;

    % -------- 8. Print boundary group information --------
    fprintf('\n========== Boundary groups ==========\n');
    for i = 1:length(mesh.boundary_groups.index)
        fprintf('Boundary group (ID=%d): %d faces\n', ...
            mesh.boundary_groups.index(i), ...
            size(mesh.boundary_groups.faces{i}, 1));
    end
end


%===================stator design tools===================

function [leading_edge, pressure_surface, trailing_edge, suction_surface, inlet_blade_cut, outlet_blade_cut] = stator_blade(alpha, beta, gamma, chord_length, leading_factor, trailing_factor, relative_maxium_thickness, relative_t_maxium_thickness, relative_leading_radius, relative_trailing_radius)
    %       |                         *
    %  α  _ *                         |
    % _ _/_ |__________o              *
    %       *                o        |
    %       |  *              \       *
    %       *   / *            \      |
    %       |—— γ    *          \  β  *
    %       *           *        \  __|
    %       |              *      \/  *
    %       *                 *    \  |
    %       |                    *  \ *
    %       *                       *\|
    %       |                         *
    %       *                         |
    %       |                         *
    %       *                         |
    %
    %
    %
    %   all the angle units are in radians, not degrees.
    %   some restrictions are needed:
    %   0 <= α <= π / 2
    %   0 < γ < π / 2
    %   -π / 8 <= β <= π / 2
    %




    A = [0,0];
    B = leading_factor * [sin(alpha), -cos(alpha)];
    D = [sin(gamma), -cos(gamma)];
    C = D + trailing_factor * [-sin(beta), cos(beta)];
    P = [A;B;C;D];
    pressure_surface_t = get_lower_surface_func(P, relative_maxium_thickness, relative_t_maxium_thickness, relative_leading_radius, relative_trailing_radius);
    suction_surface_t = get_upper_surface_func(P, relative_maxium_thickness, relative_t_maxium_thickness, relative_leading_radius, relative_trailing_radius);
    leading_edge_t = get_leading_edge_func(P, relative_maxium_thickness, relative_t_maxium_thickness, relative_leading_radius, relative_trailing_radius);
    trailing_edge_t = get_trailing_edge_func(P, relative_maxium_thickness, relative_t_maxium_thickness, relative_leading_radius, relative_trailing_radius);
    
    pressure_surface = @(t) chord_length * pressure_surface_t(t);
    suction_surface = @(t) chord_length * suction_surface_t(t);
    leading_edge = @(t) chord_length * leading_edge_t(t);
    trailing_edge = @(t) chord_length * trailing_edge_t(t);



    % 获取前缘上下端点
    le_upper = compute_upper_point(1, P, relative_t_maxium_thickness, relative_maxium_thickness, relative_leading_radius, relative_trailing_radius);
    le_lower = compute_lower_point(0, P, relative_t_maxium_thickness, relative_maxium_thickness, relative_leading_radius, relative_trailing_radius);
    
    % 计算前缘切线方向
    eps_val = 1e-6;
    le_upper_near = compute_upper_point(1 - eps_val, P, relative_t_maxium_thickness, relative_maxium_thickness, relative_leading_radius, relative_trailing_radius);
    tangent_le = (le_upper_near - le_upper) / eps_val;
    
    % 计算前缘圆弧的圆心和半径
    A = [2*(le_upper(1)-le_lower(1)), 2*(le_upper(2)-le_lower(2)); 
         tangent_le(1), tangent_le(2)];
    b = [le_upper(1)^2 + le_upper(2)^2 - le_lower(1)^2 - le_lower(2)^2;
         le_upper(1)*tangent_le(1) + le_upper(2)*tangent_le(2)];
    le_center = (A \ b)';
    le_radius = norm(le_upper - le_center);
    
    % 计算角度范围
    angle_upper = atan2(le_upper(2) - le_center(2), le_upper(1) - le_center(1));
    angle_lower = atan2(le_lower(2) - le_center(2), le_lower(1) - le_center(1));
    if angle_lower < angle_upper
        angle_lower = angle_lower + 2*pi;
    end
    
    inlet_blade_cut = (max(pi / 2, angle_upper) - angle_upper) / (angle_lower - angle_upper);
    
    



    % 获取后缘上下端点
    te_upper = compute_upper_point(0, P, relative_t_maxium_thickness, relative_maxium_thickness, relative_leading_radius, relative_trailing_radius);
    te_lower = compute_lower_point(1, P, relative_t_maxium_thickness, relative_maxium_thickness, relative_leading_radius, relative_trailing_radius);
    
    % 计算后缘切线方向
    eps_val = 1e-6;
    te_upper_near = compute_upper_point(eps_val, P, relative_t_maxium_thickness, relative_maxium_thickness, relative_leading_radius, relative_trailing_radius);
    tangent_te = (te_upper - te_upper_near) / eps_val;
    
    % 计算后缘圆弧的圆心和半径
    A = [2*(te_upper(1)-te_lower(1)), 2*(te_upper(2)-te_lower(2)); 
         tangent_te(1), tangent_te(2)];
    b = [te_upper(1)^2 + te_upper(2)^2 - te_lower(1)^2 - te_lower(2)^2;
         te_upper(1)*tangent_te(1) + te_upper(2)*tangent_te(2)];
    te_center = (A \ b)';
    te_radius = norm(te_upper - te_center);
    
    % 计算角度范围
    angle_upper = atan2(te_upper(2) - te_center(2), te_upper(1) - te_center(1));
    angle_lower = atan2(te_lower(2) - te_center(2), te_lower(1) - te_center(1));
    if angle_lower > angle_upper
        angle_lower = angle_lower - 2*pi;
    end
    
    outlet_blade_cut = (max( (-pi) / 2, angle_lower) - angle_lower) / (angle_upper - angle_lower);


end

function [suction_wall, pressure_wall, le_round, te_suction, te_pressure, inlet_period_1, inlet_period_2, outlet_period_1, outlet_period_2, inlet, outlet] = fluid_domain_boundry(leading_edge, pressure_surface, trailing_edge, suction_surface, inlet_blade_cut, outlet_blade_cut, pitch, inlet_extension, outlet_extension )
    %       |                         *
    %  α  _ *                         |
    % _ _/_ |__________o              *
    %       *                o        |
    %       |  *              \       *
    %       *   / *            \      |
    %       |—— γ    *          \  β  *
    %       *           *        \  __|
    %       |              *      \/  *
    %       *                 *    \  |
    %       |                    *  \ *
    %       *                       *\|
    %       |                         *
    %       *                         |
    %       |                         *
    %       *                         |
    %
    %
    %
    %   all the angle units are in radians, not degrees.
    %   some restrictions are needed:
    %   0 <= α <= π / 2
    %   0 < γ < π / 2
    %   -π / 8 <= β <= π / 2
    %


    
    inlet_cut = leading_edge(inlet_blade_cut);
    outlet_cut = trailing_edge(outlet_blade_cut);


    inlet_period_1 = @(t) (1 - t) * (inlet_cut + inlet_extension * [0,1]) + t * inlet_cut;
    outlet_period_1 =  @(t) t * (outlet_cut - outlet_extension * [0,1]) + (1 - t) * outlet_cut;
    inlet_period_2 = @(t) inlet_period_1(1 - t) + [pitch,0];
    outlet_period_2 = @(t) outlet_period_1(1 - t) + [pitch,0];
    inlet = @(t) (1 - t) * (inlet_cut + inlet_extension * [0,1] + [pitch,0]) + t * (inlet_cut + inlet_extension * [0,1]);
    outlet = @(t) t * (outlet_cut - outlet_extension * [0,1] + [pitch,0]) + (1 - t) * (outlet_cut - outlet_extension * [0,1]);

    le_round = @(t) leading_edge(inlet_blade_cut + (1 - t) * (1 - inlet_blade_cut)) + [pitch,0];



    % Length of the LE-arc piece that belongs to the suction wall
    le_length = integral(@(t) norm(inlet_blade_cut * ...
             (leading_edge(t*inlet_blade_cut + 1e-6) - leading_edge(t*inlet_blade_cut))/1e-6), ...
             0, 1, 'ArrayValued', true);

    % Length of the suction surface
    suction_surface_length = integral(@(t) norm( (suction_surface(min(t+1e-6,1)) - suction_surface(t))/1e-6 ), ...
             0, 1, 'ArrayValued', true);

    r = le_length / suction_surface_length;
    suction_surface_start_point = r/(1 + r);




    function point = suction_wall_function(t)
        
        if t <= suction_surface_start_point
            point = leading_edge(inlet_blade_cut * (1 - t / suction_surface_start_point));
        else
            point = suction_surface((1 - t) / (1 - suction_surface_start_point));
        end
    end
    


    if r == 0
        suction_wall = @(t) suction_surface(1 - t);
    else
        suction_wall = @suction_wall_function;
    end
    




    te_suction = @(t) trailing_edge(1 - t * (1 - outlet_blade_cut));
    te_pressure = @(t) trailing_edge((1 - t) * outlet_blade_cut) + [pitch,0];
    pressure_wall = @(t) pressure_surface(1 - t) + [pitch,0];
    

end

function mesh = generate_fluid_domain_mesh(alpha, beta, gamma, chord_length, leading_factor, trailing_factor, relative_maxium_thickness, relative_t_maxium_thickness, relative_leading_edge_radius, relative_trailing_edge_radius, pitch, inlet_extension, outlet_extension)
    % output mesh:
    %
    %   mesh.nodes                  : number_of_cells x 2 table, node
    %                                   position
    %
    %   mesh.cell_node_pointer      : number_of_cells x 4 table, refer to
    %                                   node index
    %
    %   mesh.cell_neighbour_pointer : number_of_cells x 4 table, refer to
    %                                   the index of neighbourhood cells, in same order as the node
    %
    %   mesh.boundary_faces         : K×3，for each row [node_a, node_b, owner_cell_global]
    %                                 （golobal node index + owner cell index）
    %   mesh.boundary_index         : K×1，owner boundry group's boundary_index for each boundry faces 
    %
    %   mesh.boundary_groups        : cell list，sorted or classification by boundary_index for output into OpenFOAM afterwards
    %                                   notice that the index in here should be
    %                                   from mesh.boundary_index, the rest part of
    %                                   this cell are empty(undefined). In those
    %                                   defined part, it is a table of face
    %                                   information, in each row it represent a
    %                                   face, as 
    %                                   [node_a, node_b, owner_cell_global], 
    %                                   or say, it is a regroup
    %                                   of mesh.boundary_faces
    %
    
    



    % mesh density settings
    principle_laterial_division = 61;
    inlet_division = 38;
    outlet_division = 151;
    pressure_wall_division = 76;
    suction_wall_pre_mid_part_1_division = 16;
    suction_wall_pre_mid_part_2_division = 16;
    secondary_laterial_division = 11;
    trailing_edge_division = 16;




    % boundary layer settings
    boundary_layer = true;
    absolute_boundary_layer_thickness = 0.000001;
    growth_rate = 1.3;
    boundry_layer_division = 29;








    
    pressure_boundry_layer_thickness = min(0.75 * relative_leading_edge_radius * chord_length, min(0.5 * cos(alpha) * pitch, chord_length * 0.1));
    suction_boundry_layer_thickness = min(0.75 * relative_leading_edge_radius * chord_length, min(0.5 * cos(alpha) * pitch, chord_length * 0.1));
    suction_wall_mid_node_position_1 = (2 * relative_leading_edge_radius * chord_length + sin(alpha) * (pressure_boundry_layer_thickness)) / chord_length;
    suction_wall_mid_node_position_2 = (2.5 * relative_leading_edge_radius * chord_length + sin(alpha) * (pressure_boundry_layer_thickness + suction_boundry_layer_thickness)) / chord_length;


    


    edge_control_option = [];
    if boundary_layer
        edge_control_option = [2 ,5 , absolute_boundary_layer_thickness / pressure_boundry_layer_thickness , growth_rate, 1; ...
                           20,28, absolute_boundary_layer_thickness / suction_boundry_layer_thickness  , growth_rate, 1];
    end




    


    




    % generate the blade
    [le, lower, te, upper, inlet_blade_cut, outlet_blade_cut] = ...
    stator_blade(alpha, beta, gamma, chord_length, leading_factor, trailing_factor, relative_maxium_thickness, relative_t_maxium_thickness, relative_leading_edge_radius, relative_trailing_edge_radius);

    % generate the fluid domain boundry
    [suction_wall, pressure_wall, le_round, te_suction, te_pressure, ...
    inlet_period_1, inlet_period_2, outlet_period_1, outlet_period_2, ...
    inlet, outlet] = fluid_domain_boundry(le, lower, te, upper, ...
                        inlet_blade_cut, outlet_blade_cut, ...
                        pitch, inlet_extension, outlet_extension);
    
    vertice_list = zeros(29,2);

    vertice_list(1,:) = inlet_period_1(0);
    vertice_list(2,:) = inlet_period_1(1);
    vertice_list(3,:) = outlet_period_1(0);
    vertice_list(4,:) = outlet_period_1(1);


    vertice_list(19,:) = inlet_period_2(1);
    vertice_list(20,:) = inlet_period_2(0);
    vertice_list(24,:) = outlet_period_2(1);
    vertice_list(26,:) = outlet_period_2(0);

    vertice_list(21,:) = le_round(0.5);
    vertice_list(22,:) = le_round(0);
    vertice_list(23,:) = te_pressure(1);

    vertice_list(11,:) = te_suction(0);


    

    vertice_list(27,:) = suction_wall(suction_wall_mid_node_position_1);
    vertice_list(6,:) = suction_wall(suction_wall_mid_node_position_2);
    
    c_rotation_mat = [0,1;-1,0];
    %ac_rotation_mat = [0,-1;1,0]

    vertice_list(7,:) = vertice_list(27,:) + suction_boundry_layer_thickness * curve_tangent_vector(suction_wall, suction_wall_mid_node_position_1) * c_rotation_mat;
    vertice_list(8,:) = vertice_list(6,:) + suction_boundry_layer_thickness * curve_tangent_vector(suction_wall, suction_wall_mid_node_position_2) * c_rotation_mat;
    
    
    vertice_list(5,:) = vertice_list(2,:) + suction_boundry_layer_thickness * [0,1];
    
    vertice_list(9,:) = vertice_list(11,:) + suction_boundry_layer_thickness * curve_tangent_vector(suction_wall, 1) * c_rotation_mat;
    vertice_list(10,:) = vertice_list(3,:) - suction_boundry_layer_thickness * [0,1];



    
    
    vertice_list(15,:) = vertice_list(21,:) + pressure_boundry_layer_thickness * curve_tangent_vector(le_round, 0.5) * c_rotation_mat;
    
    vertice_list(16,:) = vertice_list(22,:) + pressure_boundry_layer_thickness * curve_tangent_vector(pressure_wall, 1) * c_rotation_mat;
    vertice_list(17,:) = vertice_list(23,:) + pressure_boundry_layer_thickness * curve_tangent_vector(pressure_wall, 0) * c_rotation_mat;
    vertice_list(25,:) = vertice_list(24,:) - pressure_boundry_layer_thickness * [0,1];


    vertice_list(12,:) = 0.7 * vertice_list(4,:) + 0.3 * vertice_list(26,:);
    vertice_list(18,:) = 0.3 * vertice_list(4,:) + 0.7 * vertice_list(26,:);


    vertice_list(29,:) = le_round(0.75);
    vertice_list(14,:) = vertice_list(29,:) + pressure_boundry_layer_thickness * curve_tangent_vector(le_round, 0.75) * c_rotation_mat;
    vertice_list(28,:) = vertice_list(20,:) + pressure_boundry_layer_thickness * [0,1];
    vertice_list(13,:) = [vertice_list(14,1), vertice_list(19,2)];


    
    block_list(17) = struct('A', [], ...
                    'B', [], ...
                    'C', [], ...
                    'D', [], ...
                    'C1', function_handle.empty, ...
                    'C2', function_handle.empty, ...
                    'C3', function_handle.empty, ...
                    'C4', function_handle.empty);

    block_list(1).A = 1;
    block_list(1).B = 5;
    block_list(1).C = 14;
    block_list(1).D = 13;

    block_list(2).A = 13;
    block_list(2).B = 14;
    block_list(2).C = 28;
    block_list(2).D = 19;

    block_list(3).A = 5;
    block_list(3).B = 7;
    block_list(3).C = 15;
    block_list(3).D = 14;
    
    block_list(4).A = 7;
    block_list(4).B = 8;
    block_list(4).C = 16;
    block_list(4).D = 15;

    block_list(5).A = 8;
    block_list(5).B = 9;
    block_list(5).C = 17;
    block_list(5).D = 16;

    block_list(6).A = 10;
    block_list(6).B = 4;
    block_list(6).C = 12;
    block_list(6).D = 9;

    block_list(7).A = 9;
    block_list(7).B = 12;
    block_list(7).C = 18;
    block_list(7).D = 17;

    block_list(8).A = 17;
    block_list(8).B = 18;
    block_list(8).C = 26;
    block_list(8).D = 25;

    block_list(9).A = 5;
    block_list(9).B = 2;
    block_list(9).C = 27;
    block_list(9).D = 7;

    block_list(10).A = 27;
    block_list(10).B = 6;
    block_list(10).C = 8;
    block_list(10).D = 7;

    block_list(11).A = 6;
    block_list(11).B = 11;
    block_list(11).C = 9;
    block_list(11).D = 8;

    block_list(12).A = 3;
    block_list(12).B = 10;
    block_list(12).C = 9;
    block_list(12).D = 11;

    block_list(13).A = 14;
    block_list(13).B = 29;
    block_list(13).C = 20;
    block_list(13).D = 28;

    block_list(14).A = 14;
    block_list(14).B = 15;
    block_list(14).C = 21;
    block_list(14).D = 29;

    block_list(15).A = 15;
    block_list(15).B = 16;
    block_list(15).C = 22;
    block_list(15).D = 21;

    block_list(16).A = 16;
    block_list(16).B = 17;
    block_list(16).C = 23;
    block_list(16).D = 22;

    block_list(17).A = 23;
    block_list(17).B = 17;
    block_list(17).C = 25;
    block_list(17).D = 24;


    block_list(9).C2 = @(t) suction_wall(suction_wall_mid_node_position_1 * t);
    block_list(10).C1 = @(t) suction_wall((suction_wall_mid_node_position_2 - suction_wall_mid_node_position_1) * t + suction_wall_mid_node_position_1);
    block_list(11).C1 = @(t) suction_wall((1 - suction_wall_mid_node_position_2) * t + suction_wall_mid_node_position_2);
    block_list(12).C4 = @(t) te_suction(t);

    block_list(13).C2 = @(t) le_round((1 - 0.75) * t + 0.75);
    block_list(14).C3 = @(t) le_round((0.75 - 0.5) * t + 0.5);
    block_list(15).C3 = @(t) le_round((0.5 - 0) * t + 0);

    block_list(16).C3 = pressure_wall;
    block_list(17).C4 = te_pressure;



    
    seg = segment(vertice_list(block_list(9).D,:), vertice_list(block_list(9).A,:));
    b = (vertice_list(block_list(9).D,:) - vertice_list(block_list(9).C,:));
    base_point_pair = (vertice_list(block_list(9).B,:) - vertice_list(block_list(9).C,:));
    target_point_pair = (vertice_list(block_list(9).A,:) - vertice_list(block_list(9).D,:));
    dot_p = dot(base_point_pair, target_point_pair);
    cross = base_point_pair(1)*target_point_pair(2) - base_point_pair(2)*target_point_pair(1);
    norm_product = norm(base_point_pair) * norm(target_point_pair);
    scaling = norm(target_point_pair) / norm(base_point_pair);
    cos_theta = dot_p / norm_product;
    sin_theta = cross / norm_product;



    adjust_rotation_mat = [cos_theta, sin_theta; -sin_theta, cos_theta];
    cur = @(t) 0.4 * (scaling * (block_list(9).C2(1 - t) - block_list(9).C2(1)) * adjust_rotation_mat + b + block_list(9).C2(1)) + 0.6 * seg(t); 
    block_list(9).C4 = cur;
    block_list(3).C1 = @(t) cur(1 - t);





    seg = segment(vertice_list(block_list(12).B,:), vertice_list(block_list(12).C,:));
    b = (vertice_list(block_list(12).B,:) - vertice_list(block_list(12).A,:));
    base_point_pair = (vertice_list(block_list(12).D,:) - vertice_list(block_list(12).A,:));
    target_point_pair = (vertice_list(block_list(12).C,:) - vertice_list(block_list(12).B,:));
    dot_p = dot(base_point_pair, target_point_pair);
    cross = base_point_pair(1)*target_point_pair(2) - base_point_pair(2)*target_point_pair(1);
    norm_product = norm(base_point_pair) * norm(target_point_pair);
    scaling = norm(target_point_pair) / norm(base_point_pair);
    cos_theta = dot_p / norm_product;
    sin_theta = cross / norm_product;



    adjust_rotation_mat = [cos_theta, sin_theta; -sin_theta, cos_theta];
    cur = @(t) 0.5 * (scaling * (te_suction(1 - t) - te_suction(1)) * adjust_rotation_mat + b + te_suction(1)) + 0.5 * seg(t); 
    block_list(12).C2 = cur;
    block_list(6).C4 = @(t) cur(1-t);


    seg = segment(vertice_list(block_list(13).D,:), vertice_list(block_list(13).A,:));
    b = (vertice_list(block_list(13).D,:) - vertice_list(block_list(13).C,:));
    base_point_pair = (vertice_list(block_list(13).B,:) - vertice_list(block_list(13).C,:));
    target_point_pair = (vertice_list(block_list(13).A,:) - vertice_list(block_list(13).D,:));
    dot_p = dot(base_point_pair, target_point_pair);
    cross = base_point_pair(1)*target_point_pair(2) - base_point_pair(2)*target_point_pair(1);
    norm_product = norm(base_point_pair) * norm(target_point_pair);
    scaling = norm(target_point_pair) / norm(base_point_pair);
    cos_theta = dot_p / norm_product;
    sin_theta = cross / norm_product;



    adjust_rotation_mat = [cos_theta, sin_theta; -sin_theta, cos_theta];
    cur = @(t) 0.5 * (scaling * (block_list(13).C2(1 - t) - block_list(13).C2(1)) * adjust_rotation_mat + b + block_list(13).C2(1)) + 0.5 * seg(t); 
    block_list(13).C4 = cur;
    block_list(2).C2 = @(t) cur(1 - t);



    
    seg = segment(vertice_list(block_list(10).C,:), vertice_list(block_list(10).D,:));
    b = (vertice_list(block_list(10).C,:) - vertice_list(block_list(10).B,:));
    base_point_pair = (vertice_list(block_list(10).A,:) - vertice_list(block_list(10).B,:));
    target_point_pair = (vertice_list(block_list(10).D,:) - vertice_list(block_list(10).C,:));
    dot_p = dot(base_point_pair, target_point_pair);
    cross = base_point_pair(1)*target_point_pair(2) - base_point_pair(2)*target_point_pair(1);
    norm_product = norm(base_point_pair) * norm(target_point_pair);
    scaling = norm(target_point_pair) / norm(base_point_pair);
    cos_theta = dot_p / norm_product;
    sin_theta = cross / norm_product;



    adjust_rotation_mat = [cos_theta, sin_theta; -sin_theta, cos_theta];
    cur = @(t) 0.5 * (scaling * (block_list(10).C1(1 - t) - block_list(10).C1(1)) * adjust_rotation_mat + b + block_list(10).C1(1)) + 0.5 * seg(t); 
    block_list(10).C3 = cur;
    block_list(4).C1 = @(t) cur(1 - t);


    % =========================================================================
    % the fallowing part is for adding a curvature to boundary layer at
    % suction side, the third, or main part, this is for avoid result a
    % mesh with negative volum cells
    %
    %
    %
    % =========================================================================
    seg = segment(vertice_list(block_list(11).C,:), vertice_list(block_list(11).D,:));
    b = (vertice_list(block_list(11).C,:) - vertice_list(block_list(11).B,:));
    base_point_pair = (vertice_list(block_list(11).A,:) - vertice_list(block_list(11).B,:));
    target_point_pair = (vertice_list(block_list(11).D,:) - vertice_list(block_list(11).C,:));
    dot_p = dot(base_point_pair, target_point_pair);
    cross = base_point_pair(1)*target_point_pair(2) - base_point_pair(2)*target_point_pair(1);
    norm_product = norm(base_point_pair) * norm(target_point_pair);
    scaling = norm(target_point_pair) / norm(base_point_pair);
    cos_theta = dot_p / norm_product;
    sin_theta = cross / norm_product;



    adjust_rotation_mat = [cos_theta, sin_theta; -sin_theta, cos_theta];
    % adjust the factor (it is a linear interpolation of a curve and a
    % segment) in order to avoid negative volum or highly skewed cells,
    % mostly it is to increase it, recommened range is 0.3-0.8
    cur = @(t) 0.95 * (scaling * (block_list(11).C1(1 - t) - block_list(11).C1(1)) * adjust_rotation_mat + b + block_list(11).C1(1)) + 0.05 * seg(t); 
    block_list(11).C3 = cur;
    block_list(5).C1 = @(t) cur(1 - t);




    










    % =========================================================================
    % the fallowing part is for adding a curvature to boundary layer at
    % the pressure side leading edge, the third, or main part, this is for avoid result a
    % mesh with negative volum cells
    %
    %
    %
    % =========================================================================

    seg = segment(vertice_list(block_list(14).A,:), vertice_list(block_list(14).B,:));
    b = (vertice_list(block_list(14).A,:) - vertice_list(block_list(14).D,:));
    base_point_pair = (vertice_list(block_list(14).C,:) - vertice_list(block_list(14).D,:));
    target_point_pair = (vertice_list(block_list(14).B,:) - vertice_list(block_list(14).A,:));
    dot_p = dot(base_point_pair, target_point_pair);
    cross = base_point_pair(1)*target_point_pair(2) - base_point_pair(2)*target_point_pair(1);
    norm_product = norm(base_point_pair) * norm(target_point_pair);
    scaling = norm(target_point_pair) / norm(base_point_pair);
    cos_theta = dot_p / norm_product;
    sin_theta = cross / norm_product;



    adjust_rotation_mat = [cos_theta, sin_theta; -sin_theta, cos_theta];
    cur = @(t) 0.5 * (scaling * (block_list(14).C3(1 - t) - block_list(14).C3(1)) * adjust_rotation_mat + b + block_list(14).C3(1)) + 0.5 * seg(t); 
    block_list(14).C1 = cur;
    block_list(3).C3 = @(t) cur(1 - t);


 




    seg = segment(vertice_list(block_list(15).A,:), vertice_list(block_list(15).B,:));
    b = (vertice_list(block_list(15).A,:) - vertice_list(block_list(15).D,:));
    base_point_pair = (vertice_list(block_list(15).C,:) - vertice_list(block_list(15).D,:));
    target_point_pair = (vertice_list(block_list(15).B,:) - vertice_list(block_list(15).A,:));
    dot_p = dot(base_point_pair, target_point_pair);
    cross = base_point_pair(1)*target_point_pair(2) - base_point_pair(2)*target_point_pair(1);
    norm_product = norm(base_point_pair) * norm(target_point_pair);
    scaling = norm(target_point_pair) / norm(base_point_pair);
    cos_theta = dot_p / norm_product;
    sin_theta = cross / norm_product;



    adjust_rotation_mat = [cos_theta, sin_theta; -sin_theta, cos_theta];
    cur = @(t) 0.5 * (scaling * (block_list(15).C3(1 - t) - block_list(15).C3(1)) * adjust_rotation_mat + b + block_list(15).C3(1)) + 0.5 * seg(t); 
    block_list(15).C1 = cur;
    block_list(4).C3 = @(t) cur(1 - t);
















    % =========================================================================
    % the fallowing part is for adding a curvature to boundary layer at
    % pressure side, main part, this is for avoid result a
    % mesh with negative volum cells
    %
    %
    %
    % =========================================================================
    seg = segment(vertice_list(block_list(16).A,:), vertice_list(block_list(16).B,:));
    b = (vertice_list(block_list(16).A,:) - vertice_list(block_list(16).D,:));
    base_point_pair = (vertice_list(block_list(16).C,:) - vertice_list(block_list(16).D,:));
    target_point_pair = (vertice_list(block_list(16).B,:) - vertice_list(block_list(16).A,:));
    dot_p = dot(base_point_pair, target_point_pair);
    cross = base_point_pair(1)*target_point_pair(2) - base_point_pair(2)*target_point_pair(1);
    norm_product = norm(base_point_pair) * norm(target_point_pair);
    scaling = norm(target_point_pair) / norm(base_point_pair);
    cos_theta = dot_p / norm_product;
    sin_theta = cross / norm_product;



    adjust_rotation_mat = [cos_theta, sin_theta; -sin_theta, cos_theta];

    % adjust the factor (it is a linear interpolation of a curve and a
    % segment) in order to avoid negative volum or highly skewed cells,
    % mostly it is to increase it, recommened range is 0.5-0.8
    cur = @(t) 0.8 * (scaling * (block_list(16).C3(1 - t) - block_list(16).C3(1)) * adjust_rotation_mat + b + block_list(16).C3(1)) + 0.2 * seg(t); 
    block_list(16).C1 = cur;
    block_list(5).C3 = @(t) cur(1 - t);


    seg = segment(vertice_list(block_list(17).B,:), vertice_list(block_list(17).C,:));
    b = (vertice_list(block_list(17).A,:) - vertice_list(block_list(17).D,:));
    base_point_pair = (vertice_list(block_list(17).D,:) - vertice_list(block_list(17).A,:));
    target_point_pair = (vertice_list(block_list(17).C,:) - vertice_list(block_list(17).B,:));
    dot_p = dot(base_point_pair, target_point_pair);
    cross = base_point_pair(1)*target_point_pair(2) - base_point_pair(2)*target_point_pair(1);
    norm_product = norm(base_point_pair) * norm(target_point_pair);
    scaling = norm(target_point_pair) / norm(base_point_pair);
    cos_theta = dot_p / norm_product;
    sin_theta = cross / norm_product;



    adjust_rotation_mat = [cos_theta, sin_theta; -sin_theta, cos_theta];
    cur = @(t) 0.5 * (scaling * (block_list(17).C4(1 - t) - block_list(17).C4(1)) * adjust_rotation_mat + b + block_list(17).C4(1)) + 0.5 * seg(t); 
    block_list(17).C2 = cur;
    block_list(8).C4 = @(t) cur(1 - t);
    





    seg = segment(vertice_list(block_list(1).B,:), vertice_list(block_list(1).C,:));
    cur = @(t) seg(t) + min((inlet_extension - min(pressure_boundry_layer_thickness, suction_boundry_layer_thickness)) / 2,  min(pitch / 4, 0.4 * chord_length * leading_factor * sin(alpha))) * [0,1 - 4 * (t - 0.5)^2];
    block_list(1).C2 = @(t) cur(t);
    block_list(3).C4 = @(t) cur(1 - t);








    block_cluster = block_cluster_initialize(block_list, vertice_list);
            
    boundry_list = [
        1 , 1, 1;   % block1   C1 -> boundry 1  (inlet_period_1)
        9 , 1, 1;   % block9   C1 -> boundry 1  (inlet_period_1)
        9 , 2, 2;   % block9   C2 -> boundry 2  (suction_wall)
        10, 1, 2;   % block10  C1 -> boundry 2  (suction_wall)
        11, 1, 2;   % block11  C1 -> boundry 2  (suction_wall)
        12, 4, 3;   % block12  C4 -> boundry 3  (te_suction)
        12, 1, 4;   % block12  C1 -> boundry 4  (outlet_period_1)
        6 , 1, 4;   % block6   C1 -> boundry 4  (outlet_period_1)
        1 , 4, 5;   % block1   C4 -> boundry 5  (inlet)
        2 , 4, 5;   % block2   C4 -> boundry 5  (inlet)
        6 , 2, 6;   % block6   C1 -> boundry 6  (outlet)
        7 , 2, 6;   % block7   C2 -> boundry 6  (outlet)
        8 , 2, 6;   % block8   C2 -> boundry 6  (outlet)
        2 , 3, 7;   % block2   C3 -> boundry 7  (inlet_period_2)
        13, 3, 7;   % block13  C3 -> boundry 7  (inlet_period_2)
        17, 3, 8;   % block17  C3 -> boundry 8  (outlet_period_1)
        8 , 3, 8;   % block8   C3 -> boundry 8  (outlet_period_1)
        13, 2, 9;   % block13  C2 -> boundry 9  (le_round)
        14, 3, 9;   % block14  C3 -> boundry 9  (le_round)
        15, 3, 9;   % block15  C3 -> boundry 9  (le_round)
        16, 3, 10;  % block16  C3 -> boundry 10 (pressure_wall)
        17, 4, 11;  % block17  C4 -> boundry 11 (te_pressure)
    ];

    
    block_cluster.division_class = 2 * block_cluster.division_class;
    


    for k = 1:length(block_cluster.parallel_class)
        edges = block_cluster.parallel_class{k};
        
        for j = 1:size(edges, 1)
            %fprintf('(%d,%d) ', edges(j,1), edges(j,2));
        end

        if ismember([1,13], edges,"rows") || ismember([13,1], edges,"rows")
            block_cluster.division_class(k) = principle_laterial_division;
        end
        if ismember([2,5], edges,"rows") || ismember([5,2], edges,"rows")
            block_cluster.division_class(k) = boundry_layer_division;
        end
        if ismember([16, 22], edges,"rows") || ismember([22,16], edges,"rows")
            block_cluster.division_class(k) = boundry_layer_division;
        end
        if ismember([1,5], edges,"rows") || ismember([5,1], edges,"rows")
           block_cluster.division_class(k) = inlet_division;
        end
        if ismember([12,9], edges,"rows") || ismember([9,12], edges,"rows")
           block_cluster.division_class(k) = outlet_division;
        end
        if ismember([8,9], edges,"rows") || ismember([9,8], edges,"rows")
           block_cluster.division_class(k) = pressure_wall_division;
        end
        if ismember([5,7], edges,"rows") || ismember([7,5], edges,"rows")
           block_cluster.division_class(k) = suction_wall_pre_mid_part_1_division;
        end
        if ismember([8,7], edges,"rows") || ismember([7,8], edges,"rows")
           block_cluster.division_class(k) = suction_wall_pre_mid_part_2_division;
        end
        if ismember([13,19], edges,"rows") || ismember([19,13], edges,"rows")
           block_cluster.division_class(k) = secondary_laterial_division;
        end
        if ismember([4,12], edges,"rows") || ismember([12,4], edges,"rows")
           %block_cluster.division_class(k) = secondary_laterial_division;
        end
        if ismember([17,25], edges,"rows") || ismember([25,17], edges,"rows")
           block_cluster.division_class(k) = secondary_laterial_division;
        end
        if ismember([10,9], edges,"rows") || ismember([9,11], edges,"rows")
            block_cluster.division_class(k) = trailing_edge_division;
        end
        if ismember([17,25], edges,"rows") || ismember([25,17], edges,"rows")
            block_cluster.division_class(k) = trailing_edge_division;
        end
        %fprintf('\n');
    end

    
    %plot_block_cluster(block_cluster);
    mesh = generate_block_cluster_mesh(block_cluster,boundry_list, edge_control_option);
    %plot_mesh(mesh);
    
end

function export_mesh_to_starccm(mesh, output_file, boundary_names)
%EXPORT_MESH_TO_STARCCM Export 2D quad mesh to ASCII Fluent .msh format.
%
%   export_mesh_to_starccm(mesh, output_file, boundary_names)
%
% This version rebuilds face connectivity from cell_node_pointer instead of
% trusting cell_neighbour_pointer. Therefore, if an edge is shared by two
% cells, it is exported as an interior face even if one cell wrongly marks
% the edge as one-sided.
%
% Required mesh fields:
%   mesh.nodes                  : Nnodes x 2
%   mesh.cell_node_pointer      : Ncells x 4
%   mesh.cell_neighbour_pointer : Ncells x 4
%
% Optional mesh fields:
%   mesh.boundary_faces         : K x 3, [node_a, node_b, owner_cell]
%   mesh.boundary_index         : K x 1
%
% boundary_names:
%   cell array, boundary_names{i} gives name of boundary index i.

    if nargin < 3
        boundary_names = {};
    end

    % ================================================================
    % 1. Basic validation
    % ================================================================
    required_fields = {'nodes', 'cell_node_pointer', 'cell_neighbour_pointer'};
    for k = 1:numel(required_fields)
        if ~isfield(mesh, required_fields{k})
            error('Mesh is missing required field: %s', required_fields{k});
        end
    end

    nodes_original = mesh.nodes;
    cnp_original   = mesh.cell_node_pointer;
    cnb            = mesh.cell_neighbour_pointer;

    if size(nodes_original, 2) ~= 2
        error('mesh.nodes must be N x 2.');
    end

    if size(cnp_original, 2) ~= 4
        error('mesh.cell_node_pointer must be Ncells x 4.');
    end

    if size(cnb, 2) ~= 4
        error('mesh.cell_neighbour_pointer must be Ncells x 4.');
    end

    n_cells = size(cnp_original, 1);

    if size(cnb, 1) ~= n_cells
        error('mesh.cell_neighbour_pointer must have same number of rows as cell_node_pointer.');
    end

    if any(cnp_original(:) < 1) || any(cnp_original(:) > size(nodes_original, 1))
        error('mesh.cell_node_pointer contains node indices outside valid range.');
    end

    % ================================================================
    % 2. Merge duplicate or nearly duplicate nodes
    % ================================================================
    %
    % This helps when an internal interface is geometrically coincident but
    % the two sides still use different node IDs.
    %
    % If you do not want tolerance-based node merging, set merge_tol = 0.
    %
    bbox_size = max(max(nodes_original, [], 1) - min(nodes_original, [], 1));
    if bbox_size == 0
        merge_tol = 0;
    else
        merge_tol = 1.0e-10 * bbox_size;
    end

    [nodes, cnp, old_to_new, n_merged] = merge_duplicate_nodes( ...
        nodes_original, cnp_original, merge_tol);

    n_nodes = size(nodes, 1);

    if n_merged > 0
        %fprintf('Node merge: %d duplicate / near-duplicate nodes merged using tolerance %.3e.\n', n_merged, merge_tol);
    end

    % Check for degenerate cells after node merging.
    for c = 1:n_cells
        if numel(unique(cnp(c, :))) < 4
            error('Cell %d became degenerate after node merging.', c);
        end
    end

    % ================================================================
    % 3. Build boundary-index lookup from mesh.boundary_faces
    % ================================================================
    %
    % Important:
    %   This lookup is used only for naming/grouping external faces.
    %   It is NOT used to decide whether a face is boundary or internal.
    %
    boundary_lookup = containers.Map('KeyType', 'char', 'ValueType', 'double');

    if isfield(mesh, 'boundary_faces') && isfield(mesh, 'boundary_index') && ...
       ~isempty(mesh.boundary_faces) && ~isempty(mesh.boundary_index)

        if size(mesh.boundary_faces, 2) < 2
            error('mesh.boundary_faces must have at least two node columns.');
        end

        if size(mesh.boundary_faces, 1) ~= numel(mesh.boundary_index)
            error('mesh.boundary_faces and mesh.boundary_index sizes do not match.');
        end

        for k = 1:size(mesh.boundary_faces, 1)
            old_n1 = mesh.boundary_faces(k, 1);
            old_n2 = mesh.boundary_faces(k, 2);

            if old_n1 < 1 || old_n1 > numel(old_to_new) || ...
               old_n2 < 1 || old_n2 > numel(old_to_new)
                error('mesh.boundary_faces contains invalid node index at row %d.', k);
            end

            n1 = old_to_new(old_n1);
            n2 = old_to_new(old_n2);

            key = undirected_edge_key(n1, n2);
            boundary_lookup(key) = mesh.boundary_index(k);
        end

        auto_boundary_index = max(mesh.boundary_index(:)) + 1;
    else
        auto_boundary_index = 1;
    end

    % ================================================================
    % 4. Rebuild faces directly from cell_neighbour_pointer
    % ================================================================
    %
    % cell_neighbour_pointer(c,e) directly identifies the cell on the
    % other side of local edge e.
    %
    % Therefore:
    %   nb > 0  -> internal face, neighbour cell = nb
    %   nb <= 0 -> external boundary face
    %
    % Each internal face is emitted only once using c < nb.
    %
    % This avoids:
    %   1. building a global face/edge map
    %   2. storing every cell-edge as a string key
    %   3. scanning all map keys afterwards
    %
    % Complexity is O(Ncells), with only 4 edges examined per cell.

    % ------------------------------------------------
    % Preallocate
    % ------------------------------------------------
    n_internal_alloc = nnz(cnb > 0);
    n_boundary_alloc = nnz(cnb <= 0);

    internal_faces = zeros(n_internal_alloc, 5);
    boundary_faces = zeros(n_boundary_alloc, 5);

    n_internal_faces = 0;
    n_boundary_faces = 0;

    repaired_one_sided_internal = 0;
    neighbour_mismatch_count    = 0;
    boundary_with_positive_nb   = 0;

    % ------------------------------------------------
    % Traverse cells and use neighbour pointer directly
    % ------------------------------------------------
    for c = 1:n_cells

        cell_nodes = cnp(c, :);

        % Local quad edges:
        %   1: node 1 -> node 2
        %   2: node 2 -> node 3
        %   3: node 3 -> node 4
        %   4: node 4 -> node 1

        for e = 1:4

            [n1, n2] = local_edge_nodes(cell_nodes, e);

            if n1 == n2
                error('Degenerate edge found at cell %d, edge %d.', c, e);
            end

            nb = cnb(c, e);

            % ========================================================
            % Internal face
            % ========================================================
            if nb > 0

                if nb > n_cells
                    error( ...
                        'Invalid neighbour pointer: cell %d edge %d -> %d.', ...
                        c, e, nb);
                end

                % The same internal face is present from both cells.
                % Write it only once.
                if c < nb

                    % Orient face so that owner cell is on the left.
                    [n1o, n2o] = orient_face_with_cell_on_left( ...
                        nodes, cnp, c, n1, n2);

                    n_internal_faces = n_internal_faces + 1;

                    internal_faces(n_internal_faces, :) = ...
                        [n1o, n2o, c, nb, 0];

                    % Optional consistency check only.
                    % Remove this if maximum speed is required.
                    if ~ismember(c, cnb(nb, :))
                        neighbour_mismatch_count = ...
                            neighbour_mismatch_count + 1;
                        repaired_one_sided_internal = ...
                            repaired_one_sided_internal + 1;
                    end
                end

            % ========================================================
            % Boundary face
            % ========================================================
            else

                if nb > 0
                    boundary_with_positive_nb = ...
                        boundary_with_positive_nb + 1;
                end

                key = undirected_edge_key(n1, n2);

                % Orient boundary face.
                [n1o, n2o] = orient_face_with_cell_on_left( ...
                    nodes, cnp, c, n1, n2);

                % Determine boundary group.
                if isKey(boundary_lookup, key)

                    bidx = boundary_lookup(key);

                elseif nb < -1

                    % Optional convention:
                    % negative neighbour stores boundary index.
                    bidx = abs(nb);

                else

                    % Not explicitly listed in boundary_faces.
                    bidx = auto_boundary_index;

                end

                n_boundary_faces = n_boundary_faces + 1;

                boundary_faces(n_boundary_faces, :) = ...
                    [n1o, n2o, c, 0, bidx];
            end
        end
    end

    % ------------------------------------------------
    % Trim preallocated arrays
    % ------------------------------------------------
    internal_faces = internal_faces(1:n_internal_faces, :);
    boundary_faces = boundary_faces(1:n_boundary_faces, :);

    n_total_faces = n_internal_faces + n_boundary_faces;

    % ------------------------------------------------
    % Diagnostics
    % ------------------------------------------------
    if repaired_one_sided_internal > 0
        fprintf( ...
            'Warning: %d internal faces have inconsistent one-sided neighbour pointers.\n', ...
            repaired_one_sided_internal);
    end

    if neighbour_mismatch_count > 0
        fprintf( ...
            'Warning: %d neighbour-pointer inconsistencies detected.\n', ...
            neighbour_mismatch_count);
    end

    if boundary_with_positive_nb > 0
        fprintf( ...
            'Warning: %d boundary faces had positive neighbour pointers.\n', ...
            boundary_with_positive_nb);
    end

    if isempty(boundary_faces)

        warning('No boundary faces found after rebuilding topology.');

        boundary_indices = [];

    else

        boundary_indices = unique(boundary_faces(:, 5));
        boundary_indices = boundary_indices(:)';

    end

    % ================================================================
    % 5. Zone IDs
    % ================================================================
    cell_zone_id     = 1;
    interior_zone_id = 2;

    boundary_zone_ids = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for k = 1:numel(boundary_indices)
        boundary_zone_ids(boundary_indices(k)) = 2 + k;
    end

    % ================================================================
    % 6. Open file
    % ================================================================
    [fid, msg] = fopen(output_file, 'w');
    if fid < 0
        error('Cannot open output file: %s\n%s', output_file, msg);
    end

    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    % ================================================================
    % 7. Write Fluent ASCII mesh
    % ================================================================
    fprintf(fid, '(0 "Generated by export_mesh_to_starccm")\n');
    fprintf(fid, '(2 2)\n\n');

    % -------------------- Nodes --------------------
    fprintf(fid, '(10 (0 1 %s 0))\n', hx(n_nodes));
    fprintf(fid, '(10 (%s 1 %s 1 2)(\n', hx(1), hx(n_nodes));

    for i = 1:n_nodes
        fprintf(fid, ' %.16g %.16g\n', nodes(i, 1), nodes(i, 2));
    end

    fprintf(fid, '))\n\n');

    % -------------------- Cells --------------------
    %
    % Cell section:
    %   type = 1 means fluid
    %   element-type = 3 means quadrilateral
    %
    fprintf(fid, '(12 (0 1 %s 0))\n', hx(n_cells));
    fprintf(fid, '(12 (%s 1 %s 1 3))\n\n', hx(cell_zone_id), hx(n_cells));

    % -------------------- Faces --------------------
    fprintf(fid, '(13 (0 1 %s 0))\n', hx(n_total_faces));

    face_first = 1;

    % Interior faces.
    if n_internal_faces > 0
        face_last = face_first + n_internal_faces - 1;

        % bc-type = 2 interior
        % face-type = 2 linear edge
        fprintf(fid, '(13 (%s %s %s 2 2)(\n', ...
            hx(interior_zone_id), hx(face_first), hx(face_last));

        for i = 1:n_internal_faces
            f = internal_faces(i, :);
            fprintf(fid, ' %s %s %s %s\n', ...
                hx(f(1)), hx(f(2)), hx(f(3)), hx(f(4)));
        end

        fprintf(fid, '))\n\n');

        face_first = face_last + 1;
    end

    % Boundary faces.
    for k = 1:numel(boundary_indices)
        bidx = boundary_indices(k);
        zid  = boundary_zone_ids(bidx);

        mask = boundary_faces(:, 5) == bidx;
        faces_this_zone = boundary_faces(mask, :);

        face_last = face_first + size(faces_this_zone, 1) - 1;

        zone_name = get_boundary_name(boundary_names, bidx, auto_boundary_index);
        [bc_code, ~] = fluent_boundary_type(zone_name);

        fprintf(fid, '(13 (%s %s %s %s 2)(\n', ...
            hx(zid), hx(face_first), hx(face_last), hx(bc_code));

        for i = 1:size(faces_this_zone, 1)
            f = faces_this_zone(i, :);

            % Boundary face has neighbour cell = 0.
            fprintf(fid, ' %s %s %s 0\n', ...
                hx(f(1)), hx(f(2)), hx(f(3)));
        end

        fprintf(fid, '))\n\n');

        face_first = face_last + 1;
    end

    
    % -------------------- Zone names --------------------
    % Important:
    %   In Fluent .msh section 45, zone IDs should be DECIMAL,
    %   not hexadecimal. The usual format is:
    %
    %       (45 (zone-id zone-type zone-name domain-id)())
    %
    %   Grid sections use hexadecimal indices, but section 45 uses decimal IDs.

    domain_id = 1;

    fprintf(fid, '(45 (%d fluid fluid %d)())\n', cell_zone_id, domain_id);

    fprintf(fid, '(45 (%d interior interior %d)())\n', interior_zone_id, domain_id);

    for k = 1:numel(boundary_indices)
        bidx = boundary_indices(k);
        zid  = boundary_zone_ids(bidx);

        zone_name = get_boundary_name(boundary_names, bidx, auto_boundary_index);
        zone_name = sanitize_zone_name(zone_name);

        [~, bc_type_name] = fluent_boundary_type(zone_name);

        fprintf(fid, '(45 (%d %s %s %d)())\n', zid, bc_type_name, zone_name, domain_id);
    end
    % ================================================================
    % 8. Report
    % ================================================================
    %fprintf('\nFluent .msh export complete:\n');
    %fprintf('  File                  : %s\n', output_file);
    %fprintf('  Nodes original        : %d\n', size(nodes_original, 1));
    %fprintf('  Nodes exported        : %d\n', n_nodes);
    %fprintf('  Cells                 : %d\n', n_cells);
    %fprintf('  Internal faces        : %d\n', n_internal_faces);
    %fprintf('  Boundary faces        : %d\n', n_boundary_faces);
    %fprintf('  Boundary zones        : %d\n', numel(boundary_indices));

    if ~isempty(boundary_faces) && any(boundary_faces(:, 5) == auto_boundary_index)
        fprintf('  Note                  : Some true external boundary faces were not found in mesh.boundary_faces.\n');
        fprintf('                          They were exported as "%s".\n', ...
            get_boundary_name(boundary_names, auto_boundary_index, auto_boundary_index));
    end
end

% ========================================================================
% Helper functions for export mesh function
% ========================================================================

function [nodes_new, cnp_new, old_to_new, n_merged] = merge_duplicate_nodes(nodes, cnp, tol)
%MERGE_DUPLICATE_NODES Merge nodes with identical or near-identical coordinates.
%
% Uses a coordinate hash. This is intentionally simple and deterministic.

    n_old = size(nodes, 1);
    old_to_new = zeros(n_old, 1);

    if tol <= 0
        key_scale = 0;
    else
        key_scale = tol;
    end

    origin = min(nodes, [], 1);

    node_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    nodes_new = zeros(n_old, 2);
    n_new = 0;

    for i = 1:n_old
        p = nodes(i, :);

        if key_scale > 0
            q = round((p - origin) ./ key_scale);
            key = sprintf('%.0f,%.0f', q(1), q(2));
        else
            key = sprintf('%.17g,%.17g', p(1), p(2));
        end

        if isKey(node_map, key)
            old_to_new(i) = node_map(key);
        else
            n_new = n_new + 1;
            node_map(key) = n_new;
            old_to_new(i) = n_new;
            nodes_new(n_new, :) = p;
        end
    end

    nodes_new = nodes_new(1:n_new, :);
    cnp_new = old_to_new(cnp);

    n_merged = n_old - n_new;
end

function [n1, n2] = local_edge_nodes(cell_nodes, edge_id)
%LOCAL_EDGE_NODES Return local quad edge nodes.
%
% Quad local edges:
%   1 : node 1 -> node 2
%   2 : node 2 -> node 3
%   3 : node 3 -> node 4
%   4 : node 4 -> node 1

    switch edge_id
        case 1
            n1 = cell_nodes(1);
            n2 = cell_nodes(2);
        case 2
            n1 = cell_nodes(2);
            n2 = cell_nodes(3);
        case 3
            n1 = cell_nodes(3);
            n2 = cell_nodes(4);
        case 4
            n1 = cell_nodes(4);
            n2 = cell_nodes(1);
        otherwise
            error('Invalid local edge id: %d', edge_id);
    end
end

function key = undirected_edge_key(n1, n2)
%UNDIRECTED_EDGE_KEY Hash key for an edge independent of orientation.

    if n1 < n2
        key = sprintf('%d,%d', n1, n2);
    else
        key = sprintf('%d,%d', n2, n1);
    end
end

function [n1o, n2o] = orient_face_with_cell_on_left(nodes, cnp, cell_id, n1, n2)
%ORIENT_FACE_WITH_CELL_ON_LEFT Orient edge so owner cell centroid is left of edge.

    p1 = nodes(n1, :);
    p2 = nodes(n2, :);

    cell_pts = nodes(cnp(cell_id, :), :);
    centroid = mean(cell_pts, 1);

    edge_vec = p2 - p1;
    mid_pt   = 0.5 * (p1 + p2);
    to_cell  = centroid - mid_pt;

    cross_val = edge_vec(1) * to_cell(2) - edge_vec(2) * to_cell(1);

    if cross_val >= 0
        n1o = n1;
        n2o = n2;
    else
        n1o = n2;
        n2o = n1;
    end
end

function s = hx(v)
%HX Convert integer to Fluent-style hexadecimal string.

    v = round(double(v));

    if v < 0
        error('Cannot convert negative value to Fluent hex index.');
    end

    s = lower(dec2hex(v));
end

function name = get_boundary_name(boundary_names, bidx, auto_boundary_index)
%GET_BOUNDARY_NAME Return safe boundary name from list or fallback.

    if isstring(boundary_names)
        boundary_names = cellstr(boundary_names);
    end

    if ~isempty(boundary_names) && ...
       bidx >= 1 && ...
       bidx <= numel(boundary_names) && ...
       ~isempty(boundary_names{bidx})

        name = boundary_names{bidx};

    elseif bidx == auto_boundary_index
        name = 'unassigned_boundary';

    else
        name = sprintf('boundary_%d', bidx);
    end

    name = sanitize_zone_name(name);
end

function name = sanitize_zone_name(name)
%SANITIZE_ZONE_NAME Make a STAR-CCM+/Fluent-safe zone name.

    if isstring(name)
        name = char(name);
    end

    if isempty(name)
        name = 'boundary';
    end

    name = regexprep(name, '\s+', '_');
    name = regexprep(name, '[^A-Za-z0-9_]', '_');

    if isempty(regexp(name(1), '[A-Za-z_]', 'once'))
        name = ['zone_', name];
    end
end

function [bc_code, bc_type_name] = fluent_boundary_type(zone_name)
%FLUENT_BOUNDARY_TYPE Pick a conservative Fluent boundary type.
%
% Common Fluent face-zone boundary codes:
%   2  interior
%   3  wall
%   5  pressure-outlet
%   7  symmetry
%   10 velocity-inlet
%
% Periodic patches are exported as wall because Fluent periodic zones need
% paired periodic/shadow data. Convert them inside STAR-CCM+ if needed.

    z = lower(zone_name);

    if strcmp(z, 'inlet') || contains(z, 'velocity_inlet')
        bc_code = 10;
        bc_type_name = 'velocity-inlet';

    elseif strcmp(z, 'outlet') || contains(z, 'pressure_outlet')
        bc_code = 5;
        bc_type_name = 'pressure-outlet';

    elseif contains(z, 'symmetry')
        bc_code = 7;
        bc_type_name = 'symmetry';

    else
        bc_code = 3;
        bc_type_name = 'wall';
    end
end




























% ========================================================================
% integrated functions
% ========================================================================

function [mesh, output_file] = generate_stator_mesh( ...
    alpha, beta, gamma, chord_length, ...
    leading_factor, trailing_factor, ...
    relative_maxium_thickness, relative_t_maxium_thickness, ...
    relative_leading_edge_radius, relative_trailing_edge_radius, ...
    blade_number)

    % Physical parameters
    mass_flow = 0.4211;
    hub_diam = 0.1437;
    passage_radial_size = 0.0058;
    shroud_diam = hub_diam + passage_radial_size * 2;
    outlet_area = pi * (shroud_diam^2 - hub_diam ^2) / 4;
    outlet_mass_flux = mass_flow / outlet_area;

    mid_diam = hub_diam + passage_radial_size;
    inlet_velocity_angle = deg2rad(71.9959);

    % Stator parameters are supplied as function inputs

    pitch = pi * mid_diam / blade_number; % distance between the blades
    inlet_extension = 0.015;
    outlet_extension = 0.09;

    % Generate mesh
    mesh = generate_fluid_domain_mesh( ...
       alpha, beta, gamma, chord_length, ...
        leading_factor, trailing_factor, ...
        relative_maxium_thickness, relative_t_maxium_thickness, ...
        relative_leading_edge_radius, relative_trailing_edge_radius, ...
        pitch, inlet_extension, outlet_extension);

    % Boundary names
    boundary_names = {
        'inlet_period_1'
        'suction_wall'
        'te_suction'
        'outlet_period_1'
        'inlet'
        'outlet'
        'inlet_period_2'
        'outlet_period_2'
        'le_round'
        'pressure_wall'
        'te_pressure'
    };

    % Export mesh
    output_file = fullfile(pwd, 'stator_passage.msh');

    export_mesh_to_starccm(mesh, output_file, boundary_names);

    %fprintf('Mesh generated successfully.\n');
    %fprintf('Blade number: %d\n', blade_number);
    %fprintf('Pitch: %.8f m\n', pitch);
    %fprintf('Output file: %s\n', output_file);

end














% ========================================================================
% optmizer functions
% ========================================================================

function val = getLastMonitorValue(filename)
% GETLASTMONITORVALUE  Return the second column of the last row of a
% monitor CSV as a scalar.
%
%   val = getLastMonitorValue('Mass Flow Averaged 1 Monitor Plot.csv')
%
%   The file is expected to have a header line and rows of the form:
%       iteration, value
%   Only the tail of the file is read, so this is cheap even for large files.

    if nargin < 1 || isempty(filename)
        filename = 'Mass Flow Averaged 1 Monitor Plot.csv';
    end

    fid = fopen(filename, 'r');
    if fid < 0
        error('getLastMonitorValue:openFailed', ...
              'Could not open file: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));   %#ok<NASGU>

    % Jump to end and read a small tail chunk
    fseek(fid, 0, 'eof');
    fileSize = ftell(fid);
    chunkSize = min(4096, fileSize);
    fseek(fid, -chunkSize, 'eof');
    chunk = fread(fid, chunkSize, '*char').';

    % Split into lines, discard empty ones
    lines = regexp(chunk, '\r?\n', 'split');
    lines = lines(~cellfun(@(s) isempty(strtrim(s)), lines));
    if isempty(lines)
        error('getLastMonitorValue:emptyFile', 'File is empty: %s', filename);
    end

    % Parse the last line
    lastLine = strtrim(lines{end});
    vals = str2double(strsplit(lastLine, ','));
    if numel(vals) < 2 || any(isnan(vals(1:2)))
        error('getLastMonitorValue:parseFailed', ...
              'Could not parse numeric values from last line: "%s"', lastLine);
    end

    val = vals(2);   % scalar
end

function solution = start_simulation(starExe, macroFile)

    if nargin < 2 || isempty(starExe) ||isempty(macroFile)
        starExe   = 'D:\Program Files\Siemens\20.06.010-R8\STAR-CCM+20.06.010-R8\star\lib\win64\clang17.0vc14.2-r8\lib\starccmw.exe';
        macroFile   = 'E:\temp\blue\difusor\optimizer test\setting_up.java';
    end



    % Number of parallel processes 
    nProcs = 18;
    cmd = sprintf('"%s" -new -np %d -macro "%s"',starExe, nProcs, macroFile);
    %cmd = sprintf('"%s" -new -macro "%s"', starExe, macroFile);

    [status, cmdout] = system(cmd);
    disp(cmdout);
    if status ~= 0
        disp("abnormal exiting");
    end
    solution = [getLastMonitorValue('Mass Flow Averaged 1 Monitor Plot.csv'); getLastMonitorValue('Mass Flow 1 Monitor Plot.csv')];
end

function x = optmizer_vector_translation(...
    alpha, beta, gamma, chord_length, ...
    leading_factor, trailing_factor, ...
    relative_maxium_thickness, relative_t_maxium_thickness, ...
    relative_leading_edge_radius, relative_trailing_edge_radius, ...
    blade_number)

    x = [ ...
    alpha, beta, gamma, chord_length, ...
    leading_factor, trailing_factor, ...
    relative_maxium_thickness, relative_t_maxium_thickness, ...
    relative_leading_edge_radius, relative_trailing_edge_radius, ...
    blade_number];
end

function objective = objective_blade(x)

    % Extract design variables
    %disp(x);
    alpha = x(1);
    beta = x(2);
    gamma = x(3);

    chord_length = x(4);
    leading_factor = x(5);
    trailing_factor = x(6);

    relative_maxium_thickness = x(7);
    relative_t_maxium_thickness = x(8);

    relative_leading_edge_radius = x(9);
    relative_trailing_edge_radius = x(10);

    blade_number = x(11);








    outlet_extension = 0.09;








    generate_stator_mesh( ...
        alpha, beta, gamma, chord_length, ...
        leading_factor, trailing_factor, ...
        relative_maxium_thickness, ...
        relative_t_maxium_thickness, ...
        relative_leading_edge_radius, ...
        relative_trailing_edge_radius, ...
        blade_number);


        



        


        solution = start_simulation();
        evaluation = solution(1) * 154^2 / ((solution(2) / outlet_extension)^2 + 154^2); 
        fprintf('%g, %g\n', transpose(solution));










        objective.Fval = 350000 - evaluation;
        objective.Ineq = (x(4) * x(7)) / (0.1495 * pi / x(11)) - 0.5;
end


clc;
clear;

% Design variable vector:
% x = [alpha;
%      beta;
%      gamma;
%      chord_length;
%      leading_factor;
%      trailing_factor;
%      relative_maxium_thickness;
%      relative_t_maxium_thickness;
%      relative_leading_edge_radius;
%      relative_trailing_edge_radius;
%      blade_number];

%% Individual lower and upper bounds

lb = [ ...
    deg2rad(65);              % alpha
    -pi/36;          % beta
    deg2rad(15);    % gamma
    0.04;           % chord_length
    0.3;            % leading_factor
    0.2;            % trailing_factor
    0.06;           % relative_maxium_thickness
    0.28;            % relative_t_maxium_thickness
    0.05;           % relative_leading_edge_radius
    0.025;           % relative_trailing_edge_radius
    16];            % blade_number

ub = [ ...
    deg2rad(75);           % alpha
    0;              % beta
    pi/4;           % gamma
    0.05;           % chord_length
    0.45;            % leading_factor
    0.34;            % trailing_factor
    0.12;            % relative_maxium_thickness
    0.5;            % relative_t_maxium_thickness
    0.05;            % relative_leading_edge_radius
    0.025;           % relative_trailing_edge_radius
    24];            % blade_number

%% Integer constraint

intcon = 11;

%% Linear inequality constraints: A*x <= b

A = zeros(3,11);

% gamma <= alpha
A(1,3) = 1;
A(1,1) = -1;

% gamma >= 0
% Already enforced by lb(3) = 0

% beta <= 0
% Already enforced by ub(2) = 0

% gamma >= 0 and alpha >= gamma are enforced by:
% gamma <= alpha, lb(3) = 0

% beta >= -pi/8 is enforced by lb(2)

b = zeros(3,1);

%% Optimization options

options = optimoptions('surrogateopt', ...
    'MaxFunctionEvaluations', 40, ...
    'MaxTime', Inf, ...
    'Display', 'iter', ...
    'PlotFcn', [], ...
    'UseParallel', false, ...
    'UseVectorized', false);

%% Run optimization

[x_opt, fval, exitflag, output] = surrogateopt( ...
    @objective_blade, ...
    lb, ub, intcon, A, b, [], [], options);

%% Display results

disp('Optimal design variables:');
disp(x_opt);

disp('Minimum objective value:');
disp(fval);