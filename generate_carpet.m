%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem parameters
L = 0.08;                 % computational domain linear dimension (cm)
ASPECT_RATIO_X = 1.0;     % ratio of computational domain length and characteristic linear dimension (dimensionless)
ASPECT_RATIO_Y = 1.0;     % ratio of computational domain width  and characteristic linear dimension (dimensionless)
ASPECT_RATIO_Z = 0.25;    % ratio of computational domain height and characteristic linear dimension (dimensionless)

post_density = 100000.0;  % posts per unit area (cm^{-2})
post_height = 0.005;      % post height (cm)
post_r0 = 0.0002;         % radius of unstressed rod (cm)
E = 1.0e4;                % elastic modulus (dyne cm^{-2}) [random choice!]
A = pi*post_r0^2;         % cross-sectional area (cm^2)
I = 0.25*pi*post_r0^4;    % second moment of area (cm^4)
kappa_s = E*A;            % spring stiffness (dyne)
kappa_b = E*I;            % bending stiffness (dyne cm^2)

n_posts = ceil(post_density * (ASPECT_RATIO_X * L) * (ASPECT_RATIO_Y * L));

N = 128;                  % number of grid cells on finest grid level
h = L/N;                  % Cartesian grid spacing

n_ib_post = ceil(0.75 * post_height / h);  % number of IB points per post
n_ib = n_ib_post * n_posts;                % total number of IB points
dX = post_height / (n_ib_post-1);          % IB point spacing along the post

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertex_fid = fopen(['carpet_' num2str(N) '.vertex'], 'w');
spring_fid = fopen(['carpet_' num2str(N) '.spring'], 'w');
beam_fid   = fopen(['carpet_' num2str(N) '.beam']  , 'w');
anchor_fid = fopen(['carpet_' num2str(N) '.anchor'], 'w');

fprintf(vertex_fid, '%d\n', n_posts * n_ib_post);
fprintf(spring_fid, '%d\n', n_posts * (n_ib_post - 1));
fprintf(beam_fid, '%d\n', n_posts * (n_ib_post - 2));
fprintf(anchor_fid, '%d\n', n_posts);

for p = 0:n_posts-1

  % vertices:
  X(1) = rand() * ASPECT_RATIO_X * L;
  X(2) = rand() * ASPECT_RATIO_Y * L;
  for l = 0:n_ib_post-1
    fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n', X(1), X(2), l*dX);
  end

  % springs:
  for l = 0:n_ib_post-2
    fprintf(spring_fid, '%6d %6d %1.16e %1.16e\n', p*n_ib_post + l, p*n_ib_post + l + 1, kappa_s/dX, dX);
  end

  % beams:
  for l = 0:n_ib_post-3
    fprintf(beam_fid, '%6d %6d %6d %1.16e\n', p*n_ib_post + l, p*n_ib_post + l + 1, p*n_ib_post + l + 2, kappa_b/dX^3);
  end

  % anchor points:
  fprintf(anchor_fid, '%6d\n', p*n_ib_post);
end

fclose(vertex_fid);
fclose(spring_fid);
fclose(beam_fid);
fclose(anchor_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
