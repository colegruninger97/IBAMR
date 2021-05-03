%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem parameters
L = 0.08;                 % computational domain linear dimension (cm)
ASPECT_RATIO_X = 1.0;     % ratio of computational domain length and characteristic linear dimension (dimensionless)
ASPECT_RATIO_Y = 1.0;     % ratio of computational domain width  and characteristic linear dimension (dimensionless)
ASPECT_RATIO_Z = 0.25;    % ratio of computational domain height and characteristic linear dimension (dimensionless)

particle_density = 1.0e5; % particles per unit volume (cm^{-3})
n_particles = ceil(particle_density * (ASPECT_RATIO_X * L) * (ASPECT_RATIO_Y * L) * (ASPECT_RATIO_Z * L));

NFINEST = 128;            % number of grid cells on finest grid level
h = L/NFINEST;            % Cartesian grid spacing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertex_fid = fopen(['particles_' num2str(NFINEST) '.vertex'], 'w');

fprintf(vertex_fid, '%d\n', n_particles);

for p = 0:n_particles-1

  % vertices:
  X(1) = rand() * ASPECT_RATIO_X * L;
  X(2) = rand() * ASPECT_RATIO_Y * L;
  X(3) = rand() * ASPECT_RATIO_Z * L;
  fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n', X(1), X(2), X(3));

end

fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
