function [scattering_region] = define_map_region(Nz, Nx, num, rad, rad_inclusion, center_inclusion)
scattering_region = zeros(Nz,Nx);
[X,Z] = meshgrid(1:Nx,1:Nz);
distance = sqrt((Z - center_inclusion(1)).^2 + (X - center_inclusion(2)).^2);
available_region = distance <= rad_inclusion - rad;
centers = generate_random_positions(available_region, num, rad, [Nx,Nz]);
scattering_region = assign_circles(scattering_region, centers, rad, 1);
end
function positions = generate_random_positions(available_region, num_circles, radius, grid_size)
    % Encontrar las coordenadas disponibles
    % axial lateral
    [y, x] = find(available_region);
    available_coords = [x, y];
    
    positions = zeros(num_circles, 2);
    count = 1;

    % Intentar generar posiciones aleatorias para cada círculo
    while count <= num_circles
        % Seleccionar una posición aleatoria de las coordenadas disponibles
        random_index = randi(size(available_coords, 1));
        position = available_coords(random_index, :);

        % Verificar que la nueva posición no se superponga con círculos existentes
        if count == 1 || all(sqrt((position(1) - positions(1:count - 1, 1)).^2 + (position(2) - positions(1:count - 1, 2)).^2) >= 2 * radius)
            % Verificar que la nueva posición esté dentro de las regiones disponibles
            if position(1) - radius > 0 && position(1) + radius <= grid_size(1) && position(2) - radius > 0 && position(2) + radius <= grid_size(2)
                % Asignar la nueva posición si cumple todas las condiciones
                positions(count, :) = position;
                count = count + 1;
            end
        end
    end
end

function scattering_region = assign_circles(scattering_region, centers, radius, value)
    % Crear una cuadrícula de coordenadas
    %                         lateral                         axial
    [X, Z] = meshgrid(1:size(scattering_region, 2), 1:size(scattering_region, 1));

    % Asignar círculos en la región de dispersión
    for i = 1:size(centers, 1)
        x0 = centers(i, 1); %lateral
        z0 = centers(i, 2); %axial

        % Calcular la distancia desde el centro del círculo
        dist = sqrt((Z - z0).^2 + (X - x0).^2);
        % Asignar valores dentro del círculo
        scattering_region(dist <= radius) = value;
    end
end