function temperature = get_temperature(density, pressure)
    %GET_TEMPERATURE Calcula la temperatura a partir de la ecuacion de los gases ideales

    temperature = pressure / Air.R / density;
end
