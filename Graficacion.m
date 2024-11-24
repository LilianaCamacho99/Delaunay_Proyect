% Leer archivo .NODE
fileNode = 'acapulco.node'; % Nombre del archivo
fidNode = fopen(fileNode, 'r');
nodeHeader = fscanf(fidNode, '%d %d %d %d', [1, 4]); % Leer encabezado
dataNode = fscanf(fidNode, '%d %f %f %d', [4, nodeHeader(1)])'; % Leer datos
fclose(fidNode);

nodeIDs = dataNode(:, 1); % ID de los nodos
xCoords = dataNode(:, 2); % Coordenadas X
yCoords = dataNode(:, 3); % Coordenadas Y

% Leer archivo .ELE
fileEle = 'acapulco.ele'; % Nombre del archivo
fidEle = fopen(fileEle, 'r');
eleHeader = fscanf(fidEle, '%d %d %d', [1, 3]); % Leer encabezado
dataEle = fscanf(fidEle, '%d %d %d %d', [4, eleHeader(1)])'; % Leer datos
fclose(fidEle);

elementIDs = dataEle(:, 1); % ID de los elementos
connectivity = dataEle(:, 2:4); % Conectividad de los triángulos

% Graficar la malla triangular en 2D
figure;
triplot(connectivity, xCoords, yCoords, 'k'); % Función para graficar triángulos en 2D
xlabel('X');
ylabel('Y');
title('Malla triangular de Acapulco (2D)');
axis equal;
