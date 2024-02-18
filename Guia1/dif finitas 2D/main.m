close all; clear all; more off;

xnode = [
  -4, -3;
  -2, -3;
   0, -3;
   2, -3;
   4, -3;
   -4, -1;
   -2, -1;
   0, -1;
   2, -1;
   4, -1;
   -4, 1;
   -2, 1;
   0, 1;
   2, 1;
   4, 1;
   -4, 3;
   -2, 3;
   0, 3;
   2, 3;
   4, 3
];

icone = [
       1,2,7,6;
       2,3,8,7;
       3,4,9,8;
       4,5,10,9;
       6,7,12,11;
       7,8,13,12;
       8,9,14,13;
       9,10,15,14;
       11,12,17,16;
       12,13,18,17;
       13,14,19,18;
       14,15,20,19;
];

DIR = [
       1, 20.0000000000000000;
       6, 20.0000000000000000;
       11, 20.0000000000000000;
       16, 20.0000000000000000;
];

NEU = [
       16, 0, 3;
       17, 0, 3;
       18, 0, 3;
       19, 0, 3;
       20, 0, 3;
       1, 0, 1;
       2, 0, 1;
       3, 0, 1;
       4, 0, 1;
       5, 0, 1;
];

ROB = [
       20, 2, 10, 2;
       15, 2, 10, 2;
       10, 2, 10, 2;
        5, 2, 10, 2;
];

disp('---------------------------------------------------------------');
disp('Inicializando modelo de datos...');

model.nnodes = size(xnode,1);

model.k = [
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
    3.0000000000000000;
];

model.c = [
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
];

model.G = [
    100*xnode(1);
    100*xnode(2);
    100*xnode(3);
    100*xnode(4);
    100*xnode(5);
    100*xnode(6);
    100*xnode(7);
    100*xnode(8);
    100*xnode(9);
    100*xnode(10);
    100*xnode(11);
    100*xnode(12);
    100*xnode(13);
    100*xnode(14);
    100*xnode(15);
    100*xnode(16);
    100*xnode(17);
    100*xnode(18);
    100*xnode(19);
    100*xnode(20);
];

% Esquema Temporal: [0] Explícito, [1] Implícito, [X] Estacionario
model.ts = 2;

% Parámetros para esquemas temporales
model.rho = 1.0000000000000000;
model.cp = 1.0000000000000000;
model.maxit =            1;
model.tol = 1.000000e-05;

% Condición inicial
model.PHI_n = mean(DIR(:,2))*ones(model.nnodes,1);

disp('Iniciando el método numérico...');

% Llamada principal al Método de Diferencias Finitas
[PHI,Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

disp('Finalizada la ejecución del método numérico.');

disp('Temperatura en el nodo 8')
disp(PHI(8))
disp('Temperatura en el nodo 13')
disp(PHI(13))

disp('---------------------------------------------------------------');
disp('Iniciando el post-procesamiento...');

% mode ---> modo de visualización:
%           [0] 2D - Con malla
%           [1] 3D - Con malla
%           [2] 2D - Sin malla
%           [3] 3D - Sin malla
% graph --> tipo de gráfica:
%           [0] Temperatura (escalar)
%           [1] Flujo de Calor (vectorial)
%           [2] Flujo de Calor eje-x (escalar)
%           [3] Flujo de Calor eje-y (escalar)
%           [4] Magnitud de Flujo de Calor (escalar)
mode = 0;
graph = 0;
fdm2d_graph_mesh(full(PHI),Q,xnode,icone,mode,graph);

disp('Finalizado el post-procesamiento.');
