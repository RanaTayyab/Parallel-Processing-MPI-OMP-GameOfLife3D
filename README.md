# GameOfLife3D

The Game of Life was at start build in two dimensional grid of squares cells, each of which is in one of two possible states, live or dead (or populated and unpopulated, respectively). Every cell interacts with its eight neighbours, which are the cells that are horizontally, vertically, or diagonally adjacent.

# Rules
• Any live cell with fewer than two live neighbours dies, as if by underpopulation.

• Any live cell with two or three live neighbours lives on to the next generation.

• Any live cell with more than three live neighbours dies, as if by overpopulation.

• Any dead cell with exactly three live neighbours becomes a live cell, as if by reproduction.

Implementation involves the MPI and Openmp functions
