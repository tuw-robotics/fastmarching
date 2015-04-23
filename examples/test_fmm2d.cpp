/* An example to run FMM on a 3D grid loaded from a given text file. */

#include <iostream>
#include <array>

#include "../fmm/fmdata/fmcell.h"
#include "../ndgridmap/ndgridmap.hpp"

#include "../fmm/fmm.hpp"
#include "../fmm/fim.hpp"
#include "../io/maploader.hpp"
#include "../io/gridplotter.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, char **argv)
{
    // A bit of shorthand.
    typedef nDGridMap<FMCell, 2> FMGrid2D;
    typedef array<unsigned int, 2> Coord2D;

    // Grid, start and goal definition.
    FMGrid2D grid_fmm;
    MapLoader::loadMapFromText(argv[1], grid_fmm);
    Coord2D init_point = {1000,1000};

    // Solvers declaration.
    std::vector<Solver<FMGrid2D>*> solvers;
    solvers.push_back(new FMM<FMGrid2D>);
    solvers.push_back(new FIM<FMGrid2D>);

    // Executing every solver individually over the same grid.
    for (Solver<FMGrid2D>* s :solvers)
    {
        s->setEnvironment(&grid_fmm);
        s->setInitialPoints(init_point);
        s->compute();
        cout << "\tElapsed "<< s->getName() <<" time: " << s->getTime() << " ms" << '\n';
        //GridWriter::saveGridValues("3dresult.grid", grid_fmm);
        GridPlotter::plotArrivalTimes(grid_fmm);
    }

    return 0;
}
