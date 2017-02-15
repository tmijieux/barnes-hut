#include <mpi.h>
#include "proc.hpp"

using namespace barnes_hut;

proc::proc()
{
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_group_size);
}

proc::~proc()
{
    MPI_Finalize();
}

void proc::print(std::ostream& out) const
{
    out << "proc: rank=" << _rank << " | group_size=" << _group_size<<std::endl;
}

void proc::print() const
{
    print(std::cout);
}
