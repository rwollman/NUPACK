#pragma once
#include "Bind.h"
#include <nupack/math/BoundSolve.h>
#include <nupack/math/Sparse.h>
#include <nupack/types/Fenwick.h>
#include <nupack/types/Matrix.h>
#include <nupack/concentration/Solve.h>

namespace nupack {

/******************************************************************************/

namespace concentration {

rebind::Variable response(std::type_index, Options const &o);
std::optional<Options> request(Type<Options>, rebind::Variable const &, rebind::Dispatch &msg);

}

/******************************************************************************/

}
