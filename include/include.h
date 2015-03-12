//nog enkele definities:
#include <btas/common/blas_cxx_interface.h>

#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; }; // Define SpinQuantum as default quantum class

#include <btas/QSPARSE/QSDArray.h>

#include "MPSblas.h"

#include "global.h"
#include "dmrg.h"
#include "driver.h"
#include "davidson.h"

#include "SpinHamiltonian.h"
#include "Random.h"
#include "coupling.h"
