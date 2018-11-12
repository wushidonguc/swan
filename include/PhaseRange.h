//

#include "header.h"
#include "OneDimRange.h"
#include "Parameters.h"

/* PhaseRange class */

class PhaseRange : public OneDimRange
{
  public:
    PhaseRange(double start, double end, int n): OneDimRange(start, end, n){}
    PhaseRange(double start, double end, double spacing): OneDimRange(start, end, spacing){}    
    PhaseRange(Parameters parameters, int rank, int nprocs): 
      OneDimRange(0.5 * M_PI * ((double)parameters.mNKPointSampling / (double)nprocs * (double)rank), 0.5 * M_PI * ((double)parameters.mNKPointSampling / (double)nprocs * (double)(rank + 1)), (parameters.mNKPointSampling)){}
    ~PhaseRange();
};
