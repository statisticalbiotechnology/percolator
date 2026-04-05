#include "MonotoneRegressor.h"
#include "PAVARegressor.h"
#include "ISplineTRRRegressor.h"
// Factory
// Forward-defined concrete classes

std::unique_ptr<MonotoneRegressor> make_monotone_regressor(RegressorType type) {
  std::unique_ptr<MonotoneRegressor> ptr;
  switch (type) {
    case RegressorType::PAVA:
      ptr = std::unique_ptr<MonotoneRegressor>(new PAVARegressor());
      break;
    case RegressorType::ISPLINE_TRR:
      ptr = std::unique_ptr<MonotoneRegressor>(new ISplineTRRRegressor());
      break;
  }
  ptr->set_params(MonotoneParams());
  return ptr;
}
