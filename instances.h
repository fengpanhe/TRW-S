#ifndef __INSTANCES_H__
#define __INSTANCES_H__

#if defined(_MSC_VER)

// C4661: '...' : no suitable definition provided for explicit template instantiation request
#pragma warning(disable : 4661)

#endif

#include "typeBinary.h"
#include "typeBinaryFast.h"
#include "typeGeneral.h"
#include "typePotts.h"
#include "typePotts2.h"
#include "typeTruncatedLinear.h"
#include "typeTruncatedLinear2D.h"
#include "typeTruncatedQuadratic.h"
#include "typeTruncatedQuadratic2D.h"

#endif
