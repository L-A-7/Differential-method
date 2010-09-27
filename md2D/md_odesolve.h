/* 	\file  md_odesolve.h
 *  \brief Header file for md_odesolve
 */

#ifndef _MD_ODESOLVE_H
#define _MD_ODESOLVE_H

#include "std_include.h"

#include "md_io_utils.h"
#include "md_utils.h"
#include "md_maths.h"

int fun (double z, const double *F_reel, double *dF_reel, void *param_void);



#endif /* _MD_ODESOLVE_H */

