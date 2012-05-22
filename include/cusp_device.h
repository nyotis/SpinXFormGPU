// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- cusp_device.h
// Nikos Yiotis
//

#ifndef CUSP_DEVICE_H
#define	CUSP_DEVICE_H

//#pragma once

#include <cusp/coo_matrix.h>
    
// function prototype
void solve_on_device( cusp::coo_matrix<int, float, cusp::host_memory>& coo_host, 
                      cusp::array1d<float, cusp::host_memory>&         rhs_host,
                      cusp::array1d<float, cusp::host_memory>&         result_host );


#endif	/* CUSP_DEVICE_H */

