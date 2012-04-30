#include <cusp/coo_matrix.h>
#include <cusp/hyb_matrix.h>
#include <cusp/dia_matrix.h>
#include <cusp/ell_matrix.h>
#include <cusp/print.h>
#include <cusp/monitor.h>
#include <cusp/krylov/cg.h>
#include <cusp/precond/diagonal.h>

#include <cusp/io/matrix_market.h> // save matrix to disk in MatrixMarket format
#include "./include/cusp_device.h"

//cudaError_t error; 

void solve_on_device( cusp::coo_matrix<int, float, cusp::host_memory>& coo_host, 
                      cusp::array1d<float, cusp::host_memory>& rhs_host, 
                      cusp::array1d<float, cusp::host_memory>& result_host ) {
    	 															
    // convert COO->HYB on the host and transfer HYB to the device
    //cusp::hyb_matrix<int, float, cusp::device_memory> coo_cusp_device = coo_host;
    
    // transfer COO to the device
    cusp::coo_matrix<int, float, cusp::device_memory> coo_cusp_device = coo_host;
    //cusp::io::write_matrix_market_file( coo_cusp_device, "coo_cusp_device.mtx" );
    //cusp::print(coo_cusp_device);

    // convert COO->DIA on the host and transfer DIA to the device
    //cusp::dia_matrix<int, float, cusp::device_memory> dia_cusp_device = coo_host;
 
    // convert COO->DIA on the device
    //cusp::dia_matrix<int, float, cusp::device_memory> dia_cusp_device;
    //cusp::convert(dia_cusp_device, coo_cusp_device);
	
	 // transfer rhs_host to the device
    cusp::array1d<float, cusp::device_memory> rhs_device = rhs_host;
    //cusp::print(rhs_device);
    //cusp::io::write_matrix_market_file( rhs_device, "rhs_device.mtx" );

	 // transfer result_host to the device	
    cusp::array1d<float, cusp::device_memory> result_device = result_host;
         
    // set stopping criteria (iteration_limit = 2000, relative_tolerance = 1e-7)
    cusp::verbose_monitor<float> monitor(rhs_device, 100, 1e-2);
    
    // set preconditioner (identity) doesn't affect the speed of convergence
    cusp::identity_operator<float, cusp::device_memory> M( coo_cusp_device.num_rows, 
																			  coo_cusp_device.num_rows );

    // solve the linear system A * x = b -> coo_cusp_device * result_device = rhs_device 
    cusp::krylov::cg(coo_cusp_device, result_device, rhs_device, monitor, M);	 
    //cusp::print(x);
    //cusp::io::write_matrix_market_file(result_device, "result_device.mtx");
    
    // return the result on the host 
    // behind each '=' there is a call to cudamalloc
    result_host = result_device; 
    
    //cusp::io::write_matrix_market_file(result_host, "result_host_final.mtx");

}
    // diagonal preconditioner results in NaN
    // cusp::precond::diagonal<float, cusp::device_memory> M( coo_cusp_device ); 
																																					  																		  																		  
