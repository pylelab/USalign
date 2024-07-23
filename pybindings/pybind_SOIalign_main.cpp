//essential preamble for pybind11 casting
#include <pybind11/pybind11.h>

namespace py = pybind11;
//

#include "SOIalign_main.hpp" // point to the SOIalign_main header file

// defining the pybind11 wrapper for SOIalign_main
PYBIND11_MODULE(SOIalign_main, m) {
    m.doc() = "pybind11 port of SOIalign_main function from USalign code base"; 
    m.def("SOIalign_main", // define the function name
          &SOIalign_main,  // set return values, I think (?)
	  "function to calculate the Sequence-Order-Independent alignment between two sets of coordinates", // blurb about what the function does
	  py::arg("mobile_coords") = , // maps to xa variable
	  py::arg("target_coords") = , // maps to ya variable
	  py::arg("mobile_neighbor_list") = ,	// maps to xk variable
	  py::arg("target_neighbor_list") = ,	// maps to yk variable
	  py::arg("n_nearest_neighbors") = 3,	// maps to closeK_opt variable
	  py::arg("mobile_seq") = , // maps to seqx variable
	  py::arg("target_seq") = , // maps to seqy variable
	  py::arg("mobile_sec") = , // maps to secx variable
	  py::arg("target_sec") = , // maps to secy variable
	  // variables to be changed internal to the function... why are we passing these into the code???
	  py::arg("translation_vector") = , // maps to t0 variable
	  py::arg("rotation_matrix") = , // maps to u0 variable
	  py::arg("TM1") = , // maps to TM1 variable
	  py::arg("TM2") = , // maps to TM1 variable
	  py::arg("d0_0") = , // maps to d0_0 variable
	  py::arg("TM_0") = , // maps to TM_0 variable
	  py::arg("d0A") = , // maps to d0A variable
	  py::arg("d0B") = , // maps to d0A variable
	  py::arg("aligned_seq") = , // maps to seqM variable
	  py::arg("mobile_aligned_seq") = , // maps to seqxA variable
	  py::arg("target_aligned_seq") = , // maps to seqyA variable
	  py::arg("inverse_map") = , // maps to invmap variable
	  py::arg("final_rmds") = , // maps to rmsd0 variable
	  py::arg("n_aligned_residues") = , // maps to n_ali8 variable
	  py::arg("mobile_len") = , // maps to xlen variable
	  py::arg("target_len") = , // maps to ylen variable
	  py::arg("distance_list") = , // maps to dist_list variable
	  py::arg("mobile_sec_bond") = , // maps to secx_bond variable
	  py::arg("target_sec_bond") = , // maps to secy_bond variable
	  py::arg("mm_opt") = , // maps to variable mm_opt
	  //py::arg("") = , //
	  py::return_value_policy::copy);
}
