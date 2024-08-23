//essential preamble for pybind11 casting
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
//#include <cmath>
#include <vector>
#include <string>
//#include <algorithm>

//#include "SOIalign.h"	// 

namespace py = pybind11;


/****************************************************************************
 * FUNCTIONS NOT TO PORT TO PYTHON
 * these are still important to the the ports
 * they are called within the ported code
 ****************************************************************************/

/* creates an array object of type template class A with length of Narray1
 * then loops over those rows to create a subarray of type template class A
 * of length Narray2
 * PULLED DIRECTLY FROM USALIGN basic_fun.cpp file
 */ 
template <class A> void NewArray(A *** array, int Narray1, int Narray2)
{
    *array=new A* [Narray1];
    for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
}

// clean up the 2d matrix; have to loop through rows (Narray) since each row is 
// itself an array object
// PULLED DIRECTLY FROM USALIGN basic_fun.cpp file
template <class A> void DeleteArray(A *** array, int Narray)
{
    for(int i=0; i<Narray; i++)
        if(*(*array+i)) delete [] *(*array+i);
    if(Narray) delete [] (*array);
    (*array)=NULL;
}

// need to include the dist function from basic_fun.h for make_sec()
// actually returns the dist**2
double dist(double x[3], double y[3])
{
    double d1=x[0]-y[0];
    double d2=x[1]-y[1];
    double d3=x[2]-y[2];
 
    return (d1*d1 + d2*d2 + d3*d3);
}

/* secondary structure assignment for protein:
 * 1->coil, 2->helix, 3->turn, 4->strand 
 * for the set of distances provided, return a 2ndary structure character
 */
char sec_str(double dis13, double dis14, double dis15,
            double dis24, double dis25, double dis35)
{
    char s='C';
    
    double delta=2.1;
    if (fabs(dis15-6.37)<delta && fabs(dis14-5.18)<delta && 
        fabs(dis25-5.18)<delta && fabs(dis13-5.45)<delta &&
        fabs(dis24-5.45)<delta && fabs(dis35-5.45)<delta)
    {
        s='H'; //helix                        
        return s;
    }

    delta=1.42;
    if (fabs(dis15-13  )<delta && fabs(dis14-10.4)<delta &&
        fabs(dis25-10.4)<delta && fabs(dis13-6.1 )<delta &&
        fabs(dis24-6.1 )<delta && fabs(dis35-6.1 )<delta)
    {
        s='E'; //strand
        return s;
    }

    if (dis15 < 8) s='T'; //turn
    return s;
}


/****************************************************************************
 * FUNCTIONS TO PORT TO PYTHON
 * these are still important to the the ports
 ****************************************************************************/

/*
 * port the USalign function make_sec(); specifically the _three_ argument 
 * version, which is specific to protein 2ndary structure calculations.
 *
 * !!! switched to using standard library <string> objects rather than char* 
 * b/c python porting did not like to compile char* to char...
 * 
 * return the string associated with the 2ndary structure prediction.
 * how does this handle unresolved or broken backbone?
*/ 
std::string make_sec(py::array_t<double> coords, 
		     int len)
{
    // fill a USalign-like array, xa, from the input array, coords
    double **xa;
    NewArray(&xa, len, 3);
   
    // get the buffer regions for the input array object
    py::buffer_info coords_info = coords.request();
    // create array filled with the pointers for the array elements
    auto inp_ptr  = static_cast <double *>(coords_info.ptr);
   
    // loop over input array dimensions and assign values to the USalign-like
    // array object
    int ndims = len * 3; // number of elements in the coords array
    int i; // index from ndims
    int j; // row index
    int k; // column index
    for (i = 0; i < ndims; i++) 
    {
        j = i / 3;
	k = i % 3;
	// fill USalign-like array, xa
	xa[j][k] = inp_ptr[i];
    }
    
    std::string sec(len, 'C');
    int j1, j2, j3, j4, j5;
    double d13, d14, d15, d24, d25, d35;
    for(i=0; i<len; i++)
    {
        j1=i-2;
        j2=i-1;
        j3=i;
        j4=i+1;
        j5=i+2;        
        
        if(j1>=0 && j5<len)
        {
            d13=sqrt(dist(xa[j1], xa[j3]));
            d14=sqrt(dist(xa[j1], xa[j4]));
            d15=sqrt(dist(xa[j1], xa[j5]));
            d24=sqrt(dist(xa[j2], xa[j4]));
            d25=sqrt(dist(xa[j2], xa[j5]));
            d35=sqrt(dist(xa[j3], xa[j5]));
            sec[i]=sec_str(d13, d14, d15, d24, d25, d35);            
        }    
    } 

    DeleteArray(&xa, len);

    return sec;
}


/*
 * port the USalign function getCloseK()
 * !!! only used if closeK_opt > 3
*/ 
py::array_t<double> getCloseK(py::array_t<double> coords, 
		              const int len, 
			      const int closeK_opt)
{
    // fill a USalign-like array, xa, from the input array, coords
    double **xa;
    NewArray(&xa, len, 3);
   
    // get the buffer regions for the input array object
    py::buffer_info coords_info = coords.request();
    // create array filled with the pointers for the array elements
    auto inp_ptr  = static_cast <double *>(coords_info.ptr);
   
    // loop over input array dimensions and assign values to the USalign-like
    // array object
    int ndims = len * 3; // used to loop over indices of the coords array
    int i; // index from ndims
    int j; // row index
    int k; // column index
    for (i = 0; i < ndims; i++) 
    {
        j = i / 3;
	k = i % 3;
	// fill USalign-like array, xa
	xa[j][k] = inp_ptr[i];
    }

    // a temp array to be filled with the dist**2 values between atoms in the
    // structure
    double **score;
    NewArray(&score, len+1, len+1);

    // calc dist**2_{i,j} "score" array
    for (i=0; i<len; i++)
    {
        // fill diag elements of score array with zeroes
	score[i+1][i+1]=0;
	// fill off-diag elements with the dist**2 to atom i
        for (j=i+1; j<len; j++) 
	{
	    score[j+1][i+1] = score[i+1][j+1] = dist(xa[i], xa[j]);
	}
    }

    // create the USalign-like k_nearest array to return
    // each row is filled with the cartesian coordinates for one of the 
    // closeK_opt neighbors for a residue in the structure
    // !!! maybe this array object doesn't need to be created
    double **k_nearest;
    NewArray(&k_nearest, len*closeK_opt, 3);
    // create a fillable, sortable vector object to collect the distances and 
    // associated residue indices. these then get sorted to find the neighbors.
    //auto pair = std::make_pair(0.0,0);
    std::vector< std::pair< double, int > > close_idx_vec; //(len, pair);
    // loop over each residue to find its closest neighbors
    for (i=0; i<len; i++)
    {
        // loop over residues, filling the close_idx_vec with distance and res
	// idx pairs
	for (j=0; j<len; j++)
        {
            close_idx_vec[j].first=score[i+1][j+1];
            close_idx_vec[j].second=j;
        }
        // sort the vector by the dist**2 values
	sort(close_idx_vec.begin(), close_idx_vec.end());
	// loop over the k-nearest neighbor integer
        for (k=0; k<closeK_opt; k++)
        {
            // get the neighbor's res index
	    j=close_idx_vec[k % len].second;
            // assign the neighbor's cartesian coords to the k_nearest array
	    // !!! for this port, maybe just go straight into the py::array_t 
	    // object
	    k_nearest[i*closeK_opt+k][0]=xa[j][0];
            k_nearest[i*closeK_opt+k][1]=xa[j][1];
            k_nearest[i*closeK_opt+k][2]=xa[j][2];
        }
    }

    // define the result array object to be filled
    auto result = py::array_t<double>(len*closeK_opt*3);
    // get the buffer regions for the array object
    py::buffer_info res_info = result.request();
    // create array filled with the pointers for the array elements
    auto out_ptr = static_cast <double *>(res_info.ptr);
    
    // loop over dimensions of the k_nearest 2d array and update the results 
    // array values
    for (i = 0; i < len*closeK_opt; ++i) 
    {
        for (int j = 0; j < 3; ++j)
	{
	    k = i*3 + j;
	    out_ptr[k] = k_nearest[i][j];
	}
    }

    // clean up 
    
    //close_idx_vec.clear();
    std::vector< std::pair< double, int> >().swap(close_idx_vec);
    
    DeleteArray(&xa, len);
    DeleteArray(&score, len+1);
    // maybe this array object isn't even necessary in the first place
    DeleteArray(&k_nearest, len*closeK_opt);

    return result;
}


/*
 * port the USalign function assign_sec_bond()
 * !!! only used for sNS alignment (mm_opt==6)
 *
*/ 
py::array_t<int> assign_sec_bond(const std::string sec, const int len)
{
    //declare the USalign-like array of shape (len x 2)
    int **sec_bond;
    NewArray(&sec_bond, len, 2);
    
    int i, j, k;
    int starti=-1;
    int endi=-1;
    char ss;
    char prev_ss=0;
    
    // loop over all residues
    for (i=0; i<len; i++)
    {
        // get the residue's 2ndary structure character
	ss=sec[i];
        // fill the sec_bond[i] vectors with -1s...
	sec_bond[i][0]=sec_bond[i][1]=-1;
	// check if 2ndary structure changes between prev and current and the
	// change isn't from turn to coil or vice-versa
        if (ss!=prev_ss && !(ss=='C' && prev_ss=='T') 
                        && !(ss=='T' && prev_ss=='C'))
        {
            // if this isn't the first residue in the 2ndary structure range
	    if (starti>=0) // previous SSE end
            {
                // assign residue i as the end
		endi=i;
                // loop over residues in the 2ndary structure range and fill
		// the residue's elements in sec_bond with the start and end 
		// indices
		for (j=starti;j<endi;j++)
                {
                    sec_bond[j][0]=starti;
                    sec_bond[j][1]=endi;
                }
            }
	    // check to see if the residue is a helix, sheet, or some RNA 
	    // strand type (ignore the RNA stuff for now)
            if (ss=='H' || ss=='E' || ss=='<' || ss=='>') starti=i;
            // otherwise, continue as if not part of a 2ndary structure
	    else starti=-1;
        }
        // update prev_ss
	prev_ss=sec[i];
    }
    // if we hit the end of the structure and starti still == -1, then we need
    // to include this 2ndary structure feature
    if (starti>=0) // previous SSE end
    {
        endi=i;
        for (j=starti;j<endi;j++)
        {
            sec_bond[j][0]=starti;
            sec_bond[j][1]=endi;
        }
    }
    // checking for single-residue 2ndary structure features; ignoring them if
    // found
    for (i=0;i<len;i++) 
    {
	if (sec_bond[i][1]-sec_bond[i][0]==1)
	{
	    sec_bond[i][0]=sec_bond[i][1]=-1;
	}
    }

    // sec_bond is now filled with (len x 2) ints associated with start,stop 
    // residues; convert it to a py::array_t<int>
    
    // define the result array object to be filled
    auto result = py::array_t<int>(len*2);
    // get the buffer regions for the array object
    py::buffer_info res_info = result.request();
    // create array filled with the pointers for the array elements
    auto out_ptr = static_cast <double *>(res_info.ptr);

    // loop over dimensions of the k_nearest 2d array and update the results 
    // array values
    for (i = 0; i < len; ++i) 
    {
        for (int j = 0; j < 2; ++j)
	{
	    k = i*2 + j;
	    out_ptr[k] = sec_bond[i][j];
	}
    }

    // clean up
    DeleteArray(&sec_bond, len);

    return result;
}


/*******************************************************************************
 * define the classes for user-input and parameter inputs for SOIalign_main()
 ******************************************************************************/

/*
 * define the input class to hold data associated with each structure, will be
 * passed as arguments into the wrapper function for the SOIalign_main function
 * important variables for structures: 
 *     coordinates, 
 *     sequence,
 *     2ndary structures string,
 *     nResidues length,
 *     k-nearest neighbors' coordinates array,
 *     2ndary structure features boundary array
*/

//class to handle the input parameters that are structure-specific
class alnStruct {
    public:
	// flattened 2d array nAtoms x 3
	py::array_t<float> coords; 
	
	// flattened 2d k-nearest-neighbors array; expected shape is 
	// (len * closeK_opt * 3); only created if closeK_opt >= 3
	py::array_t<float> k_nearest;
	
	// flattened 2d 2ndary struct boundaries array; expected shape is
	// (len x 2); only needed if mm_opt == 6 (sNS method)
	py::array_t<int> sec_bond; 
	std::string seq; // string, length len; holds the sequence
	std::string sec; // string, length len; holds the 2ndary struct labels
	int len;  // number of atoms to be aligned from structure
};


// class to handle the input parameters that are not structure-specific;
// sets default values for _all_ parameters that USalign passes :(
class alnParameters {
    public:
        // default set to -1 and then to 0 or 5 depending on mm_opt
	int closeK_opt = 0;
        int molec_types = -2; // -2 corresponds to a protein-protein alignment
	int mm_opt; // SOI alignment switch; 5 for fNS, 6 for sNS
	//int i_opt = 0; // > 0 if wish to apply an initial alignment
	int a_opt = 0; // > 0 if wish normalize TMscore by avg len of structs
	bool u_opt = false; // true if wish normalize TMscore by Lnorm_ass value
	double Lnorm_ass = 0; // only used if u_opt = true
	bool d_opt = false; // true if wish to use d0_scale for calc of pair weights
	double d0_scale = 0; // only used if d_opt
	bool fast_opt = false; // true if wish for fast but inaccurate alnment
};


/*******************************************************************************
 * define the class for SOIalign output data
 ******************************************************************************/

class outputResults {
    public:
        py::array_t<float> translation_vector; // 1d array vector for translation
        py::array_t<float> rotation_matrix;    // flattened 3x3 matrix for rotation 
        std::string seqM;	// mapping string between aligned sequences
        std::string seqA_mobile;// aligned seq; ordered by alnment to target sequence
        std::string seqA_target;// aligned seq
        double TM1, TM2; // normed by mobile and target lens, respectively
	double TM3, d0a; // only used if a_opt
       	double TM4, d0u; // only used if u_opt
	double TM5, d0_scale; // only used if d_opt
	double rmsd0; // final rmsd from Kabsch algo
        double d0_out = 5.0; // norm d0 value; also signify mapping in sequence alignment
	double Liden = 0;    // number of mapped residues w/ identical res type
        int n_ali8; // number of residues w/in d0_out
	int L_ali;   // n residues aligned, can include d > d0 pairs...
	double d0_mobile, d0_target; // d0 values calc'd when norming by 1 struct
	//double TM_ali, rmsd_ali; // looks like i_opt uses these as output vars
	//double TM_0, d0_0; // not used in the original output_results() func
};













/*******************************************************************************
 * TO REMIND MYSELF:
 * ! -> denotes an important func argument to be INPUT
 * | -> denotes an important func argument to be OUTPUT
 * all other arguments will be ignored for these python bindings
 *
 *USalign.cpp main() calls: 
 *    SOIalign(xname, yname, fname_super, fname_lign, fname_matrix, sequence, 
 *             Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt, u_opt, d_opt, 
 *             TMcut, infmt1_opt, infmt2_opt, ter_opt, split_opt, outfmt_opt, 
 *             fast_opt, cp_opt, mirror_opt, het_opt, atom_opt, autojustify, 
 *             mol_opt, dir_opt, dirpair_opt, dir1_opt, dir2_opt, chain2parse1, 
 *             chain2parse2, model2parse1, model2parse2, chain1_list, 
 *             chain2_list, se_opt, closeK_opt ! , mm_opt ! );
 *
 *
 *
 *Then, in USalign.cpp SOIalign(), SOIalign_main() is called:
 *    SOIalign_main(xa ! , ya ! , xk ! , yk ! , closeK_opt ! , seqx ! , seqy ! ,
 *    		    secx ! , secy ! , t0 | , u0 | , TM1 | , TM2 | , TM3, TM4, 
 *    		    TM5, d0_0 | , TM_0 | , d0A | , d0B | , d0u, d0a, d0_out, 
 *    		    seqM | , seqxA | , seqyA | , invmap | , rmsd0 | , L_ali | , 
 *    		    Liden | , 
 *    		    TM_ali, rmsd_ali, n_ali, n_ali8, xlen, ylen, sequence, 
 *                  Lnorm_ass, d0_scale, i_opt, a_opt, u_opt, d_opt, 
 *                  force_fast_opt, mol_vec1[chain_i]+mol_vec2[chain_j], 
 *                  dist_list, secx_bond, secy_bond, mm_opt);
 *
 *
 *
 *Results are output in the USalign.cpp SOIalign() call via a 
 *output_results() call:
 *    output_results(xname.substr(dir1_opt.size()+dir_opt.size()+dirpair_opt.size()),
 *                   yname.substr(dir2_opt.size()+dir_opt.size()+dirpair_opt.size()),
 *                   chainID_list1[chain_i], chainID_list2[chain_j], xlen, ylen,
 *                   t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out, 
 *                   seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden, n_ali8,
 *                   L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, d0B, Lnorm_ass, 
 *                   d0_scale, d0a, d0u, (m_opt?fname_matrix:"").c_str(),
 *                   outfmt_opt, ter_opt, false, split_opt, o_opt, fname_super,
 *                   i_opt, a_opt, u_opt, d_opt, mirror_opt, resi_vec1, resi_vec2); 
 *
 * only need t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out, 
 *           seqM, seqxA, seqyA, 
 *           Liden, n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, 
 *           d0B, d0a, d0u 
 *
 *
 ******************************************************************************/

/*
 *
 * *** =  important to return from python bindings; check on these in the 
 *        SOIalign() function.
 *
 * Other parameters used in SOIalign_main() function:
 * t0, u0, ***
 * TM1, TM2, *** 
 * TM3, TM4, TM5,
 * d0_0, ***
 * TM_0,
 * d0A, d0B, d0u, d0a, d0_out,
 * seqM, seqxA, seqyA, invmap, ***
 * rmsd0, L_ali, Liden, ***
 * TM_ali, rmsd_ali, n_ali, n_ali8,
 * sequence, Lnorm_ass, 
 * d0_scale,
 * i_opt, a_opt, u_opt, d_opt, fast_opt, 
 * mol_type, 
 * dist_list ***
 *
 * THESE PARAMETERS ARE PASSED TO SOIalign_main FOR SCOPE PURPOSES AND/OR FOR 
 * MEMORY CONTROL. 
 *
 * some of these parameters are hard-coded in the parent SOIalign() function in 
 * USalign.cpp. So I think these should be hardcoded in an input parameter 
 * class or in the wrapper function.
 */

///*******************************************************************************
// * define the wrapper function that reads input, runs SOIalign_main, and 
// * returns the output data.
// ******************************************************************************/
//
//OutputResults runSOIalign_main( alnstruct& mobile_data, 
//				alnstruct& target_data, 
//				alnParameters& parameters)
//{
//    // create variables to be used/filled w/in the SOIalign_main() call
//    double t0[3], u0[3][3];
//    double TM1, TM2;
//    double TM3, TM4, TM5; 
//    double d0_0, TM_0;
//    double d0A, d0B, d0u, d0a;
//    double d0_out=5.0;
//    string seqM, seqxA, seqyA;// for output alignment
//    double rmsd0 = 0.0;
//    int L_ali;                // Aligned length in standard_TMscore
//    double Liden=0;
//    double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
//    int n_ali=0;
//    int n_ali8=0;
//    
//    
//    // hardcode some default parameters to avoid their implementation...
//    // !!! design choice to still input these variables into the 
//    //     SOIalign_main() call.
//    int i_opt = 0;
//    vector<string> sequence;
//
//    // convert input arrays to their USalign-like versions...
//    int i,j,k;
//
//    // mobile coords
//    double **mobile_coords;
//    NewArray(&mobile_coords, mobile_data.len, 3);
//    // get the buffer regions for the input array object
//    py::buffer_info mobile_coords_info = mobile_data.coords.request();
//    // create array filled with the pointers for the array elements
//    auto mobile_coords_ptr  = static_cast <float *>(mobile_coords_info.ptr);
//    // py::array_t is flattened array, port it to the 2d USalign array
//    for (int i = 0; i < mobile_data.len *  3; i ++)
//    {
//	j = i / 3;
//	k = i % 3;
//	mobile_coords[j][k] = mobile_coords_ptr[i];
//    }
//   
//    // mobile k_nearest array
//    int **mobile_sec_bond;
//    NewArray(&mobile_sec_bond, mobile_data.len, 2);
//    // get the buffer regions for the input array object
//    py::buffer_info mobile_sec_bond_info = mobile_data.sec_bond.request();
//    // create array filled with the pointers for the array elements
//    auto mobile_sec_bond_ptr  = static_cast <int *>(mobile_sec_bond_info.ptr);
//    // py::array_t is flattened array, port it to the 2d USalign array
//    for (int i = 0; i < mobile_data.len * 2; i ++)
//    {
//	j = i / 3;
//	k = i % 2;
//	mobile_sec_bond[j][k] = mobile_sec_bond_ptr[i];
//    }
//  
//    // mobile k_nearest array
//    double **mobile_k_nearest;
//    NewArray(&mobile_k_nearest, mobile_data.len * parameters.closeK_opt, 3);
//    // get the buffer regions for the input array object
//    py::buffer_info mobile_k_nearest_info = mobile_data.k_nearest.request();
//    // create array filled with the pointers for the array elements
//    auto mobile_k_nearest_ptr  = static_cast <float *>(mobile_k_nearest_info.ptr);
//    // py::array_t is flattened array, port it to the 2d USalign array
//    for (int i = 0; i < mobile_data.len * parameters.closeK_opt *  3; i ++)
//    {
//	j = i / 3;
//	k = i % 3;
//	mobile_k_nearest[j][k] = mobile_k_nearest_ptr[i];
//    }
//  
//    // convert std::string to c-like char pointers
//    char *mobile_seq = mobile_data.seq.c_str();
//    char *mobile_sec = mobile_data.sec.c_str();
//
//    
//    // target coords
//    double **target_coords;
//    NewArray(&target_coords, target_data.len, 3);
//    // get the buffer regions for the input array object
//    py::buffer_info target_coords_info = target_data.coords.request();
//    // create array filled with the pointers for the array elements
//    auto target_coords_ptr  = static_cast <float *>(target_coords_info.ptr);
//    // py::array_t is flattened array, port it to the 2d USalign array
//    for (int i = 0; i < mobile_data.len *  3; i ++)
//    {
//	j = i / 3;
//	k = i % 3;
//	target_coords[j][k] = target_coords_ptr[i];
//    }
//   
//    // target k_nearest array
//    int **target_sec_bond;
//    NewArray(&target_sec_bond, target_data.len, 2);
//    // get the buffer regions for the input array object
//    py::buffer_info target_sec_bond_info = target_data.sec_bond.request();
//    // create array filled with the pointers for the array elements
//    auto target_sec_bond_ptr  = static_cast <int *>(target_sec_bond_info.ptr);
//    // py::array_t is flattened array, port it to the 2d USalign array
//    for (int i = 0; i < target_data.len * 2; i ++)
//    {
//	j = i / 3;
//	k = i % 2;
//	target_sec_bond[j][k] = target_sec_bond_ptr[i];
//    }
//  
//    // mobile k_nearest array
//    double **target_k_nearest;
//    NewArray(&target_k_nearest, target_data.len * parameters.closeK_opt, 3);
//    // get the buffer regions for the input array object
//    py::buffer_info target_k_nearest_info = target_data.k_nearest.request();
//    // create array filled with the pointers for the array elements
//    auto target_k_nearest_ptr  = static_cast <float *>(target_k_nearest_info.ptr);
//    // py::array_t is flattened array, port it to the 2d USalign array
//    for (int i = 0; i < target_data.len * parameters.closeK_opt *  3; i ++)
//    {
//	j = i / 3;
//	k = i % 3;
//	target_k_nearest[j][k] = target_k_nearest_ptr[i];
//    }
//  
//    // convert std::string to c-like char pointers
//    char *target_seq = target_data.seq.c_str();
//    char *target_sec = target_data.sec.c_str();
//
//
//    // prep other parameters
//
//    // if min of len values are too large (>1500), then do fast alignment 
//    // method...
//    bool force_fast_opt=(std::min(mobile_data.len,target_data.len)>1500)?true:fast_opt;
//    
//    // these lines aren't gonna work...
//    int *invmap = new int[target_data.len+1];
//    double *dist_list = new double[target_data.len+1];
//   
//    // run the SOIalign_main function as defined in the SOIalign.cpp file,
//    // for all its faults
//    SOIalign_main(mobile_coords, target_coords, 
//		  mobile_k_nearest, target_k_nearest,
//		  parameters.closeK_opt,
//		  mobile_seq, target_seq, mobile_sec, target_sec,
//		  t0, u0,
//		  TM1, TM2, TM3, TM4, TM5,
//		  d0_0, TM_0,
//		  d0A, d0B,
//		  d0u, d0a, d0_out,
//		  seqM, seqxA, seqyA,
//		  invmap,
//		  rmsd0,
//		  L_ali, Liden,
//		  TM_ali, rmsd_ali, n_ali,
//		  n_ali8,
//		  mobile_data.len, target_data.len, 
//		  sequence,
//		  parameters.Lnorm_ass, parameters.d0_scale,
//		  i_opt, 
//		  parameters.a_opt, parameters.u_opt, parameters.d_opt, 
//		  parameters.force_fast_opt,
//		  parameters.molec_types, dist_list,
//		  mobile_sec_bond, target_sec_bond, parameters.mm_opt);
//
//    // gather only the necessary output and return them in the output class as
//    // python accessible data types
//
//    OutputResults out; // instantiate the output object
//
//    // handling the translation and rotation arrays. flatten and return as 
//    // py::array_t
//    
//    // define the translation array object to be filled
//    auto trans_array = py::array_t<float>(3);
//    // get the buffer regions for the array object
//    py::buffer_info trans_info = trans_array.request();
//    // create array filled with the pointers for the array elements
//    auto trans_ptr = static_cast <float *>(trans_info.ptr);
//    // fill those elements 
//    for (int i = 0; i<3; i++)
//    {
//	trans_ptr[i] = t0[i];
//    }
//    out.translation_vector = trans_array;
//
//    // define the translation array object to be filled
//    auto rot_array = py::array_t<float>(9);
//    // get the buffer regions for the array object
//    py::buffer_info rot_info = rot_array.request();
//    // create array filled with the pointers for the array elements
//    auto rot_ptr = static_cast <float *>(rot_info.ptr);
//    // fill those elements 
//    for (int i = 0; i<3; i++)
//    {
//	for (int j = 0; j<3; j++)
//	{
//	    k = i*3 + j; 
//            rot_ptr[k] = u0[i][j];
//	}
//    }
//    out.rotation_vector = rot_array;
//
//    // seq results are char* so need to convert it to a std::string
//    out.seqM = str(seqM);
//    out.seqA_mobile = str(seqxA);
//    out.seqA_target = str(seqyA);
//
//    // handling tm score values
//    out.TM1 = TM1;
//    out.TM2 = TM2;
//    if (parameters.a_opt > 0) 
//    {
//	out.TM3 = TM3;
//	out.d0_a = d0a;
//    }
//    else 
//    {
//        out.TM3 = -1;
//        out.d0a = -1;
//    }
//
//    if (parameters.u_opt > 0) 
//    {
//	out.TM4 = TM4;
//        out.d0u = d0u;
//    }
//    else 
//    {
//	out.TM4 = -1;
//        out.d0u = -1;
//    }
//
//    if (parameters.d_opt > 0) 
//    {
//	out.TM5 = TM5;
//        out.d0_scale = parameters.d0_scale;
//    }
//    else 
//    {
//	out.TM5 = -1;
//        out.d0_scale = -1;
//    }
//
//    // handling other  value
//    out.rmsd0 = rmsd0;
//    out.d0_out = d0_out;
//    out.Liden = Liden;
//    out.n_ali8 = n_ali8;
//    out.L_ali = L_ali;
//    out.d0_mobile = d0A;
//    out.d0_target = d0B;
//
//    // cleaning up
//    DeleteArray(&mobile_coords, mobile_data.len);
//    DeleteArray(&target_coords, target_data.len);
//    DeleteArray(&mobile_sec_bond, mobile_data.len);
//    DeleteArray(&target_sec_bond, target_data.len);
//    DeleteArray(&mobile_k_nearest, mobile_data.len * parameters.closeK_opt);
//    DeleteArray(&target_k_nearest, target_data.len * parameters.closeK_opt);
//
//    return out;
//}

/*******************************************************************************
 * defining the pybind11 wrapper for SOIalign_main
 ******************************************************************************/

PYBIND11_MODULE(SOIalign_main, m) {
    m.doc() = "pybind11 port of SOIalign_main and related functions from USalign codebase"; 
    m.def("make_sec",
	  &make_sec,
	  "function to assign 2ndary structure character to each residue in the structure",
	  py::arg("coords"), 
	  py::arg("len"));
    
    m.def("getCloseK",
	  &getCloseK,
	  "function to determine nearest neighbors for each residue",
	  py::arg("coords"), 
	  py::arg("len"),
	  py::arg("closeK_opt"));

    m.def("assign_sec_bond",
	  &assign_sec_bond,
	  "function to identify the boundaries of large helix/sheet 2ndary structure elements",
	  py::arg("sec"), 
	  py::arg("len"));

    py::class_<alnStruct>(m, "alnStruct")
	.def(py::init<py::array_t<float>, 
		      py::array_t<float>, 
		      py::array_t<int>, 
		      std::string, 
		      std::string, 
		      int>())
	.def_readwrite("coords", &alnStruct::coords)
	.def_readwrite("k_nearest", &alnStruct::k_nearest)
	.def_readwrite("sec_bond", &alnStruct::sec_bond)
	.def_readwrite("seq", &alnStruct::seq)
	.def_readwrite("sec", &alnStruct::sec)
	.def_readwrite("len", &alnStruct::len);

    py::class_<alnParameters>(m, "alnParameters")
	.def(py::init<int, 
		      int, 
		      int, 
		      int, 
		      bool, 
		      double, 
		      bool, 
		      double, 
		      bool>())
	.def_readwrite("closeK_opt", &alnParameters::closeK_opt)
	.def_readwrite("molec_types", &alnParameters::molec_types)
	.def_readwrite("mm_opt", &alnParameters::mm_opt)
	.def_readwrite("a_opt", &alnParameters::a_opt)
	.def_readwrite("u_opt", &alnParameters::u_opt)
	.def_readwrite("Lnorm_ass", &alnParameters::Lnorm_ass)
	.def_readwrite("d_opt", &alnParameters::d_opt)
	.def_readwrite("d0_scale", &alnParameters::d0_scale)
	.def_readwrite("fast_opt", &alnParameters::fast_opt);

    py::class_<outputResults>(m, "outputResults")
	.def(py::init<py::array_t<float>,
		      py::array_t<float>,
		      std::string,
		      std::string,
		      std::string,
		      double,
		      double,
		      double,
		      double,
		      double,
		      double,
		      double,
		      double,
		      double,
		      double,
		      double,
		      int,
		      int,
		      double,
		      double>())
	.def_readwrite("translation_vector", &outputResults::translation_vector)
	.def_readwrite("rotation_matrix", &outputResults::rotation_matrix)
	.def_readwrite("seqM", &outputResults::seqM)
	.def_readwrite("seqA_mobile", &outputResults::seqA_mobile)
	.def_readwrite("seqA_target", &outputResults::seqA_target)
	.def_readwrite("TM1", &outputResults::TM1)
	.def_readwrite("TM2", &outputResults::TM2)
	.def_readwrite("TM3", &outputResults::TM3)
	.def_readwrite("d0a", &outputResults::d0a)
	.def_readwrite("TM4", &outputResults::TM4)
	.def_readwrite("d0u", &outputResults::d0u)
	.def_readwrite("TM5", &outputResults::TM5)
	.def_readwrite("d0_scale", &outputResults::d0_scale)
	.def_readwrite("rmsd0", &outputResults::rmsd0)
	.def_readwrite("d0_out", &outputResults::d0_out)
	.def_readwrite("Liden", &outputResults::Liden)
	.def_readwrite("n_ali8", &outputResults::n_ali8)
	.def_readwrite("L_ali", &outputResults::L_ali)
	.def_readwrite("d0_mobile", &outputResults::d0_mobile)
	.def_readwrite("d0_target", &outputResults::d0_target);
    
//    m.def("SOIalign_main", // define the function name
//          &SOIalign_main,  // set pointer to hold this function
//	  "function to calculate the Sequence-Order-Independent alignment between two sets of coordinates", // blurb about what the function does
//	  py::arg("mobile_coords") = , // maps to xa variable
//	  py::arg("target_coords") = , // maps to ya variable
//	  py::arg("mobile_neighbor_list") = ,	// maps to xk variable
//	  py::arg("target_neighbor_list") = ,	// maps to yk variable
//	  py::arg("n_nearest_neighbors") = 3,	// maps to closeK_opt variable
//	  py::arg("mobile_seq") = , // maps to seqx variable
//	  py::arg("target_seq") = , // maps to seqy variable
//	  py::arg("mobile_sec") = , // maps to secx variable
//	  py::arg("target_sec") = , // maps to secy variable
//	  // variables to be changed internal to the function... why are we passing these into the code???
//	  py::arg("translation_vector") = , // maps to t0 variable
//	  py::arg("rotation_matrix") = , // maps to u0 variable
//	  py::arg("TM1") = , // maps to TM1 variable
//	  py::arg("TM2") = , // maps to TM1 variable
//	  py::arg("d0_0") = , // maps to d0_0 variable
//	  py::arg("TM_0") = , // maps to TM_0 variable
//	  py::arg("d0A") = , // maps to d0A variable
//	  py::arg("d0B") = , // maps to d0A variable
//	  py::arg("aligned_seq") = , // maps to seqM variable
//	  py::arg("mobile_aligned_seq") = , // maps to seqxA variable
//	  py::arg("target_aligned_seq") = , // maps to seqyA variable
//	  py::arg("inverse_map") = , // maps to invmap variable
//	  py::arg("final_rmds") = , // maps to rmsd0 variable
//	  py::arg("n_aligned_residues") = , // maps to n_ali8 variable
//	  py::arg("mobile_len") = , // maps to xlen variable
//	  py::arg("target_len") = , // maps to ylen variable
//	  py::arg("distance_list") = , // maps to dist_list variable
//	  py::arg("mobile_sec_bond") = , // maps to secx_bond variable
//	  py::arg("target_sec_bond") = , // maps to secy_bond variable
//	  py::arg("mm_opt") = , // maps to variable mm_opt
//	  //py::arg("") = , //
//	  py::return_value_policy::copy);
}
