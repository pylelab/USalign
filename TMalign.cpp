/*
===============================================================================
   This is a re-implementation of TM-align algorithm in C/C++. The code was 
   is written by Jianyi Yang and later updated by Jianjie Wu at The Yang Zhang 
   lab, Department of Computational Medicine and Bioinformatics, University of 
   Michigan, 100 Washtenaw Avenue, Ann Arbor, MI 48109-2218. Please report bugs 
   and questions to jianjiew@umich.edu or zhng@umich.edu

   DISCLAIMER:
     Permission to use, copy, modify, and distribute this program for 
     any purpose, with or without fee, is hereby granted, provided that
     the notices on the head, the reference information, and this
     copyright notice appear in all copies or substantial portions of 
     the Software. It is provided "as is" without express or implied 
     warranty.
   *************** updating history ********************************
   2012/01/24: A C/C++ code of TM-align was constructed by Jianyi Yang
   2016/05/21: Several updates of this program were made by Jianjie Wu, including
              (1) fixed several compiling bugs
	      (2) made I/O of C/C++ version consistent with the Fortran version
	      (3) added outputs including full-atom and ligand structures
	      (4) added options of '-i', '-I' and '-m'
   2016/05/25: fixed a bug on PDB file reading
===============================================================================
*/
#include "basic_define.h"
#include "global_var.h"
#include "param_set.h"

using namespace std;


#include "basic_fun.h"
#include "NW.h"
#include "Kabsch.h"
#include "TMalign.h"

void print_extra_help()
{
	cout <<
"Additional options: \n"
"    -fast    Fast but slightly inaccurate alignment\n"
"\n"
//"    -dir1    Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
//"             under 'chain1_folder'\n"
//"             $ TMalign -dir1 chain1_folder/ chain1_list chain2\n"
//"\n"
//"    -ter     Strings to mark the end of a chain\n"
//"             3: (default) 'TER', 'END', or different chain ID\n"
//"             2: 'END', or different chain ID\n"
//"             1: 'END'\n"
    <<endl;
}

void print_help(bool h_opt=false)
{
	cout <<
"\n"
"*****************************************************************************\n"
"* TM-align (Version "<< version
                      <<    "): An algorithm for protein structure alignment *\n"
"* Based on statistics:                                                      *\n"
"*          0.0 < TM-score < 0.30, random structural similarity              *\n"
"*          0.5 < TM-score < 1.00, in about the same fold                    *\n"
"* Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)       *\n"
"* Please email your comments and suggestions to Yang Zhang (zhng@umich.edu) *\n"
"*****************************************************************************\n"
"\n"
"Usage: TMalign PDB1.pdb PDB2.pdb [Options]\n\n"
"Options:\n"
"    -u    TM-score normalized by user assigned length\n"
"          warning: it should be >= minimum length of the two structures\n"
"          otherwise, TM-score may be >1\n\n"
"    -a    TM-score normalized by the average length of two structures\n"
"          T or F, (default F)\n\n"
"    -i    Ask TM-align to start with an alignment, specified in fasta\n"
"          file 'align.txt'\n\n"
"    -I    Ask TM-align to stick the alignment to 'align.txt'\n\n"
"    -m    Output TM-align rotation matrix:\n\n"
"    -d    TM-score scaled by an assigned d0, e.g. 5 Angstroms\n\n"
"    -o    output the superposition to TM.sup, TM.sup_all and TM.sup_atm\n"
"          $ TMalign chain1 chain2 -o TM.sup\n"
"          To view superimposed C-alpha traces of aligned regions by rasmol:\n"
"          $ rasmol -script TM.sup\n"
"          To view superimposed C-alpha traces of all regions:\n"
"          $ rasmol -script TM.sup_all\n\n"
"          To view superimposed full-atom structures of aligned regions:\n"
"          $ rasmol -script TM.sup_atm\n\n"
"          To view superimposed full-atom structures of all regions:\n"
"          $ rasmol -script TM.sup_all_atm\n\n"
"          To view superimposed full-atom structures of all regions with ligands:\n"
"          $ rasmol -script TM.sup_all_atm_lig\n\n"
"    -v    print the version of TM-align\n\n"
"    -h    print this help\n\n"
"    (Options -u, -a, -d -o won't change the final structure alignment)\n\n"
"Example usages:\n"
"    TMalign PDB1.pdb PDB2.pdb\n"
"    TMalign PDB1.pdb PDB2.pdb -u 100 -d 5.0\n"
"    TMalign PDB1.pdb PDB2.pdb -a T -o PDB1.sup\n"
"    TMalign PDB1.pdb PDB2.pdb -i align.txt\n"
"    TMalign PDB1.pdb PDB2.pdb -m matrix.txt\n"
    <<endl;

    if (h_opt) print_extra_help();
       
    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();
	
	
	clock_t t1, t2;
	t1 = clock();

    /**********************/
    /*    get argument    */ 
    /**********************/
    char xname[MAXLEN], yname[MAXLEN],  Lnorm_ave[MAXLEN];
	bool A_opt, B_opt, h_opt=false;
	A_opt = B_opt = o_opt = a_opt = u_opt = d_opt = v_opt = false;
	i_opt = false;// set -i flag to be false
	m_opt = false;// set -m flag to be false
	char fname_lign[MAXLEN] = "";
	char fname_matrix[MAXLEN] = "";// set names to ""
	I_opt = false;// set -I flag to be false
    fast_opt = false;// set -fast flag to be false

	int nameIdx = 0;
	for(int i = 1; i < argc; i++)
	{
		if ( !strcmp(argv[i],"-o") && i < (argc-1) ) 
		{
			strcpy(out_reg, argv[i + 1]);      o_opt = true; i++;
		}
		else if ( !strcmp(argv[i],"-u") && i < (argc-1) ) 
		{
			Lnorm_ass = atof(argv[i + 1]); u_opt = true; i++;
		}
		else if ( !strcmp(argv[i],"-a") && i < (argc-1) ) 
		{
			strcpy(Lnorm_ave, argv[i + 1]);     a_opt = true; i++;
		}
		else if ( !strcmp(argv[i],"-d") && i < (argc-1) ) 
		{
			d0_scale = atof(argv[i + 1]); d_opt = true; i++;
		}
		else if ( !strcmp(argv[i],"-v") ) 
		{
			v_opt = true; 
		}
		else if ( !strcmp(argv[i],"-h") ) 
		{ 
			h_opt = true; 
		}
		else if ( !strcmp(argv[i],"-i") && i < (argc-1) ) 
		{
			strcpy(fname_lign, argv[i + 1]);      i_opt = true; i++;
		}
		else if (!strcmp(argv[i], "-m") && i < (argc-1) ) 
		{
			strcpy(fname_matrix, argv[i + 1]);      m_opt = true; i++;
		}// get filename for rotation matrix
		else if (!strcmp(argv[i], "-I") && i < (argc-1) ) 
		{
			strcpy(fname_lign, argv[i + 1]);      I_opt = true; i++;
		}
		else if (!strcmp(argv[i], "-fast"))
		{
			fast_opt = true;
		}
		else
		{
			if (nameIdx == 0)
			{
				strcpy(xname, argv[i]);      B_opt = true;
			}
			else if (nameIdx == 1)
			{
				strcpy(yname, argv[i]);      A_opt = true;
			}
			nameIdx++;
		}
	}

	if(!B_opt || !A_opt)
	{

		if( h_opt ) print_help(h_opt);
		
		if(v_opt)
		{
			cout<<endl<<
"*****************************************************************************\n"
"* TM-align (Version "<< version <<"): A protein structural alignment algorithm *\n"
"* Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)       *\n"
"* Please email your comments and suggestions to Yang Zhang (zhng@umich.edu) *\n"
"*****************************************************************************\n"
                <<endl;
			exit(EXIT_FAILURE);
		}
	}

	if( !B_opt )
	{
		cout << "Please provide structure B" << endl;
		exit(EXIT_FAILURE);
	}		
	if( !A_opt )
	{
		cout << "Please provide structure A" << endl;
		exit(EXIT_FAILURE);
	}


	if( a_opt )
	{
		if(!strcmp(Lnorm_ave, "T"))
		{
		}
		else if(!strcmp(Lnorm_ave, "F"))
		{
			a_opt=false;
		}
		else
		{
			cout << "Wrong value for option -a!  It should be T or F" << endl;
			exit(EXIT_FAILURE);
		}
	}
	if( u_opt )
	{
		if(Lnorm_ass<=0)
		{
			cout << "Wrong value for option -u!  It should be >0" << endl;
			exit(EXIT_FAILURE);
		}
	}
	if( d_opt )
	{
		if(d0_scale<=0)
		{
			cout << "Wrong value for option -d!  It should be >0" << endl;
			exit(EXIT_FAILURE);
		}
	}
	////// read initial alignment file from 'alignment.txt' //////
	string basename = string(argv[0]);
	int idx = basename.find_last_of("\\");
	basename = basename.substr(0, idx + 1);
	if (i_opt || I_opt)// Ask TM-align to start with an alignment,specified in fasta file 'align.txt'
	{
		if (fname_lign == "")
		{
			cout << "Please provide a file name for option -i!" << endl;
			exit(EXIT_FAILURE);
		}
		// open alignment file
		int n_p = 0;// number of structures in alignment file
		string line;

		string fullpath = basename + fname_lign;
		char path1[1000];
		strcpy(path1, fullpath.c_str());
		ifstream fileIn(path1);
		if (fileIn.is_open())
		{
			bool bContinue = true;
			while (fileIn.good() && bContinue)
			{
				getline(fileIn, line);
				if (line.compare(0, 1, ">") == 0)// Flag for a new structure
				{
					strcpy(sequence[n_p], "");
					n_p++;
					if (n_p > 2)
						bContinue = false;
				}
				else// Read data
				{
					if (n_p > 0 && line!="")
					{
						strcat(sequence[n_p-1], line.c_str());
					}
				}
			}
			fileIn.close();
		}
		else
		{
			cout << "\nAlignment file does not exist.\n";
			exit(EXIT_FAILURE);
		}
	
		if (n_p < 2)
		{
			cout << "\nERROR: Fasta format is wrong, two proteins should be included.\n";
			exit(EXIT_FAILURE);
		}
		else
		{
			if (strlen(sequence[0]) != strlen(sequence[1]))
			{
				cout << "\nWarning: FASTA format may be wrong, the length in alignment should be equal respectively to the aligned proteins.\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	if (m_opt)// Output TM - align rotation matrix: matrix.txt
	{
		if (fname_matrix == "")
		{
			cout << "Please provide a file name for option -m!" << endl;
			exit(EXIT_FAILURE);
		}
	}


    /*******************/
    /*    load data    */ 
    /*******************/
    load_PDB_allocate_memory(xname, yname);


    /***********************/
    /*    parameter set    */ 
    /***********************/
	parameter_set4search(xlen, ylen);          //please set parameters in the function
    int simplify_step     = 40;               //for similified search engine
    int score_sum_method  = 8;                //for scoring method, whether only sum over pairs with dis<score_d8
        
	int i;
    int *invmap0          = new int[ylen+1]; 
    int *invmap           = new int[ylen+1]; 
    double TM, TMmax=-1;
	for(i=0; i<ylen; i++)
	{
		invmap0[i]=-1;
	}	


	double ddcc=0.4;
	if(Lnorm <= 40) ddcc=0.1;   //Lnorm was setted in parameter_set4search
	double local_d0_search = d0_search;

	//************************************************//
	//    get initial alignment from user's input:    //
	//    Stick to the initial alignment              //
	//************************************************//
	char dest[1000];
	bool bAlignStick = false;
	if (I_opt)// if input has set parameter for "-I"
	{
		for (int j = 1; j < ylen; j++)// Set aligned position to be "-1"
			invmap[j] = -1;

		int i1 = -1;// in C version, index starts from zero, not from one
		int i2 = -1;
		int L1 = strlen(sequence[0]);
		int L2 = strlen(sequence[1]);
		int L = min(L1, L2);// Get positions for aligned residues
		for (int kk1 = 0; kk1 < L; kk1++)
		{
			if (sequence[0][kk1] != '-')
				i1++;
			if (sequence[1][kk1] != '-')
			{
				i2++;
				if (i2 >= ylen || i1 >= xlen)
					kk1 = L;
				else
				{
					if (sequence[0][kk1] != '-')
					{
						invmap[i2] = i1;
					}
				}
			}
		}
		
		//--------------- 2. Align proteins from original alignment
		double prevD0_MIN = D0_MIN;// stored for later use
		int prevLnorm = Lnorm;
		double prevd0 = d0;
		TM_ali = standard_TMscore(xa, ya, xlen, ylen, invmap, L_ali, rmsd_ali);
		D0_MIN = prevD0_MIN;
		Lnorm = prevLnorm;
		d0 = prevd0;
		TM = detailed_search_standard(xa, ya, xlen, ylen, invmap, t, u, 40, 8, local_d0_search, true);
		if (TM > TMmax)
		{
			TMmax = TM;
			for (i = 0; i<ylen; i++)
				invmap0[i] = invmap[i];
		}
		bAlignStick = true;
	}

    /******************************************************/
    /*    get initial alignment with gapless threading    */ 
    /******************************************************/
	if (!bAlignStick)
	{
		get_initial(xa, ya, xlen, ylen, invmap0);
		//find the max TMscore for this initial alignment with the simplified search_engin
		TM = detailed_search(xa, ya, xlen, ylen, invmap0, t, u, simplify_step, score_sum_method, local_d0_search);
		if (TM>TMmax)
		{
			TMmax = TM;
		}
		//run dynamic programing iteratively to find the best alignment
		TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30, local_d0_search);
		if (TM>TMmax)
		{
			TMmax = TM;
			for (int i = 0; i<ylen; i++)
			{
				invmap0[i] = invmap[i];
			}
		}


		/************************************************************/
		/*    get initial alignment based on secondary structure    */
		/************************************************************/
		get_initial_ss(xa, ya, xlen, ylen, invmap);
		TM = detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method, local_d0_search);
		if (TM>TMmax)
		{
			TMmax = TM;
			for (int i = 0; i<ylen; i++)
			{
				invmap0[i] = invmap[i];
			}
		}
		if (TM > TMmax*0.2)
		{
			TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30, local_d0_search);
			if (TM>TMmax)
			{
				TMmax = TM;
				for (int i = 0; i<ylen; i++)
				{
					invmap0[i] = invmap[i];
				}
			}
		}


		/************************************************************/
		/*    get initial alignment based on local superposition    */
		/************************************************************/
		//=initial5 in original TM-align
		if (get_initial5(xa, ya, xlen, ylen, invmap))
		{
			TM = detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method, local_d0_search);
			if (TM>TMmax)
			{
				TMmax = TM;
				for (int i = 0; i<ylen; i++)
				{
					invmap0[i] = invmap[i];
				}
			}
			if (TM > TMmax*ddcc)
			{
				TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 2, local_d0_search);
				if (TM>TMmax)
				{
					TMmax = TM;
					for (int i = 0; i<ylen; i++)
					{
						invmap0[i] = invmap[i];
					}
				}
			}
		}
		else
		{
			cout << endl << endl << "Warning: initial alignment from local superposition fail!" << endl << endl << endl;
		}


		/*******************************************************************************/
		/*    get initial alignment based on previous alignment+secondary structure    */
		/*******************************************************************************/
		//=initial3 in original TM-align
		get_initial_ssplus(xa, ya, xlen, ylen, invmap0, invmap);
		TM = detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method, local_d0_search);
		if (TM>TMmax)
		{
			TMmax = TM;
			for (i = 0; i<ylen; i++)
			{
				invmap0[i] = invmap[i];
			}
		}
		if (TM > TMmax*ddcc)
		{
			TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30, local_d0_search);
			if (TM>TMmax)
			{
				TMmax = TM;
				for (i = 0; i<ylen; i++)
				{
					invmap0[i] = invmap[i];
				}
			}
		}


		/*******************************************************************/
		/*    get initial alignment based on fragment gapless threading    */
		/*******************************************************************/
		//=initial4 in original TM-align
		get_initial_fgt(xa, ya, xlen, ylen, xresno, yresno, invmap);
		TM = detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method, local_d0_search);
		if (TM>TMmax)
		{
			TMmax = TM;
			for (i = 0; i<ylen; i++)
			{
				invmap0[i] = invmap[i];
			}
		}
		if (TM > TMmax*ddcc)
		{
			TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 1, 2, 2, local_d0_search);
			if (TM>TMmax)
			{
				TMmax = TM;
				for (i = 0; i<ylen; i++)
				{
					invmap0[i] = invmap[i];
				}
			}
		}


		//************************************************//
		//    get initial alignment from user's input:    //
		//************************************************//
		if (i_opt)// if input has set parameter for "-i"
		{
			for (int j = 0; j < ylen; j++)// Set aligned position to be "-1"
				invmap[j] = -1;
		
			int i1 = -1;// in C version, index starts from zero, not from one
			int i2 = -1;
			int L1 = strlen(sequence[0]);
			int L2 = strlen(sequence[1]);
			int L = min(L1, L2);// Get positions for aligned residues
			for (int kk1 = 0; kk1 < L; kk1++)
			{
				if (sequence[0][kk1] != '-')
					i1++;
				if (sequence[1][kk1] != '-')
				{
					i2++;
					if (i2 >= ylen || i1 >= xlen)
						kk1 = L;
					else
					{
						if (sequence[0][kk1] != '-')
						{
							invmap[i2] = i1;
						}
					}
				}
			}

			//--------------- 2. Align proteins from original alignment
			double prevD0_MIN = D0_MIN;// stored for later use
			int prevLnorm = Lnorm;
			double prevd0 = d0;
			TM_ali = standard_TMscore(xa, ya, xlen, ylen, invmap, L_ali, rmsd_ali);
			D0_MIN = prevD0_MIN;
			Lnorm = prevLnorm;
			d0 = prevd0;

			TM = detailed_search_standard(xa, ya, xlen, ylen, invmap, t, u, 40, 8, local_d0_search, true);
			if (TM > TMmax)
			{
				TMmax = TM;
				for (i = 0; i<ylen; i++)
					invmap0[i] = invmap[i];
			}
			TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, (fast_opt)?2:30, local_d0_search);// Different from get_initial, get_initial_ss and get_initial_ssplus
			if (TM>TMmax)
			{
				TMmax = TM;
				for (i = 0; i<ylen; i++)
				{
					invmap0[i] = invmap[i];
				}
			}
		}
	}


	
    //*******************************************************************//
    //    The alignment will not be changed any more in the following    //
    //*******************************************************************//
	//check if the initial alignment is generated approately	
	bool flag=false;
	for(i=0; i<ylen; i++)
	{
		if(invmap0[i]>=0)
		{
			flag=true;
			break;			
		}			
	}		
	if(!flag) 
	{
		cout << "There is no alignment between the two proteins!" << endl;
		cout << "Program stop with no result!" << endl;
		return 1;
	}


    //********************************************************************//
    //    Detailed TMscore search engine --> prepare for final TMscore    //
    //********************************************************************//       
    //run detailed TMscore search engine for the best alignment, and 
	//extract the best rotation matrix (t, u) for the best alginment
	simplify_step=1;
    if (fast_opt) simplify_step=40;
    score_sum_method=8;
	TM = detailed_search_standard(xa, ya, xlen, ylen, invmap0, t, u, simplify_step, score_sum_method, local_d0_search, false);
		

	//select pairs with dis<d8 for final TMscore computation and output alignment
	int n_ali8, k=0;
	int n_ali=0;
	int *m1, *m2;
	double d;
	m1=new int[xlen]; //alignd index in x
	m2=new int[ylen]; //alignd index in y
	do_rotation(xa, xt, xlen, t, u);
	k=0;
    for(int j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)//aligned
        {
			n_ali++;        
            d=sqrt(dist(&xt[i][0], &ya[j][0]));
			if (d <= score_d8 || (I_opt == true))
			{
				m1[k]=i;
				m2[k]=j;

				xtm[k][0]=xa[i][0];
                xtm[k][1]=xa[i][1];
                xtm[k][2]=xa[i][2];
                    
                ytm[k][0]=ya[j][0];
                ytm[k][1]=ya[j][1];
                ytm[k][2]=ya[j][2];			

				r1[k][0] = xt[i][0];
				r1[k][1] = xt[i][1];
				r1[k][2] = xt[i][2];
				r2[k][0] = ya[j][0];
				r2[k][1] = ya[j][1];
				r2[k][2] = ya[j][2];
				
				k++;
			}
		}
	}
	n_ali8=k;

	double rmsd0 = 0.0;
	Kabsch(r1, r2, n_ali8, 0, &rmsd0, t, u);// rmsd0 is used for final output, only recalculate rmsd0, not t & u
	rmsd0 = sqrt(rmsd0 / n_ali8);


    //****************************************//
    //              Final TMscore             //
    //    Please set parameters for output    //
    //****************************************//
	double rmsd;
	double t0[3], u0[3][3];
	double TM1, TM2;
	double d0_out=5.0;  
    simplify_step=1;
    score_sum_method=0;

	double d0_0, TM_0;
	double Lnorm_0=ylen;
	
	
	//normalized by length of structure A
	parameter_set4final(Lnorm_0);
	d0A=d0;
	d0_0=d0A;
	local_d0_search = d0_search;
	TM1 = TMscore8_search(xtm, ytm, n_ali8, t0, u0, simplify_step, score_sum_method, &rmsd, local_d0_search);
	TM_0 = TM1;

	//normalized by length of structure B
	parameter_set4final(xlen+0.0);
	d0B=d0;
	local_d0_search = d0_search;
	TM2 = TMscore8_search(xtm, ytm, n_ali8, t, u, simplify_step, score_sum_method, &rmsd, local_d0_search);


	if(a_opt)
	{
		//normalized by average length of structures A, B
		Lnorm_0=(xlen+ylen)*0.5;
		parameter_set4final(Lnorm_0);
		d0a=d0;
		d0_0=d0a;
		local_d0_search = d0_search;

		TM3 = TMscore8_search(xtm, ytm, n_ali8, t0, u0, simplify_step, score_sum_method, &rmsd, local_d0_search);
		TM_0=TM3;
	}
	if(u_opt)
	{	
		//normalized by user assigned length		
		parameter_set4final(Lnorm_ass);		
		d0u=d0;		
		d0_0=d0u;
		Lnorm_0=Lnorm_ass;
		local_d0_search = d0_search;
		TM4 = TMscore8_search(xtm, ytm, n_ali8, t0, u0, simplify_step, score_sum_method, &rmsd, local_d0_search);
		TM_0=TM4;
	}
	if(d_opt)
	{	
		//scaled by user assigned d0
		parameter_set4scale(ylen, d0_scale);
		d0_out=d0_scale;
		d0_0=d0_scale;
		//Lnorm_0=ylen;
		Lnorm_d0=Lnorm_0;
		local_d0_search = d0_search;
		TM5 = TMscore8_search(xtm, ytm, n_ali8, t0, u0, simplify_step, score_sum_method, &rmsd, local_d0_search);
		TM_0=TM5;
	}
        
	output_results(xname, yname, xlen, ylen, t0, u0, TM1, TM2, rmsd0, d0_out, m1, m2, n_ali8, n_ali, TM_0, Lnorm_0, d0_0, fname_matrix);


    //*************************//
    //    Done! Free memory    //
    //*************************//
	free_memory();

    delete [] invmap0;
    delete [] invmap;
	delete [] m1;
	delete [] m2;

    t2 = clock();    
    float diff = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    printf("\nTotal running time is %5.2f seconds\n", diff);        
 	return 0;	
}
