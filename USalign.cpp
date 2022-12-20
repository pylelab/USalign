/* command line argument parsing and document of US-align main program */

#include "MMalign.h"
#include "SOIalign.h"
#include "flexalign.h"

using namespace std;

void print_version()
{
    cout << 
"\n"
" ********************************************************************\n"
" * US-align (Version 20220924)                                      *\n"
" * Universal Structure Alignment of Proteins and Nucleic Acids      *\n"
" * Reference: C Zhang, M Shine, AM Pyle, Y Zhang. (2022) Nat Methods*\n"
" * Please email comments and suggestions to zhang@zhanggroup.org    *\n"
" ********************************************************************"
    << endl;
}

void print_extra_help()
{
    cout <<
"Additional options:\n"
"      -v  Print the version of US-align\n"
"\n"
"      -a  TM-score normalized by the average length of two structures\n"
"          T or F, (default F). -a does not change the final alignment.\n"
"\n"
"   -fast  Fast but slightly inaccurate alignment\n"
"\n"
"    -dir  Perform all-against-all alignment among the list of PDB\n"
"          chains listed by 'chain_list' under 'chain_folder'. Note\n"
"          that the slash is necessary.\n"
"          $ USalign -dir chain_folder/ chain_list\n"
"\n"
"   -dir1  Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
"          under 'chain1_folder'. Note that the slash is necessary.\n"
"          $ USalign -dir1 chain1_folder/ chain1_list chain2\n"
"\n"
"   -dir2  Use chain1 to search a list of PDB chains listed by 'chain2_list'\n"
"          under 'chain2_folder'\n"
"          $ USalign chain1 -dir2 chain2_folder/ chain2_list\n"
"\n"
" -suffix  (Only when -dir1 and/or -dir2 are set, default is empty)\n"
"          add file name suffix to files listed by chain1_list or chain2_list\n"
"\n"
"   -atom  4-character atom name used to represent a residue.\n"
"          Default is \" C3'\" for RNA/DNA and \" CA \" for proteins\n"
"          (note the spaces before and after CA).\n"
"\n"
"  -split  Whether to split PDB file into multiple chains\n"
"           0: treat the whole structure as one single chain\n"
"           1: treat each MODEL as a separate chain\n"
"           2: (default) treat each chain as a separate chain\n"
"\n"
" -outfmt  Output format\n"
"           0: (default) full output\n"
"           1: fasta format compact output\n"
"           2: tabular format very compact output\n"
"          -1: full output, but without version or citation information\n"
"\n"
"  -TMcut  -1: (default) do not consider TMcut\n"
"          Values in [0.5,1): Do not proceed with TM-align for this\n"
"          structure pair if TM-score is unlikely to reach TMcut.\n"
"          TMcut is normalized as set by -a option:\n"
"          -2: normalized by longer structure length\n"
"          -1: normalized by shorter structure length\n"
"           0: (default, same as F) normalized by second structure\n"
"           1: same as T, normalized by average structure length\n"
"\n"
" -mirror  Whether to align the mirror image of input structure\n"
"           0: (default) do not align mirrored structure\n"
"           1: align mirror of Structure_1 to origin Structure_2,\n"
"              which usually requires the '-het 1' option:\n"
"              $ USalign 4glu.pdb 3p9w.pdb -mirror 1 -het 1\n"
"\n"
"    -het  Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"           0: (default) only align 'ATOM  ' residues\n"
"           1: align both 'ATOM  ' and 'HETATM' residues\n"
"           2: align both 'ATOM  ' and MSE residues\n"
"\n"
"   -full  Whether to show full pairwise alignment of individual chains for\n"
"          -mm 2 or 4. T or F, (default F)\n"
//"\n"
//" -closeK  Number of closest atoms used for sequence order independent\n"
//"          initial alignment. default: 5\n"
//"\n"
//" -hinge   Maximum number of hinge allowed in flexible alignment. default: 9\n"
"\n"
"   -se    Do not perform superposition. Useful for extracting alignment from\n"
"          superposed structure pairs\n"
"\n"
" -infmt1  Input format for structure_11\n"
" -infmt2  Input format for structure_2\n"
"          -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
"           0: PDB format\n"
"           1: SPICKER format\n"
//"           2: xyz format\n"
"           3: PDBx/mmCIF format\n"
"\n"
"Advanced usage 1 (generate an image for a pair of superposed structures):\n"
"    USalign 1cpc.pdb 1mba.pdb -o sup\n"
"    pymol -c -d @sup_all_atm.pml -g sup_all_atm.png\n"
"\n"
"Advanced usage 2 (a quick search of query.pdb against I-TASSER PDB library):\n"
"    wget https://zhanggroup.org/library/PDB.tar.bz2\n"
"    tar -xjvf PDB.tar.bz2\n"
"    USalign query.pdb -dir2 PDB/ PDB/list -suffix .pdb -outfmt 2 -fast\n"
    <<endl;
}

void print_help(bool h_opt=false)
{
    print_version();
    cout <<
"\n"
"Usage: USalign PDB1.pdb PDB2.pdb [Options]\n"
"\n"
"Options:\n"
"    -mol  Type of molecule(s) to align.\n"
"          auto: (default) align both protein and nucleic acids.\n"
"          prot: only align proteins in a structure.\n"
"          RNA : only align RNA and DNA in a structure.\n"
"\n"
"     -mm  Multimeric alignment option:\n"
"          0: (default) alignment of two monomeric structures\n"
"          1: alignment of two multi-chain oligomeric structures\n"
"          2: alignment of individual chains to an oligomeric structure\n"
"             $ USalign -dir1 monomers/ list oligomer.pdb -ter 0 -mm 2\n"
"          3: alignment of circularly permuted structure\n"
"          4: alignment of multiple monomeric chains into a consensus alignment\n"
"             $ USalign -dir chains/ list -suffix .pdb -mm 4\n"
"          5: fully non-sequential (fNS) alignment\n"
"          6: semi-non-sequential (sNS) alignment\n"
"          To use -mm 1 or -mm 2, '-ter' option must be 0 or 1.\n"
"\n"
"    -ter  Number of chains to align.\n"
"          3: only align the first chain, or the first segment of the\n"
"             first chain as marked by the 'TER' string in PDB file\n"
"          2: (default) only align the first chain\n"
"          1: align all chains of the first model (recommended for aligning\n"
"             asymmetric units)\n"
"          0: align all chains from all models (recommended for aligning\n"
"             biological assemblies, i.e. biounits)\n"
"\n"
" -TMscore Whether to perform TM-score superposition without structure-based\n"
"          alignment. The same as -byresi.\n"
"          0: (default) sequence independent structure alignment\n"
"          1: superpose two structures by assuming that a pair of residues\n"
"             with the same residue index are equivalent between the two\n"
"             structures\n"
"          2: superpose two complex structures, assuming that a pair of\n"
"             residues with the same residue index and the same chain ID\n"
"             are equivalent between the two structures\n"
//"          3: (similar to TMscore '-c' option; used with -ter 0 or 1)\n"
//"             align by residue index and order of chain\n"
//"          4: sequence dependent alignment: perform Needleman-Wunsch\n"
//"             global sequence alignment, followed by TM-score superposition\n"
"          5: sequence dependent alignment: perform glocal sequence\n"
"             alignment followed by TM-score superposition.\n"
"             -byresi 5 is the same as -seq\n"
"          6: superpose two complex structures by first deriving optimal\n"
"             chain mapping, followed by TM-score superposition for residues\n"
"             with the same residue ID\n"
"\n"
"      -I  Use the final alignment specified by FASTA file 'align.txt'\n"
"\n"
"      -i  Use alignment specified by 'align.txt' as an initial alignment\n"
"\n"
"      -m  Output rotation matrix for superposition\n"
"\n"
"      -d  TM-score scaled by an assigned d0, e.g., '-d 3.5' reports MaxSub\n"
"          score, where d0 is 3.5 Angstrom. -d does not change final alignment.\n"
"\n"
"      -u  TM-score normalized by an assigned length. It should be >= length\n"
"          of protein to avoid TM-score >1. -u does not change final alignment.\n"
"\n"
"      -o  Output superposed structure1 to sup.* for PyMOL viewing.\n"
"          $ USalign structure1.pdb structure2.pdb -o sup\n"
"          $ pymol -d @sup.pml                # C-alpha trace aligned region\n"
"          $ pymol -d @sup_all.pml            # C-alpha trace whole chain\n"
"          $ pymol -d @sup_atm.pml            # full-atom aligned region\n"
"          $ pymol -d @sup_all_atm.pml        # full-atom whole chain\n"
"          $ pymol -d @sup_all_atm_lig.pml    # full-atom with all molecules\n"
"\n"
" -rasmol  Output superposed structure1 to sup.* for RasMol viewing.\n"
"          $ USalign structure1.pdb structure2.pdb -rasmol sup\n"
"          $ rasmol -script sup               # C-alpha trace aligned region\n"
"          $ rasmol -script sup_all           # C-alpha trace whole chain\n"
"          $ rasmol -script sup_atm           # full-atom aligned region\n"
"          $ rasmol -script sup_all_atm       # full-atom whole chain\n"
"          $ rasmol -script sup_all_atm_lig   # full-atom with all molecules\n"
"\n"
//"      -h  Print the full help message, including additional options\n"
//"\n"
"Example usages ('gunzip' program is needed to read .gz compressed files):\n"
"    USalign 101m.cif.gz 1mba.pdb             # pairwise monomeric protein alignment\n"
"    USalign 1qf6.cif 5yyn.pdb.gz -mol RNA    # pairwise monomeric RNA alignment\n"
"    USalign model.pdb native.pdb -TMscore 1  # calculate TM-score between two conformations of a monomer\n"
"    USalign 4v4a.cif 4v49.cif -mm 1 -ter 1   # oligomeric alignment for asymmetic units\n"
"    USalign 3ksc.pdb1 4lej.pdb1 -mm 1 -ter 0 # oligomeric alignment for biological units\n"
"    USalign 1ajk.pdb.gz 2ayh.pdb.gz -mm 3    # circular permutation alignment\n"
    <<endl;

    //if (h_opt) 
        print_extra_help();

    exit(EXIT_SUCCESS);
}

/* TMalign, RNAalign, CPalign, TMscore */
int TMalign(string &xname, string &yname, const string &fname_super,
    const string &fname_lign, const string &fname_matrix,
    vector<string> &sequence, const double Lnorm_ass, const double d0_scale,
    const bool m_opt, const int  i_opt, const int o_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const double TMcut,
    const int infmt1_opt, const int infmt2_opt, const int ter_opt,
    const int split_opt, const int outfmt_opt, const bool fast_opt,
    const int cp_opt, const int mirror_opt, const int het_opt,
    const string &atom_opt, const string &mol_opt, const string &dir_opt,
    const string &dir1_opt, const string &dir2_opt, const int byresi_opt,
    const vector<string> &chain1_list, const vector<string> &chain2_list,
    const bool se_opt)
{
    /* declare previously global variables */
    vector<vector<string> >PDB_lines1; // text of chain1
    vector<vector<string> >PDB_lines2; // text of chain2
    vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    vector<string> chainID_list1;      // list of chainID1
    vector<string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    r;                  // residue index
    int    xlen, ylen;         // chain length
    int    xchainnum,ychainnum;// number of chains in a PDB file
    char   *seqx, *seqy;       // for the protein sequence 
    char   *secx, *secy;       // for the secondary structure 
    double **xa, **ya;         // for input vectors xa[0...xlen-1][0..2] and
                               // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure 
                               // --> superpose xa onto ya
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2
    int read_resi=byresi_opt;  // whether to read residue index
    if (byresi_opt==0 && o_opt) read_resi=2;

    /* loop over file names */
    for (i=0;i<chain1_list.size();i++)
    {
        /* parse chain 1 */
        xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1,
            mol_vec1, ter_opt, infmt1_opt, atom_opt, split_opt, het_opt);
        if (!xchainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=PDB_lines1[chain_i].size();
            if (mol_opt=="RNA") mol_vec1[chain_i]=1;
            else if (mol_opt=="protein") mol_vec1[chain_i]=-1;
            if (!xlen)
            {
                cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<endl;
                continue;
            }
            else if (xlen<3)
            {
                cerr<<"Sequence is too short <3!: "<<xname<<endl;
                continue;
            }
            NewArray(&xa, xlen, 3);
            seqx = new char[xlen + 1];
            secx = new char[xlen + 1];
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx, 
                resi_vec1, read_resi);
            if (mirror_opt) for (r=0;r<xlen;r++) xa[r][2]=-xa[r][2];
            if (mol_vec1[chain_i]>0) make_sec(seqx,xa, xlen, secx,atom_opt);
            else make_sec(xa, xlen, secx); // secondary structure assignment

            for (j=(dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                /* parse chain 2 */
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2, chainID_list2,
                        mol_vec2, ter_opt, infmt2_opt, atom_opt, split_opt,
                        het_opt);
                    if (!ychainnum)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain number 0."<<endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ylen=PDB_lines2[chain_j].size();
                    if (mol_opt=="RNA") mol_vec2[chain_j]=1;
                    else if (mol_opt=="protein") mol_vec2[chain_j]=-1;
                    if (!ylen)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain length 0."<<endl;
                        continue;
                    }
                    else if (ylen<3)
                    {
                        cerr<<"Sequence is too short <3!: "<<yname<<endl;
                        continue;
                    }
                    NewArray(&ya, ylen, 3);
                    seqy = new char[ylen + 1];
                    secy = new char[ylen + 1];
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy,
                        resi_vec2, read_resi);
                    if (mol_vec2[chain_j]>0)
                         make_sec(seqy, ya, ylen, secy, atom_opt);
                    else make_sec(ya, ylen, secy);

                    if (byresi_opt) extract_aln_from_resi(sequence,
                        seqx,seqy,resi_vec1,resi_vec2,byresi_opt);

                    /* declare variable specific to this pair of TMalign */
                    double t0[3], u0[3][3];
                    double TM1, TM2;
                    double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                    double d0_0, TM_0;
                    double d0A, d0B, d0u, d0a;
                    double d0_out=5.0;
                    string seqM, seqxA, seqyA;// for output alignment
                    double rmsd0 = 0.0;
                    int L_ali;                // Aligned length in standard_TMscore
                    double Liden=0;
                    double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                    int n_ali=0;
                    int n_ali8=0;
                    bool force_fast_opt=(getmin(xlen,ylen)>1500)?true:fast_opt;

                    /* entry function for structure alignment */
                    if (cp_opt) CPalign_main(
                        xa, ya, seqx, seqy, secx, secy,
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                        mol_vec1[chain_i]+mol_vec2[chain_j],TMcut);
                    else if (se_opt)
                    {
                        int *invmap = new int[ylen+1];
                        u0[0][0]=u0[1][1]=u0[2][2]=1;
                        u0[0][1]=         u0[0][2]=
                        u0[1][0]=         u0[1][2]=
                        u0[2][0]=         u0[2][1]=
                        t0[0]   =t0[1]   =t0[2]   =0;
                        se_main(
                            xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                            seqM, seqxA, seqyA,
                            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                            xlen, ylen, sequence, Lnorm_ass, d0_scale,
                            i_opt, a_opt, u_opt, d_opt,
                            mol_vec1[chain_i]+mol_vec2[chain_j], 
                            outfmt_opt, invmap);
                        if (outfmt_opt>=2) 
                        {
                            Liden=L_ali=0;
                            int r1,r2;
                            for (r2=0;r2<ylen;r2++)
                            {
                                r1=invmap[r2];
                                if (r1<0) continue;
                                L_ali+=1;
                                Liden+=(seqx[r1]==seqy[r2]);
                            }
                        }
                        delete [] invmap;
                    }
                    else TMalign_main(
                        xa, ya, seqx, seqy, secx, secy,
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                        mol_vec1[chain_i]+mol_vec2[chain_j],TMcut);

                    /* print result */
                    if (outfmt_opt==0) print_version();
                    int left_num=0;
                    int right_num=0;
                    int left_aln_num=0;
                    int right_aln_num=0;
                    bool after_cp=false;
                    if (cp_opt) after_cp=output_cp(
                        xname.substr(dir1_opt.size()+dir_opt.size()),
                        yname.substr(dir2_opt.size()+dir_opt.size()),
                        seqxA,seqyA,outfmt_opt,left_num,right_num,
                        left_aln_num,right_aln_num);
                    output_results(
                        xname.substr(dir1_opt.size()+dir_opt.size()),
                        yname.substr(dir2_opt.size()+dir_opt.size()),
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5,
                        rmsd0, d0_out, seqM.c_str(),
                        seqxA.c_str(), seqyA.c_str(), Liden,
                        n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
                        d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u, 
                        (m_opt?fname_matrix:"").c_str(),
                        outfmt_opt, ter_opt, false, split_opt, o_opt,
                        fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
                        resi_vec1, resi_vec2);
                    if (cp_opt && outfmt_opt<=0)
                    {
                        cout<<"###############\t###############\n"
                            <<"#Aligned atom 1\tAligned atom 2#\n";
                        size_t r1=right_num;
                        size_t r2=0;
                        size_t r;
                        for (r=0;r<seqxA.size();r++)
                        {
                            r1+=seqxA[r]!='-';
                            r2+=seqyA[r]!='-';
                            if (seqxA[r]=='*')
                            {
                                cout<<"###### Circular\tPermutation ###\n";
                                r1=0;
                            }
                            else if (seqxA[r]!='-' && seqyA[r]!='-')
                            {
                                cout<<PDB_lines1[chain_i][r1-1].substr(12,15)<<'\t'
                                    <<PDB_lines2[chain_j][r2-1].substr(12,15)<<'\n';
                            }
                        }
                        cout<<"###############\t###############"<<endl;
                    }

                    /* Done! Free memory */
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();
                    DeleteArray(&ya, ylen);
                    delete [] seqy;
                    delete [] secy;
                    resi_vec2.clear();
                } // chain_j
                if (chain2_list.size()>1)
                {
                    yname.clear();
                    for (chain_j=0;chain_j<ychainnum;chain_j++)
                        PDB_lines2[chain_j].clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
            PDB_lines1[chain_i].clear();
            DeleteArray(&xa, xlen);
            delete [] seqx;
            delete [] secx;
            resi_vec1.clear();
        } // chain_i
        xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
    if (chain2_list.size()==1)
    {
        yname.clear();
        for (chain_j=0;chain_j<ychainnum;chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        resi_vec2.clear();
        chainID_list2.clear();
        mol_vec2.clear();
    }
    return 0;
}

/* MMalign if more than two chains. TMalign if only one chain */
int MMalign(const string &xname, const string &yname,
    const string &fname_super, const string &fname_lign,
    const string &fname_matrix, vector<string> &sequence,
    const double d0_scale, const bool m_opt, const int o_opt,
    const int a_opt, const bool d_opt, const bool full_opt,
    const double TMcut, const int infmt1_opt, const int infmt2_opt,
    const int ter_opt, const int split_opt, const int outfmt_opt,
    bool fast_opt, const int mirror_opt, const int het_opt,
    const string &atom_opt, const string &mol_opt,
    const string &dir1_opt, const string &dir2_opt,
    const vector<string> &chain1_list, const vector<string> &chain2_list,
    const int byresi_opt)
{
    /* declare previously global variables */
    vector<vector<vector<double> > > xa_vec; // structure of complex1
    vector<vector<vector<double> > > ya_vec; // structure of complex2
    vector<vector<char> >seqx_vec; // sequence of complex1
    vector<vector<char> >seqy_vec; // sequence of complex2
    vector<vector<char> >secx_vec; // secondary structure of complex1
    vector<vector<char> >secy_vec; // secondary structure of complex2
    vector<int> mol_vec1;          // molecule type of complex1, RNA if >0
    vector<int> mol_vec2;          // molecule type of complex2, RNA if >0
    vector<string> chainID_list1;  // list of chainID1
    vector<string> chainID_list2;  // list of chainID2
    vector<int> xlen_vec;          // length of complex1
    vector<int> ylen_vec;          // length of complex2
    int    i,j;                    // chain index
    int    xlen, ylen;             // chain length
    double **xa, **ya;             // structure of single chain
    char   *seqx, *seqy;           // for the protein sequence 
    char   *secx, *secy;           // for the secondary structure 
    int    xlen_aa,ylen_aa;        // total length of protein
    int    xlen_na,ylen_na;        // total length of RNA/DNA
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2

    /* parse complex */
    parse_chain_list(chain1_list, xa_vec, seqx_vec, secx_vec, mol_vec1,
        xlen_vec, chainID_list1, ter_opt, split_opt, mol_opt, infmt1_opt,
        atom_opt, mirror_opt, het_opt, xlen_aa, xlen_na, o_opt, resi_vec1);
    if (xa_vec.size()==0) PrintErrorAndQuit("ERROR! 0 chain in complex 1");
    parse_chain_list(chain2_list, ya_vec, seqy_vec, secy_vec, mol_vec2,
        ylen_vec, chainID_list2, ter_opt, split_opt, mol_opt, infmt2_opt,
        atom_opt, 0, het_opt, ylen_aa, ylen_na, o_opt, resi_vec2);
    if (ya_vec.size()==0) PrintErrorAndQuit("ERROR! 0 chain in complex 2");
    int len_aa=getmin(xlen_aa,ylen_aa);
    int len_na=getmin(xlen_na,ylen_na);
    if (a_opt)
    {
        len_aa=(xlen_aa+ylen_aa)/2;
        len_na=(xlen_na+ylen_na)/2;
    }
    int i_opt=0;
    if (byresi_opt) i_opt=3;

    /* perform monomer alignment if there is only one chain */
    if (xa_vec.size()==1 && ya_vec.size()==1)
    {
        xlen = xlen_vec[0];
        ylen = ylen_vec[0];
        seqx = new char[xlen+1];
        seqy = new char[ylen+1];
        secx = new char[xlen+1];
        secy = new char[ylen+1];
        NewArray(&xa, xlen, 3);
        NewArray(&ya, ylen, 3);
        copy_chain_data(xa_vec[0],seqx_vec[0],secx_vec[0], xlen,xa,seqx,secx);
        copy_chain_data(ya_vec[0],seqy_vec[0],secy_vec[0], ylen,ya,seqy,secy);
        
        /* declare variable specific to this pair of TMalign */
        double t0[3], u0[3][3];
        double TM1, TM2;
        double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
        double d0_0, TM_0;
        double d0A, d0B, d0u, d0a;
        double d0_out=5.0;
        string seqM, seqxA, seqyA;// for output alignment
        double rmsd0 = 0.0;
        int L_ali;                // Aligned length in standard_TMscore
        double Liden=0;
        double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
        int n_ali=0;
        int n_ali8=0;
        
        if (byresi_opt) extract_aln_from_resi(sequence,
            seqx,seqy,resi_vec1,resi_vec2,byresi_opt);

        /* entry function for structure alignment */
        TMalign_main(xa, ya, seqx, seqy, secx, secy,
            t0, u0, TM1, TM2, TM3, TM4, TM5,
            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
            seqM, seqxA, seqyA,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, 0, d0_scale,
            i_opt, a_opt, false, d_opt, fast_opt,
            mol_vec1[0]+mol_vec2[0],TMcut);

        /* print result */
        output_results(
            xname.substr(dir1_opt.size()),
            yname.substr(dir2_opt.size()),
            chainID_list1[0], chainID_list2[0],
            xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
            seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden,
            n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, d0B,
            0, d0_scale, d0a, d0u, (m_opt?fname_matrix:"").c_str(),
            outfmt_opt, ter_opt, true, split_opt, o_opt, fname_super,
            0, a_opt, false, d_opt, mirror_opt, resi_vec1, resi_vec2);

        /* clean up */
        seqM.clear();
        seqxA.clear();
        seqyA.clear();
        delete[]seqx;
        delete[]seqy;
        delete[]secx;
        delete[]secy;
        DeleteArray(&xa,xlen);
        DeleteArray(&ya,ylen);

        vector<vector<vector<double> > >().swap(xa_vec); // structure of complex1
        vector<vector<vector<double> > >().swap(ya_vec); // structure of complex2
        vector<vector<char> >().swap(seqx_vec); // sequence of complex1
        vector<vector<char> >().swap(seqy_vec); // sequence of complex2
        vector<vector<char> >().swap(secx_vec); // secondary structure of complex1
        vector<vector<char> >().swap(secy_vec); // secondary structure of complex2
        mol_vec1.clear();       // molecule type of complex1, RNA if >0
        mol_vec2.clear();       // molecule type of complex2, RNA if >0
        chainID_list1.clear();  // list of chainID1
        chainID_list2.clear();  // list of chainID2
        xlen_vec.clear();       // length of complex1
        ylen_vec.clear();       // length of complex2
        return 0;
    }

    /* declare TM-score tables */
    int chain1_num=xa_vec.size();
    int chain2_num=ya_vec.size();
    vector<string> tmp_str_vec(chain2_num,"");
    double **TMave_mat;
    double **ut_mat; // rotation matrices for all-against-all alignment
    int ui,uj,ut_idx;
    NewArray(&TMave_mat,chain1_num,chain2_num);
    NewArray(&ut_mat,chain1_num*chain2_num,4*3);
    vector<vector<string> >seqxA_mat(chain1_num,tmp_str_vec);
    vector<vector<string> > seqM_mat(chain1_num,tmp_str_vec);
    vector<vector<string> >seqyA_mat(chain1_num,tmp_str_vec);

    double maxTMmono=-1;
    int maxTMmono_i,maxTMmono_j;

    /* get all-against-all alignment */
    if (len_aa+len_na>500) fast_opt=true;
    for (i=0;i<chain1_num;i++)
    {
        xlen=xlen_vec[i];
        if (xlen<3)
        {
            for (j=0;j<chain2_num;j++) TMave_mat[i][j]=-1;
            continue;
        }
        seqx = new char[xlen+1];
        secx = new char[xlen+1];
        NewArray(&xa, xlen, 3);
        copy_chain_data(xa_vec[i],seqx_vec[i],secx_vec[i],
            xlen,xa,seqx,secx);

        for (j=0;j<chain2_num;j++)
        {
            ut_idx=i*chain2_num+j;
            for (ui=0;ui<4;ui++)
                for (uj=0;uj<3;uj++) ut_mat[ut_idx][ui*3+uj]=0;
            ut_mat[ut_idx][0]=1;
            ut_mat[ut_idx][4]=1;
            ut_mat[ut_idx][8]=1;

            if (mol_vec1[i]*mol_vec2[j]<0) //no protein-RNA alignment
            {
                TMave_mat[i][j]=-1;
                continue;
            }

            ylen=ylen_vec[j];
            if (ylen<3)
            {
                TMave_mat[i][j]=-1;
                continue;
            }
            seqy = new char[ylen+1];
            secy = new char[ylen+1];
            NewArray(&ya, ylen, 3);
            copy_chain_data(ya_vec[j],seqy_vec[j],secy_vec[j],
                ylen,ya,seqy,secy);

            /* declare variable specific to this pair of TMalign */
            double t0[3], u0[3][3];
            double TM1, TM2;
            double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out=5.0;
            string seqM, seqxA, seqyA;// for output alignment
            double rmsd0 = 0.0;
            int L_ali;                // Aligned length in standard_TMscore
            double Liden=0;
            double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
            int n_ali=0;
            int n_ali8=0;

            int Lnorm_tmp=len_aa;
            if (mol_vec1[i]+mol_vec2[j]>0) Lnorm_tmp=len_na;
            
            if (byresi_opt)
            {
                int total_aln=extract_aln_from_resi(sequence,
                    seqx,seqy,resi_vec1,resi_vec2,xlen_vec,ylen_vec, i, j);
                seqxA_mat[i][j]=sequence[0];
                seqyA_mat[i][j]=sequence[1];
                if (total_aln>xlen+ylen-3)
                {
                    for (ui=0;ui<3;ui++) for (uj=0;uj<3;uj++) 
                        ut_mat[ut_idx][ui*3+uj]=(ui==uj)?1:0;
                    for (uj=0;uj<3;uj++) ut_mat[ut_idx][9+uj]=0;
                    TMave_mat[i][j]=0;
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();

                    delete[]seqy;
                    delete[]secy;
                    DeleteArray(&ya,ylen);
                    continue;
                }
            }

            /* entry function for structure alignment */
            TMalign_main(xa, ya, seqx, seqy, secx, secy,
                t0, u0, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                i_opt, false, true, false, fast_opt,
                mol_vec1[i]+mol_vec2[j],TMcut);

            /* store result */
            for (ui=0;ui<3;ui++)
                for (uj=0;uj<3;uj++) ut_mat[ut_idx][ui*3+uj]=u0[ui][uj];
            for (uj=0;uj<3;uj++) ut_mat[ut_idx][9+uj]=t0[uj];
            seqxA_mat[i][j]=seqxA;
            seqyA_mat[i][j]=seqyA;
            TMave_mat[i][j]=TM4*Lnorm_tmp;
            if (TMave_mat[i][j]>maxTMmono)
            {
                maxTMmono=TMave_mat[i][j];
                maxTMmono_i=i;
                maxTMmono_j=j;
            }

            /* clean up */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();

            delete[]seqy;
            delete[]secy;
            DeleteArray(&ya,ylen);
        }

        delete[]seqx;
        delete[]secx;
        DeleteArray(&xa,xlen);
    }

    /* calculate initial chain-chain assignment */
    int *assign1_list; // value is index of assigned chain2
    int *assign2_list; // value is index of assigned chain1
    assign1_list=new int[chain1_num];
    assign2_list=new int[chain2_num];
    double total_score=enhanced_greedy_search(TMave_mat, assign1_list,
        assign2_list, chain1_num, chain2_num);
    if (total_score<=0) PrintErrorAndQuit("ERROR! No assignable chain");

    /* refine alignment for large oligomers */
    int aln_chain_num=count_assign_pair(assign1_list,chain1_num);
    bool is_oligomer=(aln_chain_num>=3);
    if (aln_chain_num==2) // dimer alignment
    {
        int na_chain_num1,na_chain_num2,aa_chain_num1,aa_chain_num2;
        count_na_aa_chain_num(na_chain_num1,aa_chain_num1,mol_vec1);
        count_na_aa_chain_num(na_chain_num2,aa_chain_num2,mol_vec2);

        /* align protein-RNA hybrid dimer to another hybrid dimer */
        if (na_chain_num1==1 && na_chain_num2==1 && 
            aa_chain_num1==1 && aa_chain_num2==1) is_oligomer=false;
        /* align pure protein dimer or pure RNA dimer */
        else if ((getmin(na_chain_num1,na_chain_num2)==0 && 
                    aa_chain_num1==2 && aa_chain_num2==2) ||
                 (getmin(aa_chain_num1,aa_chain_num2)==0 && 
                    na_chain_num1==2 && na_chain_num2==2))
        {
            adjust_dimer_assignment(xa_vec,ya_vec,xlen_vec,ylen_vec,mol_vec1,
                mol_vec2,assign1_list,assign2_list,seqxA_mat,seqyA_mat);
            is_oligomer=false; // cannot refiner further
        }
        else is_oligomer=true; /* align oligomers to dimer */
    }

    if (aln_chain_num>=3 || is_oligomer) // oligomer alignment
    {
        /* extract centroid coordinates */
        double **xcentroids;
        double **ycentroids;
        NewArray(&xcentroids, chain1_num, 3);
        NewArray(&ycentroids, chain2_num, 3);
        double d0MM=getmin(
            calculate_centroids(xa_vec, chain1_num, xcentroids),
            calculate_centroids(ya_vec, chain2_num, ycentroids));

        /* refine enhanced greedy search with centroid superposition */
        //double het_deg=check_heterooligomer(TMave_mat, chain1_num, chain2_num);
        homo_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na, ut_mat);
        hetero_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na);
        
        /* clean up */
        DeleteArray(&xcentroids, chain1_num);
        DeleteArray(&ycentroids, chain2_num);
    }

    /* store initial assignment */
    int init_pair_num=count_assign_pair(assign1_list,chain1_num);
    int *assign1_init, *assign2_init;
    assign1_init=new int[chain1_num];
    assign2_init=new int[chain2_num];
    double **TMave_init;
    NewArray(&TMave_init,chain1_num,chain2_num);
    vector<vector<string> >seqxA_init(chain1_num,tmp_str_vec);
    vector<vector<string> >seqyA_init(chain1_num,tmp_str_vec);
    vector<string> sequence_init;
    copy_chain_assign_data(chain1_num, chain2_num, sequence_init,
        seqxA_mat,  seqyA_mat,  assign1_list, assign2_list, TMave_mat,
        seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init);

    /* perform iterative alignment */
    double max_total_score=0; // ignore old total_score because previous
                              // score was from monomeric chain superpositions
    int max_iter=5-(int)((len_aa+len_na)/200);
    if (max_iter<2) max_iter=2;
    if (byresi_opt==0) MMalign_iter(max_total_score, max_iter, xa_vec, ya_vec,
        seqx_vec, seqy_vec, secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec,
        ylen_vec, xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num,
        chain2_num, TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list,
        sequence, d0_scale, fast_opt);

    /* sometime MMalign_iter is even worse than monomer alignment */
    if (byresi_opt==0 && max_total_score<maxTMmono)
    {
        copy_chain_assign_data(chain1_num, chain2_num, sequence,
            seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
            seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
        for (i=0;i<chain1_num;i++)
        {
            if (i!=maxTMmono_i) assign1_list[i]=-1;
            else assign1_list[i]=maxTMmono_j;
        }
        for (j=0;j<chain2_num;j++)
        {
            if (j!=maxTMmono_j) assign2_list[j]=-1;
            else assign2_list[j]=maxTMmono_i;
        }
        sequence[0]=seqxA_mat[maxTMmono_i][maxTMmono_j];
        sequence[1]=seqyA_mat[maxTMmono_i][maxTMmono_j];
        max_total_score=maxTMmono;
        MMalign_iter(max_total_score, max_iter, xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
            xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num, chain2_num,
            TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list, sequence,
            d0_scale, fast_opt);
    }

    /* perform cross chain alignment
     * in some cases, this leads to dramatic improvement, esp for homodimer */
    int iter_pair_num=count_assign_pair(assign1_list,chain1_num);
    if (iter_pair_num>=init_pair_num) copy_chain_assign_data(
        chain1_num, chain2_num, sequence_init,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat,
        seqxA_init, seqyA_init, assign1_init,  assign2_init,  TMave_init);
    double max_total_score_cross=max_total_score;
    if (byresi_opt==0 && len_aa+len_na<10000)
    {
        MMalign_dimer(max_total_score_cross, xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
            xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num, chain2_num,
            TMave_init, seqxA_init, seqyA_init, assign1_init, assign2_init,
            sequence_init, d0_scale, fast_opt);
        if (max_total_score_cross>max_total_score) 
        {
            max_total_score=max_total_score_cross;
            copy_chain_assign_data(chain1_num, chain2_num, sequence,
                seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
                seqxA_mat,  seqyA_mat,  assign1_list, assign2_list, TMave_mat);
        }
    } 

    /* final alignment */
    if (outfmt_opt==0) print_version();
    MMalign_final(xname.substr(dir1_opt.size()), yname.substr(dir2_opt.size()),
        chainID_list1, chainID_list2,
        fname_super, fname_lign, fname_matrix,
        xa_vec, ya_vec, seqx_vec, seqy_vec,
        secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
        xa, ya, seqx, seqy, secx, secy, len_aa, len_na,
        chain1_num, chain2_num, TMave_mat,
        seqxA_mat, seqM_mat, seqyA_mat, assign1_list, assign2_list, sequence,
        d0_scale, m_opt, o_opt, outfmt_opt, ter_opt, split_opt,
        a_opt, d_opt, fast_opt, full_opt, mirror_opt, resi_vec1, resi_vec2);

    /* clean up everything */
    delete [] assign1_list;
    delete [] assign2_list;
    DeleteArray(&TMave_mat,chain1_num);
    DeleteArray(&ut_mat,   chain1_num*chain2_num);
    vector<vector<string> >().swap(seqxA_mat);
    vector<vector<string> >().swap(seqM_mat);
    vector<vector<string> >().swap(seqyA_mat);
    vector<string>().swap(tmp_str_vec);

    delete [] assign1_init;
    delete [] assign2_init;
    DeleteArray(&TMave_init,chain1_num);
    vector<vector<string> >().swap(seqxA_init);
    vector<vector<string> >().swap(seqyA_init);

    vector<vector<vector<double> > >().swap(xa_vec); // structure of complex1
    vector<vector<vector<double> > >().swap(ya_vec); // structure of complex2
    vector<vector<char> >().swap(seqx_vec); // sequence of complex1
    vector<vector<char> >().swap(seqy_vec); // sequence of complex2
    vector<vector<char> >().swap(secx_vec); // secondary structure of complex1
    vector<vector<char> >().swap(secy_vec); // secondary structure of complex2
    mol_vec1.clear();       // molecule type of complex1, RNA if >0
    mol_vec2.clear();       // molecule type of complex2, RNA if >0
    vector<string>().swap(chainID_list1);  // list of chainID1
    vector<string>().swap(chainID_list2);  // list of chainID2
    xlen_vec.clear();       // length of complex1
    ylen_vec.clear();       // length of complex2
    vector<string> ().swap(resi_vec1);  // residue index for chain1
    vector<string> ().swap(resi_vec2);  // residue index for chain2
    return 1;
}


/* alignment individual chains to a complex. */
int MMdock(const string &xname, const string &yname, const string &fname_super, 
    const string &fname_matrix, vector<string> &sequence, const double Lnorm_ass,
    const double d0_scale, const bool m_opt, const int o_opt,
    const int a_opt, const bool u_opt, const bool d_opt,
    const double TMcut, const int infmt1_opt, const int infmt2_opt,
    const int ter_opt, const int split_opt, const int outfmt_opt,
    bool fast_opt, const int mirror_opt, const int het_opt,
    const string &atom_opt, const string &mol_opt,
    const string &dir1_opt, const string &dir2_opt,
    const vector<string> &chain1_list, const vector<string> &chain2_list)
{
    /* declare previously global variables */
    vector<vector<vector<double> > > xa_vec; // structure of complex1
    vector<vector<vector<double> > > ya_vec; // structure of complex2
    vector<vector<char> >seqx_vec; // sequence of complex1
    vector<vector<char> >seqy_vec; // sequence of complex2
    vector<vector<char> >secx_vec; // secondary structure of complex1
    vector<vector<char> >secy_vec; // secondary structure of complex2
    vector<int> mol_vec1;          // molecule type of complex1, RNA if >0
    vector<int> mol_vec2;          // molecule type of complex2, RNA if >0
    vector<string> chainID_list1;  // list of chainID1
    vector<string> chainID_list2;  // list of chainID2
    vector<int> xlen_vec;          // length of complex1
    vector<int> ylen_vec;          // length of complex2
    int    i,j;                    // chain index
    int    xlen, ylen;             // chain length
    double **xa, **ya;             // structure of single chain
    char   *seqx, *seqy;           // for the protein sequence 
    char   *secx, *secy;           // for the secondary structure 
    int    xlen_aa,ylen_aa;        // total length of protein
    int    xlen_na,ylen_na;        // total length of RNA/DNA
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2

    /* parse complex */
    parse_chain_list(chain1_list, xa_vec, seqx_vec, secx_vec, mol_vec1,
        xlen_vec, chainID_list1, ter_opt, split_opt, mol_opt, infmt1_opt,
        atom_opt, mirror_opt, het_opt, xlen_aa, xlen_na, o_opt, resi_vec1);
    if (xa_vec.size()==0) PrintErrorAndQuit("ERROR! 0 individual chain");
    parse_chain_list(chain2_list, ya_vec, seqy_vec, secy_vec, mol_vec2,
        ylen_vec, chainID_list2, ter_opt, split_opt, mol_opt, infmt2_opt,
        atom_opt, 0, het_opt, ylen_aa, ylen_na, o_opt, resi_vec2);
    if (xa_vec.size()>ya_vec.size()) PrintErrorAndQuit(
        "ERROR! more individual chains to align than number of chains in complex template");
    int len_aa=getmin(xlen_aa,ylen_aa);
    int len_na=getmin(xlen_na,ylen_na);
    if (a_opt)
    {
        len_aa=(xlen_aa+ylen_aa)/2;
        len_na=(xlen_na+ylen_na)/2;
    }

    /* perform monomer alignment if there is only one chain */
    if (xa_vec.size()==1 && ya_vec.size()==1)
    {
        xlen = xlen_vec[0];
        ylen = ylen_vec[0];
        seqx = new char[xlen+1];
        seqy = new char[ylen+1];
        secx = new char[xlen+1];
        secy = new char[ylen+1];
        NewArray(&xa, xlen, 3);
        NewArray(&ya, ylen, 3);
        copy_chain_data(xa_vec[0],seqx_vec[0],secx_vec[0], xlen,xa,seqx,secx);
        copy_chain_data(ya_vec[0],seqy_vec[0],secy_vec[0], ylen,ya,seqy,secy);
        
        /* declare variable specific to this pair of TMalign */
        double t0[3], u0[3][3];
        double TM1, TM2;
        double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
        double d0_0, TM_0;
        double d0A, d0B, d0u, d0a;
        double d0_out=5.0;
        string seqM, seqxA, seqyA;// for output alignment
        double rmsd0 = 0.0;
        int L_ali;                // Aligned length in standard_TMscore
        double Liden=0;
        double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
        int n_ali=0;
        int n_ali8=0;

        /* entry function for structure alignment */
        TMalign_main(xa, ya, seqx, seqy, secx, secy,
            t0, u0, TM1, TM2, TM3, TM4, TM5,
            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
            seqM, seqxA, seqyA,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, Lnorm_ass, d0_scale,
            0, a_opt, u_opt, d_opt, fast_opt,
            mol_vec1[0]+mol_vec2[0],TMcut);

        /* print result */
        output_results(
            xname.substr(dir1_opt.size()),
            yname.substr(dir2_opt.size()),
            chainID_list1[0], chainID_list2[0],
            xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
            seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden,
            n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, d0B,
            Lnorm_ass, d0_scale, d0a, d0u, (m_opt?fname_matrix:"").c_str(),
            (outfmt_opt==2?outfmt_opt:3), ter_opt, true, split_opt, o_opt, fname_super,
            0, a_opt, false, d_opt, mirror_opt, resi_vec1, resi_vec2);
        if (outfmt_opt==2) printf("%s%s\t%s%s\t%.4f\n",
            xname.substr(dir1_opt.size()).c_str(), chainID_list1[0].c_str(), 
            yname.substr(dir2_opt.size()).c_str(), chainID_list2[0].c_str(),
            sqrt((TM1*TM1+TM2*TM2)/2));

        /* clean up */
        seqM.clear();
        seqxA.clear();
        seqyA.clear();
        delete[]seqx;
        delete[]seqy;
        delete[]secx;
        delete[]secy;
        DeleteArray(&xa,xlen);
        DeleteArray(&ya,ylen);

        vector<vector<vector<double> > >().swap(xa_vec); // structure of complex1
        vector<vector<vector<double> > >().swap(ya_vec); // structure of complex2
        vector<vector<char> >().swap(seqx_vec); // sequence of complex1
        vector<vector<char> >().swap(seqy_vec); // sequence of complex2
        vector<vector<char> >().swap(secx_vec); // secondary structure of complex1
        vector<vector<char> >().swap(secy_vec); // secondary structure of complex2
        mol_vec1.clear();       // molecule type of complex1, RNA if >0
        mol_vec2.clear();       // molecule type of complex2, RNA if >0
        chainID_list1.clear();  // list of chainID1
        chainID_list2.clear();  // list of chainID2
        xlen_vec.clear();       // length of complex1
        ylen_vec.clear();       // length of complex2
        return 0;
    }

    /* declare TM-score tables */
    int chain1_num=xa_vec.size();
    int chain2_num=ya_vec.size();
    vector<string> tmp_str_vec(chain2_num,"");
    double **TMave_mat;
    NewArray(&TMave_mat,chain1_num,chain2_num);
    vector<vector<string> >seqxA_mat(chain1_num,tmp_str_vec);
    vector<vector<string> > seqM_mat(chain1_num,tmp_str_vec);
    vector<vector<string> >seqyA_mat(chain1_num,tmp_str_vec);

    /* trimComplex */
    vector<vector<vector<double> > > ya_trim_vec; // structure of complex2
    vector<vector<char> >seqy_trim_vec; // sequence of complex2
    vector<vector<char> >secy_trim_vec; // secondary structure of complex2
    vector<int> ylen_trim_vec;          // length of complex2
    int Lchain_aa_max1=0;
    int Lchain_na_max1=0;
    for (i=0;i<chain1_num;i++)
    {
        xlen=xlen_vec[i];
        if      (mol_vec1[i]>0  && xlen>Lchain_na_max1) Lchain_na_max1=xlen;
        else if (mol_vec1[i]<=0 && xlen>Lchain_aa_max1) Lchain_aa_max1=xlen;
    }
    int trim_chain_count=trimComplex(ya_trim_vec,seqy_trim_vec,
        secy_trim_vec,ylen_trim_vec,ya_vec,seqy_vec,secy_vec,ylen_vec,
        mol_vec2,Lchain_aa_max1,Lchain_na_max1);
    int    ylen_trim;             // chain length
    double **ya_trim;             // structure of single chain
    char   *seqy_trim;           // for the protein sequence
    char   *secy_trim;           // for the secondary structure
    double **xt;

    /* get all-against-all alignment */
    if (len_aa+len_na>500) fast_opt=true;
    for (i=0;i<chain1_num;i++)
    {
        xlen=xlen_vec[i];
        if (xlen<3)
        {
            for (j=0;j<chain2_num;j++) TMave_mat[i][j]=-1;
            continue;
        }
        seqx = new char[xlen+1];
        secx = new char[xlen+1];
        NewArray(&xa, xlen, 3);
        copy_chain_data(xa_vec[i],seqx_vec[i],secx_vec[i],
            xlen,xa,seqx,secx);

        for (j=0;j<chain2_num;j++)
        {
            if (mol_vec1[i]*mol_vec2[j]<0) //no protein-RNA alignment
            {
                TMave_mat[i][j]=-1;
                continue;
            }

            ylen=ylen_vec[j];
            if (ylen<3)
            {
                TMave_mat[i][j]=-1;
                continue;
            }
            seqy = new char[ylen+1];
            secy = new char[ylen+1];
            NewArray(&ya, ylen, 3);
            copy_chain_data(ya_vec[j],seqy_vec[j],secy_vec[j],
                ylen,ya,seqy,secy);

            /* declare variable specific to this pair of TMalign */
            double t0[3], u0[3][3];
            double TM1, TM2;
            double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out=5.0;
            string seqM, seqxA, seqyA;// for output alignment
            double rmsd0 = 0.0;
            int L_ali;                // Aligned length in standard_TMscore
            double Liden=0;
            double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
            int n_ali=0;
            int n_ali8=0;

            int Lnorm_tmp=len_aa;
            if (mol_vec1[i]+mol_vec2[j]>0) Lnorm_tmp=len_na;

            /* entry function for structure alignment */
            if (trim_chain_count && ylen_trim_vec[j]<ylen)
            {
                ylen_trim = ylen_trim_vec[j];
                seqy_trim = new char[ylen_trim+1];
                secy_trim = new char[ylen_trim+1];
                NewArray(&ya_trim, ylen_trim, 3);
                copy_chain_data(ya_trim_vec[j],seqy_trim_vec[j],secy_trim_vec[j],
                    ylen_trim,ya_trim,seqy_trim,secy_trim);
                TMalign_main(xa, ya_trim, seqx, seqy_trim, secx, secy_trim,
                    t0, u0, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen_trim, sequence, Lnorm_tmp, d0_scale,
                    0, false, true, false, fast_opt,
                    mol_vec1[i]+mol_vec2[j],TMcut);
                seqxA.clear();
                seqyA.clear();
                delete[]seqy_trim;
                delete[]secy_trim;
                DeleteArray(&ya_trim,ylen_trim);

                NewArray(&xt,xlen,3);
                do_rotation(xa, xt, xlen, t0, u0);
                int *invmap = new int[ylen+1];
                se_main(xt, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                    0, false, 2, false, mol_vec1[i]+mol_vec2[j], 1, invmap);
                delete[]invmap;
                
                if (sequence.size()<2) sequence.push_back("");
                if (sequence.size()<2) sequence.push_back("");
                sequence[0]=seqxA;
                sequence[1]=seqyA;
                TMalign_main(xt, ya, seqx, seqy, secx, secy,
                    t0, u0, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                    2, false, true, false, fast_opt,
                    mol_vec1[i]+mol_vec2[j],TMcut);
                DeleteArray(&xt, xlen);
            }
            else
            {
                TMalign_main(xa, ya, seqx, seqy, secx, secy,
                    t0, u0, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                    0, false, true, false, fast_opt,
                    mol_vec1[i]+mol_vec2[j],TMcut);
            }
            
            /* store result */
            seqxA_mat[i][j]=seqxA;
            seqyA_mat[i][j]=seqyA;
            TMave_mat[i][j]=TM4*Lnorm_tmp;

            /* clean up */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();

            delete[]seqy;
            delete[]secy;
            DeleteArray(&ya,ylen);
        }

        delete[]seqx;
        delete[]secx;
        DeleteArray(&xa,xlen);
    }
    vector<vector<vector<double> > >().swap(ya_trim_vec);
    vector<vector<char> >().swap(seqy_trim_vec);
    vector<vector<char> >().swap(secy_trim_vec);
    vector<int> ().swap(ylen_trim_vec);

    /* calculate initial chain-chain assignment */
    int *assign1_list; // value is index of assigned chain2
    int *assign2_list; // value is index of assigned chain1
    assign1_list=new int[chain1_num];
    assign2_list=new int[chain2_num];
    enhanced_greedy_search(TMave_mat, assign1_list,
        assign2_list, chain1_num, chain2_num);

    /* final alignment */
    if (outfmt_opt==0) print_version();
    double **ut_mat; // rotation matrices for all-against-all alignment
    NewArray(&ut_mat,chain1_num,4*3);
    int ui,uj;
    vector<string>xname_vec;
    vector<string>yname_vec;
    vector<double>TM_vec;
    for (i=0;i<chain1_num;i++)
    {
        j=assign1_list[i];
        xname_vec.push_back(xname+chainID_list1[i]);
        if (j<0)
        {
            cerr<<"Warning! "<<chainID_list1[i]<<" cannot be alighed"<<endl;
            for (ui=0;ui<3;ui++)
            {
                for (uj=0;uj<4;uj++) ut_mat[i][ui*3+uj]=0;
                ut_mat[i][ui*3+ui]=1;
            }
            yname_vec.push_back(yname);
            continue;
        }
        yname_vec.push_back(yname+chainID_list2[j]);

        xlen =xlen_vec[i];
        seqx = new char[xlen+1];
        secx = new char[xlen+1];
        NewArray(&xa, xlen, 3);
        copy_chain_data(xa_vec[i],seqx_vec[i],secx_vec[i], xlen,xa,seqx,secx);

        ylen =ylen_vec[j];
        seqy = new char[ylen+1];
        secy = new char[ylen+1];
        NewArray(&ya, ylen, 3);
        copy_chain_data(ya_vec[j],seqy_vec[j],secy_vec[j], ylen,ya,seqy,secy);

        /* declare variable specific to this pair of TMalign */
        double t0[3], u0[3][3];
        double TM1, TM2;
        double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
        double d0_0, TM_0;
        double d0A, d0B, d0u, d0a;
        double d0_out=5.0;
        string seqM, seqxA, seqyA;// for output alignment
        double rmsd0 = 0.0;
        int L_ali;                // Aligned length in standard_TMscore
        double Liden=0;
        double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
        int n_ali=0;
        int n_ali8=0;

        int c;
        for (c=0; c<sequence.size(); c++) sequence[c].clear();
        sequence.clear();
        sequence.push_back(seqxA_mat[i][j]);
        sequence.push_back(seqyA_mat[i][j]);
            
        /* entry function for structure alignment */
        TMalign_main(xa, ya, seqx, seqy, secx, secy,
            t0, u0, TM1, TM2, TM3, TM4, TM5,
            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
            seqM, seqxA, seqyA,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, Lnorm_ass, d0_scale,
            3, a_opt, u_opt, d_opt, fast_opt,
            mol_vec1[i]+mol_vec2[j]);
        
        for (ui=0;ui<3;ui++) for (uj=0;uj<3;uj++) ut_mat[i][ui*3+uj]=u0[ui][uj];
        for (uj=0;uj<3;uj++) ut_mat[i][9+uj]=t0[uj];

        TM_vec.push_back(TM1);
        TM_vec.push_back(TM2);

        if (outfmt_opt<2) output_results(
            xname.c_str(), yname.c_str(),
            chainID_list1[i], chainID_list2[j],
            xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5,
            rmsd0, d0_out, seqM.c_str(),
            seqxA.c_str(), seqyA.c_str(), Liden,
            n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
            d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u, 
            "", outfmt_opt, ter_opt, false, split_opt, 
            false, "",//o_opt, fname_super+chainID_list1[i], 
            false, a_opt, u_opt, d_opt, mirror_opt,
            resi_vec1, resi_vec2);
        
        /* clean up */
        seqM.clear();
        seqxA.clear();
        seqyA.clear();

        delete[]seqy;
        delete[]secy;
        DeleteArray(&ya,ylen);

        delete[]seqx;
        delete[]secx;
        DeleteArray(&xa,xlen);
    }
    if (outfmt_opt==2)
    {
        double TM=0;
        for (i=0;i<TM_vec.size();i++) TM+=TM_vec[i]*TM_vec[i];
        TM=sqrt(TM/TM_vec.size());
        string query_name=xname;
        string template_name=yname;
        for (i=0;i<chain1_num;i++)
        {
            j=assign1_list[i];
            if (j<0) continue;
            query_name   +=chainID_list1[i];
            template_name+=chainID_list2[j];
        }
        printf("%s\t%s\t%.4f\n",query_name.c_str(),template_name.c_str(),TM);
        query_name.clear();
        template_name.clear();
    }

    if (m_opt) output_dock_rotation_matrix(fname_matrix.c_str(),
        xname_vec,yname_vec, ut_mat, assign1_list);

    if (o_opt) output_dock(chain1_list, ter_opt, split_opt, infmt1_opt,
        atom_opt, mirror_opt, ut_mat, fname_super);


    /* clean up everything */
    vector<double>().swap(TM_vec);
    vector<string>().swap(xname_vec);
    vector<string>().swap(yname_vec);
    delete [] assign1_list;
    delete [] assign2_list;
    DeleteArray(&TMave_mat,chain1_num);
    DeleteArray(&ut_mat,   chain1_num);
    vector<vector<string> >().swap(seqxA_mat);
    vector<vector<string> >().swap(seqM_mat);
    vector<vector<string> >().swap(seqyA_mat);
    vector<string>().swap(tmp_str_vec);

    vector<vector<vector<double> > >().swap(xa_vec); // structure of complex1
    vector<vector<vector<double> > >().swap(ya_vec); // structure of complex2
    vector<vector<char> >().swap(seqx_vec); // sequence of complex1
    vector<vector<char> >().swap(seqy_vec); // sequence of complex2
    vector<vector<char> >().swap(secx_vec); // secondary structure of complex1
    vector<vector<char> >().swap(secy_vec); // secondary structure of complex2
    mol_vec1.clear();       // molecule type of complex1, RNA if >0
    mol_vec2.clear();       // molecule type of complex2, RNA if >0
    vector<string>().swap(chainID_list1);  // list of chainID1
    vector<string>().swap(chainID_list2);  // list of chainID2
    xlen_vec.clear();       // length of complex1
    ylen_vec.clear();       // length of complex2
    return 1;
}

int mTMalign(string &xname, string &yname, const string &fname_super,
    const string &fname_matrix,
    vector<string> &sequence, double Lnorm_ass, const double d0_scale,
    const bool m_opt, const int  i_opt, const int o_opt, const int a_opt,
    bool u_opt, const bool d_opt, const bool full_opt, const double TMcut,
    const int infmt_opt, const int ter_opt,
    const int split_opt, const int outfmt_opt, bool fast_opt,
    const int het_opt,
    const string &atom_opt, const string &mol_opt, const string &dir_opt,
    const int byresi_opt,
    const vector<string> &chain_list)
{
    /* declare previously global variables */
    vector<vector<vector<double> > >a_vec;  // atomic structure
    vector<vector<vector<double> > >ua_vec; // unchanged atomic structure 
    vector<vector<char> >seq_vec;  // sequence of complex
    vector<vector<char> >sec_vec;  // secondary structure of complex
    vector<int> mol_vec;           // molecule type of complex1, RNA if >0
    vector<string> chainID_list;   // list of chainID
    vector<int> len_vec;           // length of complex
    int    i,j;                    // chain index
    int    xlen, ylen;             // chain length
    double **xa, **ya;             // structure of single chain
    char   *seqx, *seqy;           // for the protein sequence 
    char   *secx, *secy;           // for the secondary structure 
    int    len_aa,len_na;          // total length of protein and RNA/DNA
    vector<string> resi_vec;       // residue index for chain

    /* parse chain list */
    parse_chain_list(chain_list, a_vec, seq_vec, sec_vec, mol_vec,
        len_vec, chainID_list, ter_opt, split_opt, mol_opt, infmt_opt,
        atom_opt, false, het_opt, len_aa, len_na, o_opt, resi_vec);
    int chain_num=a_vec.size();
    if (chain_num<=1) PrintErrorAndQuit("ERROR! <2 chains for multiple alignment");
    if (m_opt||o_opt) for (i=0;i<chain_num;i++) ua_vec.push_back(a_vec[i]);
    int mol_type=0;
    int total_len=0;
    xlen=0;
    for (i=0; i<chain_num; i++)
    {
        if (len_vec[i]>xlen) xlen=len_vec[i];
        total_len+=len_vec[i];
        mol_type+=mol_vec[i];
    }
    if (!u_opt) Lnorm_ass=total_len/chain_num;
    u_opt=true;
    total_len-=xlen;
    if (total_len>750) fast_opt=true;

    /* get all-against-all alignment */
    double **TMave_mat;
    NewArray(&TMave_mat,chain_num,chain_num);
    vector<string> tmp_str_vec(chain_num,"");
    vector<vector<string> >seqxA_mat(chain_num,tmp_str_vec);
    vector<vector<string> >seqyA_mat(chain_num,tmp_str_vec);
    for (i=0;i<chain_num;i++) for (j=0;j<chain_num;j++) TMave_mat[i][j]=0;
    for (i=0;i<chain_num;i++)
    {
        xlen=len_vec[i];
        if (xlen<3) continue;
        seqx = new char[xlen+1];
        secx = new char[xlen+1];
        NewArray(&xa, xlen, 3);
        copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i],xlen,xa,seqx,secx);
        seqxA_mat[i][i]=seqyA_mat[i][i]=(string)(seqx);
        for (j=i+1;j<chain_num;j++)
        {
            ylen=len_vec[j];
            if (ylen<3) continue;
            seqy = new char[ylen+1];
            secy = new char[ylen+1];
            NewArray(&ya, ylen, 3);
            copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j],ylen,ya,seqy,secy);
            
            /* declare variable specific to this pair of TMalign */
            double t0[3], u0[3][3];
            double TM1, TM2;
            double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out=5.0;
            string seqM, seqxA, seqyA;// for output alignment
            double rmsd0 = 0.0;
            int L_ali;                // Aligned length in standard_TMscore
            double Liden=0;
            double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
            int n_ali=0;
            int n_ali8=0;

            /* entry function for structure alignment */
            TMalign_main(xa, ya, seqx, seqy, secx, secy,
                t0, u0, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_ass, d0_scale,
                0, false, u_opt, false, fast_opt,
                mol_type,TMcut);

            /* store result */
            TMave_mat[i][j]=TMave_mat[j][i]=TM4;
            seqxA_mat[i][j]=seqyA_mat[j][i]=seqxA;
            seqyA_mat[i][j]=seqxA_mat[j][i]=seqyA;
            //cout<<chain_list[i]<<':'<<chainID_list[i]
                //<<chain_list[j]<<':'<<chainID_list[j]<<"\tTM4="<<TM4<<endl;
            if (full_opt) output_results(
                chain_list[i],chain_list[j], chainID_list[i], chainID_list[j],
                xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
                seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden,
                n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, d0B,
                Lnorm_ass, d0_scale, d0a, d0u, "",
                outfmt_opt, ter_opt, true, split_opt, o_opt, "",
                0, a_opt, false, d_opt, false, resi_vec, resi_vec);

            /* clean up */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();

            delete[]seqy;
            delete[]secy;
            DeleteArray(&ya,ylen);
        }

        delete[]seqx;
        delete[]secx;
        DeleteArray(&xa,xlen);
    }

    /* representative related variables */   
    int r;
    int repr_idx=0;
    vector<string>xname_vec;
    for (i=0;i<chain_num;i++) xname_vec.push_back(
        chain_list[i].substr(dir_opt.size())+chainID_list[i]);
    vector<string>yname_vec;
    double *TMave_list;
    TMave_list = new double[chain_num];
    int *assign_list;
    assign_list=new int[chain_num];
    vector<string> msa(ylen,""); // row is position along msa; column is sequence

    int compare_num;
    double TM1_total, TM2_total;
    double TM3_total, TM4_total, TM5_total;     // for a_opt, u_opt, d_opt
    double d0_0_total, TM_0_total;
    double d0A_total, d0B_total, d0u_total, d0a_total;
    double d0_out_total;
    double rmsd0_total;
    int L_ali_total;                // Aligned length in standard_TMscore
    double Liden_total;
    double TM_ali_total, rmsd_ali_total;  // TMscore and rmsd in standard_TMscore
    int n_ali_total;
    int n_ali8_total;
    int xlen_total, ylen_total;
    double TM4_total_max=0;

    int max_iter=5-(int)(total_len/200);
    if (max_iter<2) max_iter=2;
    int iter=0;
    vector<double> TM_vec(chain_num,0);
    vector<double> d0_vec(chain_num,0);
    vector<double> seqID_vec(chain_num,0);
    vector<vector<double> > TM_mat(chain_num,TM_vec);
    vector<vector<double> > d0_mat(chain_num,d0_vec);
    vector<vector<double> > seqID_mat(chain_num,seqID_vec);
    for (iter=0; iter<max_iter; iter++)
    {
        /* select representative */   
        for (j=0; j<chain_num; j++) TMave_list[j]=0;
        for (i=0; i<chain_num; i++ )
        {
            for (j=0; j<chain_num; j++)
            {
                //cout<<'\t'<<setprecision(4)<<TMave_mat[i][j];
                TMave_list[j]+=TMave_mat[i][j];
            }
            //cout<<'\t'<<chain_list[i]<<':'<<chainID_list[i]<<endl;
        }
        repr_idx=0;
        double repr_TM=0;
        for (j=0; j<chain_num; j++)
        {
            //cout<<chain_list[j]<<'\t'<<len_vec[j]<<'\t'<<TMave_list[j]<<endl;
            if (TMave_list[j]<repr_TM) continue;
            repr_TM=TMave_list[j];
            repr_idx=j;
        }
        //cout<<"repr="<<repr_idx<<"; "<<chain_list[repr_idx]<<"; TM="<<repr_TM<<endl;

        /* superpose superpose */
        yname=chain_list[repr_idx].substr(dir_opt.size())+chainID_list[repr_idx];
        double **xt;
        vector<pair<double,int> >TM_pair_vec; // TM vs chain

        for (i=0; i<chain_num; i++) assign_list[i]=-1;
        assign_list[repr_idx]=repr_idx;
        //ylen = len_vec[repr_idx];
        //seqy = new char[ylen+1];
        //secy = new char[ylen+1];
        //NewArray(&ya, ylen, 3);
        //copy_chain_data(a_vec[repr_idx],seq_vec[repr_idx],sec_vec[repr_idx], ylen,ya,seqy,secy);
        for (r=0;r<sequence.size();r++) sequence[r].clear(); sequence.clear();
        sequence.push_back("");
        sequence.push_back("");
        for (i=0;i<chain_num;i++)
        {
            yname_vec.push_back(yname);
            xlen = len_vec[i];
            if (i==repr_idx || xlen<3) continue;
            TM_pair_vec.push_back(make_pair(-TMave_mat[i][repr_idx],i));
        }
        sort(TM_pair_vec.begin(),TM_pair_vec.end());
    
        int tm_idx;
        if (outfmt_opt<0) cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
                               <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<endl;
        for (tm_idx=0; tm_idx<TM_pair_vec.size(); tm_idx++)
        {
            i=TM_pair_vec[tm_idx].second;
            xlen = len_vec[i];
            seqx = new char[xlen+1];
            secx = new char[xlen+1];
            NewArray(&xa, xlen, 3);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);

            double maxTM=TMave_mat[i][repr_idx];
            int maxj=repr_idx;
            for (j=0;j<chain_num;j++)
            {
                if (i==j || assign_list[j]<0 || TMave_mat[i][j]<=maxTM) continue;
                maxj=j;
                maxTM=TMave_mat[i][j];
            }
            j=maxj;
            assign_list[i]=j;
            ylen = len_vec[j];
            seqy = new char[ylen+1];
            secy = new char[ylen+1];
            NewArray(&ya, ylen, 3);
            copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j], ylen,ya,seqy,secy);

            sequence[0]=seqxA_mat[i][j];
            sequence[1]=seqyA_mat[i][j];
            //cout<<"tm_idx="<<tm_idx<<"\ti="<<i<<"\tj="<<j<<endl;
            //cout<<"superpose "<<xname_vec[i]<<" to "<<xname_vec[j]<<endl;

            /* declare variable specific to this pair of TMalign */
            double t0[3], u0[3][3];
            double TM1, TM2;
            double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out=5.0;
            string seqM, seqxA, seqyA;// for output alignment
            double rmsd0 = 0.0;
            int L_ali;                // Aligned length in standard_TMscore
            double Liden=0;
            double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
            int n_ali=0;
            int n_ali8=0;

            /* entry function for structure alignment */
            TMalign_main(xa, ya, seqx, seqy, secx, secy,
                t0, u0, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_ass, d0_scale,
                2,  a_opt, u_opt, d_opt, fast_opt, mol_type);
        
            if (outfmt_opt<0) output_results(
                xname_vec[i].c_str(), xname_vec[j].c_str(), "", "",
                xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5,
                rmsd0, d0_out, seqM.c_str(),
                seqxA.c_str(), seqyA.c_str(), Liden,
                n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
                d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u, 
                "", 2,//outfmt_opt,
                ter_opt, false, split_opt, 
                false, "",//o_opt, fname_super+chainID_list1[i], 
                false, a_opt, u_opt, d_opt, false,
                resi_vec, resi_vec);
         
            NewArray(&xt,xlen,3);
            do_rotation(xa, xt, xlen, t0, u0);
            for (r=0;r<xlen;r++)
            {
                a_vec[i][r][0]=xt[r][0];
                a_vec[i][r][1]=xt[r][1];
                a_vec[i][r][2]=xt[r][2];
            }
            DeleteArray(&xt, xlen);
        
            /* clean up */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();
            sequence[0].clear();
            sequence[1].clear();

            delete[]seqx;
            delete[]secx;
            DeleteArray(&xa,xlen);
        
            delete[]seqy;
            delete[]secy;
            DeleteArray(&ya,ylen);
        }
        ylen = len_vec[repr_idx];
        seqy = new char[ylen+1];
        secy = new char[ylen+1];
        NewArray(&ya, ylen, 3);
        copy_chain_data(a_vec[repr_idx],seq_vec[repr_idx],sec_vec[repr_idx], ylen,ya,seqy,secy);

        /* recover alignment */ 
        int    ylen_ext=ylen;        // chain length
        double **ya_ext;             // structure of single chain
        char   *seqy_ext;            // for the protein sequence 
        char   *secy_ext;            // for the secondary structure 
        for (r=0;r<msa.size();r++) msa[r].clear(); msa.clear();
        msa.assign(ylen,""); // row is position along msa; column is sequence
        vector<string> msa_ext;      // row is position along msa; column is sequence
        for (r=0;r<ylen;r++) msa[r]=seqy[r];
        //for (r=0;r<msa.size();r++) cout<<"["<<r<<"]\t"<<msa[r]<<endl;
        //cout<<"start recover"<<endl;
        assign_list[repr_idx]=0;
        for (tm_idx=0; tm_idx<TM_pair_vec.size(); tm_idx++)
        {
            i=TM_pair_vec[tm_idx].second;
            assign_list[i]=tm_idx+1;

            xlen = len_vec[i];
            seqx = new char[xlen+1];
            secx = new char[xlen+1];
            NewArray(&xa, xlen, 3);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);
        
            /* declare variable specific to this pair of TMalign */
            double TM1, TM2;
            double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out=5.0;
            string seqM, seqxA, seqyA;// for output alignment
            double rmsd0 = 0.0;
            int L_ali;                // Aligned length in standard_TMscore
            double Liden=0;
            double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
            int n_ali=0;
            int n_ali8=0;
            int *invmap = new int[ylen+1];

            se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_ass, d0_scale,
                0, a_opt, u_opt, d_opt, mol_type, 1, invmap);

            int rx=0,ry=0;
            ylen_ext=seqxA.size();
            NewArray(&ya_ext, ylen_ext, 3);             // structure of single chain
            seqy_ext= new char[ylen_ext+1];            // for the protein sequence 
            secy_ext= new char[ylen_ext+1];            // for the secondary structure 
            string tmp_gap="";
            for (r=0;r<msa[0].size();r++) tmp_gap+='-';
            for (r=msa_ext.size();r<ylen_ext;r++) msa_ext.push_back("");
            //cout<<"x:"<<xname_vec[i]<<'\n'<<seqxA<<endl;
            //cout<<"y:"<<xname_vec[repr_idx]<<'\n'<<seqyA<<endl;
            for (r=0;r<ylen_ext;r++)
            {
                if (seqyA[r]=='-')
                {
                    msa_ext[r]=tmp_gap+seqxA[r];
                    ya_ext[r][0]=xa[rx][0];
                    ya_ext[r][1]=xa[rx][1];
                    ya_ext[r][2]=xa[rx][2];
                    seqy_ext[r]=seqx[rx];
                    secy_ext[r]=secx[rx];
                }
                else
                {
                    msa_ext[r]=msa[ry]+seqxA[r];
                    ya_ext[r][0]=ya[ry][0];
                    ya_ext[r][1]=ya[ry][1];
                    ya_ext[r][2]=ya[ry][2];
                    seqy_ext[r]=seqy[ry];
                    secy_ext[r]=secy[ry];
                }
                rx+=(seqxA[r]!='-');
                ry+=(seqyA[r]!='-');
            }

            /* copy ya_ext to ya */
            delete[]seqy;
            delete[]secy;
            DeleteArray(&ya,ylen);

            ylen=ylen_ext;
            NewArray(&ya,ylen,3);
            seqy = new char[ylen+1];
            secy = new char[ylen+1];
            for (r=0;r<ylen;r++)
            {
                ya[r][0]=ya_ext[r][0];
                ya[r][1]=ya_ext[r][1];
                ya[r][2]=ya_ext[r][2];
                seqy[r]=seqy_ext[r];
                secy[r]=secy_ext[r];
            }
            for (r=0;r<ylen;r++)
            {
                if (r<msa.size()) msa[r]=msa_ext[r];
                else msa.push_back(msa_ext[r]);
            }
            //for (r=0;r<ylen_ext;r++) cout<<"["<<r<<"]\t"<<msa_ext[r]<<'\t'<<seqy[r]<<'\t'
                    //<<ya[r][0]<<'\t'<<ya[r][1]<<'\t'<<ya[r][2]<<'\t'<<secy[r]<<endl;

            /* clean up */
            tmp_gap.clear();
            delete[]invmap;
            seqM.clear();
            seqxA.clear();
            seqyA.clear();

            delete[]seqx;
            delete[]secx;
            DeleteArray(&xa,xlen);

            delete[]seqy_ext;
            delete[]secy_ext;
            DeleteArray(&ya_ext,ylen_ext);
        }
        vector<string>().swap(msa_ext);
        vector<pair<double,int> >().swap(TM_pair_vec);
        for (i=0; i<chain_num; i++)
        {
            tm_idx=assign_list[i];
            if (tm_idx<0) continue;
            seqyA_mat[i][i]="";
            for (r=0 ;r<ylen ; r++) seqyA_mat[i][i]+=msa[r][tm_idx];
            seqxA_mat[i][i]=seqyA_mat[i][i];
            //cout<<xname_vec[i]<<'\t'<<seqxA_mat[i][i]<<endl;
        }
        for (i=0;i<chain_num; i++)
        {
            if (assign_list[i]<0) continue;
            string seqxA=seqxA_mat[i][i];
            for (j=0; j<chain_num; j++)
            {
                if (i==j || assign_list[j]<0) continue;
                string seqyA=seqyA_mat[j][j];
                seqxA_mat[i][j]=seqyA_mat[i][j]="";
                for (r=0;r<ylen;r++)
                {
                    if (seqxA[r]=='-' && seqyA[r]=='-') continue;
                    seqxA_mat[i][j]+=seqxA[r];
                    seqyA_mat[i][j]+=seqyA[r];
                }
                seqyA.clear();
            }
            seqxA.clear();
        }

        /* recover statistics such as TM-score */ 
        compare_num=0;
        TM1_total=0, TM2_total=0;
        TM3_total=0, TM4_total=0, TM5_total=0;
        d0_0_total=0, TM_0_total=0;
        d0A_total=0, d0B_total=0, d0u_total=0, d0a_total=0;
        d0_out_total=0;
        rmsd0_total = 0.0;
        L_ali_total=0;
        Liden_total=0;
        TM_ali_total=0, rmsd_ali_total=0;
        n_ali_total=0;
        n_ali8_total=0;
        xlen_total=0, ylen_total=0;
        for (i=0; i< chain_num; i++)
        {
            xlen=len_vec[i];
            if (xlen<3) continue;
            seqx = new char[xlen+1];
            secx = new char[xlen+1];
            NewArray(&xa, xlen, 3);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);
            for (j=i+1;j<chain_num;j++)
            {
                ylen=len_vec[j];
                if (ylen<3) continue;
                compare_num++;
                seqy = new char[ylen+1];
                secy = new char[ylen+1];
                NewArray(&ya, ylen, 3);
                copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j],ylen,ya,seqy,secy);
                sequence[0]=seqxA_mat[i][j];
                sequence[1]=seqyA_mat[i][j];
            
                /* declare variable specific to this pair of TMalign */
                double TM1, TM2;
                double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                double d0_0, TM_0;
                double d0A, d0B, d0u, d0a;
                double d0_out=5.0;
                string seqM, seqxA, seqyA;// for output alignment
                double rmsd0 = 0.0;
                int L_ali=0;              // Aligned length in standard_TMscore
                double Liden=0;
                double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                int n_ali=0;
                int n_ali8=0;
                int *invmap = new int[ylen+1];

                se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_ass, d0_scale,
                    true, a_opt, u_opt, d_opt, mol_type, 1, invmap);

                if (xlen<=ylen)
                {
                    xlen_total+=xlen;
                    ylen_total+=ylen;
                    TM1_total+=TM1;
                    TM2_total+=TM2;
                    d0A_total+=d0A;
                    d0B_total+=d0B;
                }
                else
                {
                    xlen_total+=ylen;
                    ylen_total+=xlen;
                    TM1_total+=TM2;
                    TM2_total+=TM1;
                    d0A_total+=d0B;
                    d0B_total+=d0A;
                }
                TM_mat[i][j]=TM2;
                TM_mat[j][i]=TM1;
                d0_mat[i][j]=d0B;
                d0_mat[j][i]=d0A;
                seqID_mat[i][j]=1.*Liden/xlen;
                seqID_mat[j][i]=1.*Liden/ylen;

                TM3_total+=TM3;
                TM4_total+=TM4;
                TM5_total+=TM5;
                d0_0_total+=d0_0;
                TM_0_total+=TM_0;
                d0u_total+=d0u;
                d0_out_total+=d0_out;
                rmsd0_total+=rmsd0;
                L_ali_total+=L_ali;        // Aligned length in standard_TMscore
                Liden_total+=Liden;
                TM_ali_total+=TM_ali;
                rmsd_ali_total+=rmsd_ali;  // TMscore and rmsd in standard_TMscore
                n_ali_total+=n_ali;
                n_ali8_total+=n_ali8;

                /* clean up */
                delete[]invmap;
                seqM.clear();
                seqxA.clear();
                seqyA.clear();

                delete[]seqy;
                delete[]secy;
                DeleteArray(&ya,ylen);
            }
            delete[]seqx;
            delete[]secx;
            DeleteArray(&xa,xlen);
        }
        if (TM4_total<=TM4_total_max) break;
        TM4_total_max=TM4_total;
    }
    for (i=0;i<chain_num;i++)
    {
        for (j=0;j<chain_num;j++)
        {
            if (i==j) continue;
            TM_vec[i]+=TM_mat[i][j];
            d0_vec[i]+=d0_mat[i][j];
            seqID_vec[i]+=seqID_mat[i][j];
        }
        TM_vec[i]/=(chain_num-1);
        d0_vec[i]/=(chain_num-1);
        seqID_vec[i]/=(chain_num-1);
    }
    xlen_total    /=compare_num;
    ylen_total    /=compare_num;
    TM1_total     /=compare_num;
    TM2_total     /=compare_num;
    d0A_total     /=compare_num;
    d0B_total     /=compare_num;
    TM3_total     /=compare_num;
    TM4_total     /=compare_num;
    TM5_total     /=compare_num;
    d0_0_total    /=compare_num;
    TM_0_total    /=compare_num;
    d0u_total     /=compare_num;
    d0_out_total  /=compare_num;
    rmsd0_total   /=compare_num;
    L_ali_total   /=compare_num;
    Liden_total   /=compare_num;
    TM_ali_total  /=compare_num;
    rmsd_ali_total/=compare_num;
    n_ali_total   /=compare_num;
    n_ali8_total  /=compare_num;
    xname="shorter";
    yname="longer";
    string seqM="";
    string seqxA="";
    string seqyA="";
    double t0[3];
    double u0[3][3];
    stringstream buf;
    for (i=0; i<chain_num; i++)
    {
        if (assign_list[i]<0) continue;
        buf <<">"<<xname_vec[i]<<"\tL="<<len_vec[i]
            <<"\td0="<<setiosflags(ios::fixed)<<setprecision(2)<<d0_vec[i]
            <<"\tseqID="<<setiosflags(ios::fixed)<<setprecision(3)<<seqID_vec[i]
            <<"\tTM-score="<<setiosflags(ios::fixed)<<setprecision(5)<<TM_vec[i];
        if (i==repr_idx) buf<<"\t*";
        buf<<'\n'<<seqxA_mat[i][i]<<endl;
    }
    seqM=buf.str();
    seqM=seqM.substr(0,seqM.size()-1);
    buf.str(string());
    //MergeAlign(seqxA_mat,seqyA_mat,repr_idx,xname_vec,chain_num,seqM);
    if (outfmt_opt==0) print_version();
    output_mTMalign_results( xname,yname, "","",
        xlen_total, ylen_total, t0, u0, TM1_total, TM2_total, 
        TM3_total, TM4_total, TM5_total, rmsd0_total, d0_out_total,
        seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden_total,
        n_ali8_total, L_ali_total, TM_ali_total, rmsd_ali_total,
        TM_0_total, d0_0_total, d0A_total, d0B_total,
        Lnorm_ass, d0_scale, d0a_total, d0u_total, 
        "", outfmt_opt, ter_opt, 0, split_opt, false,
        "", false, a_opt, u_opt, d_opt, false,
        resi_vec, resi_vec );

    if (m_opt || o_opt)
    {
        double **ut_mat; // rotation matrices for all-against-all alignment
        int ui,uj;
        double t[3], u[3][3];
        double rmsd;
        NewArray(&ut_mat,chain_num,4*3);
        for (i=0;i<chain_num;i++)
        {
            xlen=ylen=a_vec[i].size();
            NewArray(&xa,xlen,3);
            NewArray(&ya,xlen,3);
            for (r=0;r<xlen;r++)
            {
                xa[r][0]=ua_vec[i][r][0];
                xa[r][1]=ua_vec[i][r][1];
                xa[r][2]=ua_vec[i][r][2];
                ya[r][0]= a_vec[i][r][0];
                ya[r][1]= a_vec[i][r][1];
                ya[r][2]= a_vec[i][r][2];
            }
            Kabsch(xa,ya,xlen,1,&rmsd,t,u);
            for (ui=0;ui<3;ui++) for (uj=0;uj<3;uj++) ut_mat[i][ui*3+uj]=u[ui][uj];
            for (uj=0;uj<3;uj++) ut_mat[i][9+uj]=t[uj];
            DeleteArray(&xa,xlen);
            DeleteArray(&ya,xlen);
        }
        vector<vector<vector<double> > >().swap(ua_vec);

        if (m_opt)
        {
            assign_list[repr_idx]=-1;
            output_dock_rotation_matrix(fname_matrix.c_str(),
                xname_vec,yname_vec, ut_mat, assign_list);
        }

        if (o_opt) output_dock(chain_list, ter_opt, split_opt, 
                infmt_opt, atom_opt, false, ut_mat, fname_super);
        
        DeleteArray(&ut_mat,chain_num);
    }

    /* clean up */
    vector<string>().swap(msa);
    vector<string>().swap(tmp_str_vec);
    vector<vector<string> >().swap(seqxA_mat);
    vector<vector<string> >().swap(seqyA_mat);
    vector<string>().swap(xname_vec);
    vector<string>().swap(yname_vec);
    delete[]TMave_list;
    DeleteArray(&TMave_mat,chain_num);
    vector<vector<vector<double> > >().swap(a_vec); // structure of complex
    vector<vector<char> >().swap(seq_vec); // sequence of complex
    vector<vector<char> >().swap(sec_vec); // secondary structure of complex
    vector<int>().swap(mol_vec);           // molecule type of complex1, RNA if >0
    vector<string>().swap(chainID_list);   // list of chainID
    vector<int>().swap(len_vec);           // length of complex
    vector<double>().swap(TM_vec);
    vector<double>().swap(d0_vec);
    vector<double>().swap(seqID_vec);
    vector<vector<double> >().swap(TM_mat);
    vector<vector<double> >().swap(d0_mat);
    vector<vector<double> >().swap(seqID_mat);
    return 1;
}

/* sequence order independent alignment */
int SOIalign(string &xname, string &yname, const string &fname_super,
    const string &fname_lign, const string &fname_matrix,
    vector<string> &sequence, const double Lnorm_ass, const double d0_scale,
    const bool m_opt, const int  i_opt, const int o_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const double TMcut,
    const int infmt1_opt, const int infmt2_opt, const int ter_opt,
    const int split_opt, const int outfmt_opt, const bool fast_opt,
    const int cp_opt, const int mirror_opt, const int het_opt,
    const string &atom_opt, const string &mol_opt, const string &dir_opt,
    const string &dir1_opt, const string &dir2_opt, 
    const vector<string> &chain1_list, const vector<string> &chain2_list,
    const bool se_opt, const int closeK_opt, const int mm_opt)
{
    /* declare previously global variables */
    vector<vector<string> >PDB_lines1; // text of chain1
    vector<vector<string> >PDB_lines2; // text of chain2
    vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    vector<string> chainID_list1;      // list of chainID1
    vector<string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    r;                  // residue index
    int    xlen, ylen;         // chain length
    int    xchainnum,ychainnum;// number of chains in a PDB file
    char   *seqx, *seqy;       // for the protein sequence 
    char   *secx, *secy;       // for the secondary structure 
    int    **secx_bond;        // boundary of secondary structure
    int    **secy_bond;        // boundary of secondary structure
    double **xa, **ya;         // for input vectors xa[0...xlen-1][0..2] and
                               // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure 
                               // --> superpose xa onto ya
    double **xk, **yk;         // k closest residues
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2
    int read_resi=0;  // whether to read residue index
    if (o_opt) read_resi=2;

    /* loop over file names */
    for (i=0;i<chain1_list.size();i++)
    {
        /* parse chain 1 */
        xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1,
            mol_vec1, ter_opt, infmt1_opt, atom_opt, split_opt, het_opt);
        if (!xchainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=PDB_lines1[chain_i].size();
            if (mol_opt=="RNA") mol_vec1[chain_i]=1;
            else if (mol_opt=="protein") mol_vec1[chain_i]=-1;
            if (!xlen)
            {
                cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<endl;
                continue;
            }
            else if (xlen<3)
            {
                cerr<<"Sequence is too short <3!: "<<xname<<endl;
                continue;
            }
            NewArray(&xa, xlen, 3);
            if (closeK_opt>=3) NewArray(&xk, xlen*closeK_opt, 3);
            seqx = new char[xlen + 1];
            secx = new char[xlen + 1];
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx, 
                resi_vec1, read_resi);
            if (mirror_opt) for (r=0;r<xlen;r++) xa[r][2]=-xa[r][2];
            if (mol_vec1[chain_i]>0) make_sec(seqx,xa, xlen, secx,atom_opt);
            else make_sec(xa, xlen, secx); // secondary structure assignment
            if (closeK_opt>=3) getCloseK(xa, xlen, closeK_opt, xk);
            if (mm_opt==6) 
            {
                NewArray(&secx_bond, xlen, 2);
                assign_sec_bond(secx_bond, secx, xlen);
            }

            for (j=(dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                /* parse chain 2 */
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2, chainID_list2,
                        mol_vec2, ter_opt, infmt2_opt, atom_opt, split_opt,
                        het_opt);
                    if (!ychainnum)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain number 0."<<endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ylen=PDB_lines2[chain_j].size();
                    if (mol_opt=="RNA") mol_vec2[chain_j]=1;
                    else if (mol_opt=="protein") mol_vec2[chain_j]=-1;
                    if (!ylen)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain length 0."<<endl;
                        continue;
                    }
                    else if (ylen<3)
                    {
                        cerr<<"Sequence is too short <3!: "<<yname<<endl;
                        continue;
                    }
                    NewArray(&ya, ylen, 3);
                    if (closeK_opt>=3) NewArray(&yk, ylen*closeK_opt, 3);
                    seqy = new char[ylen + 1];
                    secy = new char[ylen + 1];
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy,
                        resi_vec2, read_resi);
                    if (mol_vec2[chain_j]>0)
                         make_sec(seqy, ya, ylen, secy, atom_opt);
                    else make_sec(ya, ylen, secy);
                    if (closeK_opt>=3) getCloseK(ya, ylen, closeK_opt, yk);
                    if (mm_opt==6) 
                    {
                        NewArray(&secy_bond, ylen, 2);
                        assign_sec_bond(secy_bond, secy, ylen);
                    }

                    /* declare variable specific to this pair of TMalign */
                    double t0[3], u0[3][3];
                    double TM1, TM2;
                    double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                    double d0_0, TM_0;
                    double d0A, d0B, d0u, d0a;
                    double d0_out=5.0;
                    string seqM, seqxA, seqyA;// for output alignment
                    double rmsd0 = 0.0;
                    int L_ali;                // Aligned length in standard_TMscore
                    double Liden=0;
                    double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                    int n_ali=0;
                    int n_ali8=0;
                    bool force_fast_opt=(getmin(xlen,ylen)>1500)?true:fast_opt;
                    int *invmap = new int[ylen+1];
                    double *dist_list = new double[ylen+1];

                    /* entry function for structure alignment */
                    if (se_opt) 
                    {
                        u0[0][0]=u0[1][1]=u0[2][2]=1;
                        u0[0][1]=         u0[0][2]=
                        u0[1][0]=         u0[1][2]=
                        u0[2][0]=         u0[2][1]=
                        t0[0]   =t0[1]   =t0[2]   =0;
                        soi_se_main(
                            xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                            seqM, seqxA, seqyA,
                            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                            xlen, ylen, Lnorm_ass, d0_scale,
                            i_opt, a_opt, u_opt, d_opt,
                            mol_vec1[chain_i]+mol_vec2[chain_j], 
                            outfmt_opt, invmap, dist_list,
                            secx_bond, secy_bond, mm_opt);
                        if (outfmt_opt>=2) 
                        {
                            Liden=L_ali=0;
                            int r1,r2;
                            for (r2=0;r2<ylen;r2++)
                            {
                                r1=invmap[r2];
                                if (r1<0) continue;
                                L_ali+=1;
                                Liden+=(seqx[r1]==seqy[r2]);
                            }
                        }
                    }
                    else SOIalign_main(xa, ya, xk, yk, closeK_opt,
                        seqx, seqy, secx, secy,
                        t0, u0, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA, invmap,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                        mol_vec1[chain_i]+mol_vec2[chain_j], dist_list,
                        secx_bond, secy_bond, mm_opt);

                    /* print result */
                    if (outfmt_opt==0) print_version();
                    output_results(
                        xname.substr(dir1_opt.size()+dir_opt.size()),
                        yname.substr(dir2_opt.size()+dir_opt.size()),
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5,
                        rmsd0, d0_out, seqM.c_str(),
                        seqxA.c_str(), seqyA.c_str(), Liden,
                        n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
                        d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u, 
                        (m_opt?fname_matrix:"").c_str(),
                        outfmt_opt, ter_opt, false, split_opt, o_opt,
                        fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
                        resi_vec1, resi_vec2);
                    if (outfmt_opt<=0)
                    {
                        cout<<"###############\t###############\t#########"<<endl;
                        cout<<"#Aligned atom 1\tAligned atom 2 \tDistance#"<<endl;
                        int r1,r2;
                        for (r2=0;r2<ylen;r2++)
                        {
                            r1=invmap[r2];
                            if (r1<0) continue;
                            cout<<PDB_lines1[chain_i][r1].substr(12,15)<<'\t'
                                <<PDB_lines2[chain_j][r2].substr(12,15)<<'\t'
                                <<setw(9)<<setiosflags(ios::fixed)<<setprecision(3)
                                <<dist_list[r2]<<'\n';
                        }
                        cout<<"###############\t###############\t#########"<<endl;
                    }

                    /* Done! Free memory */
                    delete [] invmap;
                    delete [] dist_list;
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();
                    DeleteArray(&ya, ylen);
                    if (closeK_opt>=3) DeleteArray(&yk, ylen*closeK_opt);
                    delete [] seqy;
                    delete [] secy;
                    resi_vec2.clear();
                    if (mm_opt==6) DeleteArray(&secy_bond, ylen);
                } // chain_j
                if (chain2_list.size()>1)
                {
                    yname.clear();
                    for (chain_j=0;chain_j<ychainnum;chain_j++)
                        PDB_lines2[chain_j].clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
            PDB_lines1[chain_i].clear();
            DeleteArray(&xa, xlen);
            if (closeK_opt>=3) DeleteArray(&xk, xlen*closeK_opt);
            delete [] seqx;
            delete [] secx;
            resi_vec1.clear();
            if (mm_opt==6) DeleteArray(&secx_bond, xlen);
        } // chain_i
        xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
    if (chain2_list.size()==1)
    {
        yname.clear();
        for (chain_j=0;chain_j<ychainnum;chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        resi_vec2.clear();
        chainID_list2.clear();
        mol_vec2.clear();
    }
    return 0;
}

int flexalign(string &xname, string &yname, const string &fname_super,
    const string &fname_lign, const string &fname_matrix,
    vector<string> &sequence, const double Lnorm_ass, const double d0_scale,
    const bool m_opt, const int  i_opt, const int o_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const double TMcut,
    const int infmt1_opt, const int infmt2_opt, const int ter_opt,
    const int split_opt, const int outfmt_opt, const bool fast_opt,
    const int mirror_opt, const int het_opt,
    const string &atom_opt, const string &mol_opt, const string &dir_opt,
    const string &dir1_opt, const string &dir2_opt, const int byresi_opt,
    const vector<string> &chain1_list, const vector<string> &chain2_list,
    const int hinge_opt)
{
    /* declare previously global variables */
    vector<vector<string> >PDB_lines1; // text of chain1
    vector<vector<string> >PDB_lines2; // text of chain2
    vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    vector<string> chainID_list1;      // list of chainID1
    vector<string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    r;                  // residue index
    int    xlen, ylen;         // chain length
    int    xchainnum,ychainnum;// number of chains in a PDB file
    char   *seqx, *seqy;       // for the protein sequence 
    char   *secx, *secy;       // for the secondary structure 
    double **xa, **ya;         // for input vectors xa[0...xlen-1][0..2] and
                               // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure 
                               // --> superpose xa onto ya
    vector<string> resi_vec1;  // residue index for chain1
    vector<string> resi_vec2;  // residue index for chain2
    int read_resi=byresi_opt;  // whether to read residue index
    if (byresi_opt==0 && o_opt) read_resi=2;

    /* loop over file names */
    for (i=0;i<chain1_list.size();i++)
    {
        /* parse chain 1 */
        xname=chain1_list[i];
        xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1,
            mol_vec1, ter_opt, infmt1_opt, atom_opt, split_opt, het_opt);
        if (!xchainnum)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain number 0."<<endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=PDB_lines1[chain_i].size();
            if (mol_opt=="RNA") mol_vec1[chain_i]=1;
            else if (mol_opt=="protein") mol_vec1[chain_i]=-1;
            if (!xlen)
            {
                cerr<<"Warning! Cannot parse file: "<<xname
                    <<". Chain length 0."<<endl;
                continue;
            }
            else if (xlen<3)
            {
                cerr<<"Sequence is too short <3!: "<<xname<<endl;
                continue;
            }
            NewArray(&xa, xlen, 3);
            seqx = new char[xlen + 1];
            secx = new char[xlen + 1];
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx, 
                resi_vec1, read_resi);
            if (mirror_opt) for (r=0;r<xlen;r++) xa[r][2]=-xa[r][2];
            if (mol_vec1[chain_i]>0) make_sec(seqx,xa, xlen, secx,atom_opt);
            else make_sec(xa, xlen, secx); // secondary structure assignment

            for (j=(dir_opt.size()>0)*(i+1);j<chain2_list.size();j++)
            {
                /* parse chain 2 */
                if (PDB_lines2.size()==0)
                {
                    yname=chain2_list[j];
                    ychainnum=get_PDB_lines(yname, PDB_lines2, chainID_list2,
                        mol_vec2, ter_opt, infmt2_opt, atom_opt, split_opt,
                        het_opt);
                    if (!ychainnum)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain number 0."<<endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ylen=PDB_lines2[chain_j].size();
                    if (mol_opt=="RNA") mol_vec2[chain_j]=1;
                    else if (mol_opt=="protein") mol_vec2[chain_j]=-1;
                    if (!ylen)
                    {
                        cerr<<"Warning! Cannot parse file: "<<yname
                            <<". Chain length 0."<<endl;
                        continue;
                    }
                    else if (ylen<3)
                    {
                        cerr<<"Sequence is too short <3!: "<<yname<<endl;
                        continue;
                    }
                    NewArray(&ya, ylen, 3);
                    seqy = new char[ylen + 1];
                    secy = new char[ylen + 1];
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy,
                        resi_vec2, read_resi);
                    if (mol_vec2[chain_j]>0)
                         make_sec(seqy, ya, ylen, secy, atom_opt);
                    else make_sec(ya, ylen, secy);

                    if (byresi_opt) extract_aln_from_resi(sequence,
                        seqx,seqy,resi_vec1,resi_vec2,byresi_opt);

                    /* declare variable specific to this pair of TMalign */
                    double t0[3], u0[3][3];
                    double TM1, TM2;
                    double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
                    double d0_0, TM_0;
                    double d0A, d0B, d0u, d0a;
                    double d0_out=5.0;
                    string seqM, seqxA, seqyA;// for output alignment
                    double rmsd0 = 0.0;
                    int L_ali;                // Aligned length in standard_TMscore
                    double Liden=0;
                    double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
                    int n_ali=0;
                    int n_ali8=0;
                    bool force_fast_opt=(getmin(xlen,ylen)>1500)?true:fast_opt;
                    vector<vector<double> >tu_vec;

                    /* entry function for structure alignment */
                    int hingeNum=flexalign_main(
                        xa, ya, seqx, seqy, secx, secy,
                        t0, u0, tu_vec, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                        mol_vec1[chain_i]+mol_vec2[chain_j],hinge_opt);
                    
                    if (hinge_opt && hingeNum<=1 &&
                        n_ali8<0.6*getmin(xlen,ylen))
                    {
                        double t0_h[3], u0_h[3][3];
                        double TM1_h, TM2_h;
                        double TM3_h, TM4_h, TM5_h;
                        double d0_0_h, TM_0_h;
                        double d0_out_h=5.0;
                        string seqM_h, seqxA_h, seqyA_h;
                        double rmsd0_h = 0.0;
                        int L_ali_h;
                        double Liden_h=0;
                        double TM_ali_h, rmsd_ali_h;
                        int n_ali_h=0;
                        int n_ali8_h=0;
                        vector<vector<double> >tu_vec_h(1,tu_vec[0]);
                        tu2t_u(tu_vec[0],t0_h,u0_h);

                        int hingeNum_h=flexalign_main(
                            xa, ya, seqx, seqy, secx, secy,
                            t0_h, u0_h, tu_vec_h,
                            TM1_h, TM2_h, TM3_h, TM4_h, TM5_h,
                            d0_0_h, TM_0_h, d0A, d0B, d0u, d0a, d0_out_h,
                            seqM_h, seqxA_h, seqyA_h, rmsd0_h, L_ali_h,
                            Liden_h, TM_ali_h, rmsd_ali_h, n_ali_h, n_ali8_h,
                            xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
                            a_opt, u_opt, d_opt, force_fast_opt,
                            mol_vec1[chain_i]+mol_vec2[chain_j],hinge_opt);
                        
                        double TM  =(TM1  >TM2  )?TM1  :TM2;
                        double TM_h=(TM1_h>TM2_h)?TM1_h:TM2_h;
                        if (TM_h>TM)
                        {
                            hingeNum=hingeNum_h;
                            tu2t_u(tu_vec_h[0],t0,u0);
                            TM1=TM1_h;
                            TM2=TM2_h;
                            TM3=TM3_h;
                            TM4=TM4_h;
                            TM5=TM5_h;
                            d0_0=d0_0_h;
                            TM_0=TM_0_h;
                            d0_out=d0_out_h;
                            seqM=seqM_h;
                            seqxA=seqxA_h;
                            seqyA=seqyA_h;
                            rmsd0=rmsd0_h;
                            L_ali=L_ali_h;
                            Liden=Liden_h;
                            TM_ali=TM_ali_h;
                            rmsd_ali=rmsd_ali_h;
                            n_ali=n_ali_h;
                            n_ali8=n_ali8_h;
                            tu_vec.clear();
                            for (int hinge=0;hinge<tu_vec_h.size();hinge++)
                                tu_vec.push_back(tu_vec_h[hinge]);
                        }
                        else tu2t_u(tu_vec[0],t0,u0);
                    }

                    /* print result */
                    if (outfmt_opt==0) print_version();
                    output_flexalign_results(
                        xname.substr(dir1_opt.size()+dir_opt.size()),
                        yname.substr(dir2_opt.size()+dir_opt.size()),
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        xlen, ylen, t0, u0, tu_vec, TM1, TM2, TM3, TM4, TM5,
                        rmsd0, d0_out, seqM.c_str(),
                        seqxA.c_str(), seqyA.c_str(), Liden,
                        n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
                        d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u, 
                        (m_opt?fname_matrix:"").c_str(),
                        outfmt_opt, ter_opt, false, split_opt, o_opt,
                        fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
                        resi_vec1, resi_vec2);

                    /* Done! Free memory */
                    tu_vec.clear();
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();
                    DeleteArray(&ya, ylen);
                    delete [] seqy;
                    delete [] secy;
                    resi_vec2.clear();
                } // chain_j
                if (chain2_list.size()>1)
                {
                    yname.clear();
                    for (chain_j=0;chain_j<ychainnum;chain_j++)
                        PDB_lines2[chain_j].clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
            PDB_lines1[chain_i].clear();
            DeleteArray(&xa, xlen);
            delete [] seqx;
            delete [] secx;
            resi_vec1.clear();
        } // chain_i
        xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
    if (chain2_list.size()==1)
    {
        yname.clear();
        for (chain_j=0;chain_j<ychainnum;chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        resi_vec2.clear();
        chainID_list2.clear();
        mol_vec2.clear();
    }
    return 0;
}


int main(int argc, char *argv[])
{
    if (argc < 2) print_help();


    clock_t t1, t2;
    t1 = clock();

    /**********************/
    /*    get argument    */
    /**********************/
    string xname       = "";
    string yname       = "";
    string fname_super = ""; // file name for superposed structure
    string fname_lign  = ""; // file name for user alignment
    string fname_matrix= ""; // file name for output matrix
    vector<string> sequence; // get value from alignment file
    double Lnorm_ass, d0_scale;

    bool h_opt = false; // print full help message
    bool v_opt = false; // print version
    bool m_opt = false; // flag for -m, output rotation matrix
    int  i_opt = 0;     // 1 for -i, 3 for -I
    int  o_opt = 0;     // 1 for -o, 2 for -rasmol
    int  a_opt = 0;     // flag for -a, do not normalized by average length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0

    bool   full_opt  = false;// do not show chain level alignment
    double TMcut     =-1;
    bool   se_opt    =false;
    int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
    int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
    int    ter_opt   =2;     // END, or different chainID
    int    split_opt =2;     // split each chains
    int    outfmt_opt=0;     // set -outfmt to full output
    bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
    int    cp_opt    =0;     // do not check circular permutation
    int    closeK_opt=-1;    // number of atoms for SOI initial alignment.
                             // 5 and 0 for -mm 5 and 6
    int    hinge_opt =9;     // maximum number of hinge allowed for flexible
    int    mirror_opt=0;     // do not align mirror
    int    het_opt=0;        // do not read HETATM residues
    int    mm_opt=0;         // do not perform MM-align
    string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    string suffix_opt="";    // set -suffix to empty
    string dir_opt   ="";    // set -dir to empty
    string dir1_opt  ="";    // set -dir1 to empty
    string dir2_opt  ="";    // set -dir2 to empty
    int    byresi_opt=0;     // set -byresi to 0
    vector<string> chain1_list; // only when -dir1 is set
    vector<string> chain2_list; // only when -dir2 is set

    for(int i = 1; i < argc; i++)
    {
        if ( !strcmp(argv[i],"-o") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -o");
            if (o_opt==2)
                cerr<<"Warning! -rasmol is already set. Ignore -o"<<endl;
            else
            {
                fname_super = argv[i + 1];
                o_opt = 1;
            }
            i++;
        }
        else if ( !strcmp(argv[i],"-rasmol") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -rasmol");
            if (o_opt==1)
                cerr<<"Warning! -o is already set. Ignore -rasmol"<<endl;
            else
            {
                fname_super = argv[i + 1];
                o_opt = 2;
            }
            i++;
        }
        else if ( !strcmp(argv[i],"-u") || !strcmp(argv[i],"-L") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -u or -L");
            Lnorm_ass = atof(argv[i + 1]); u_opt = true; i++;
            if (Lnorm_ass<=0) PrintErrorAndQuit(
                "ERROR! The value for -u or -L should be >0");
        }
        else if ( !strcmp(argv[i],"-a") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -a");
            if (!strcmp(argv[i + 1], "T"))      a_opt=true;
            else if (!strcmp(argv[i + 1], "F")) a_opt=false;
            else 
            {
                a_opt=atoi(argv[i + 1]);
                if (a_opt!=-2 && a_opt!=-1 && a_opt!=1)
                    PrintErrorAndQuit("-a must be -2, -1, 1, T or F");
            }
            i++;
        }
        else if ( !strcmp(argv[i],"-full") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -full");
            if (!strcmp(argv[i + 1], "T"))      full_opt=true;
            else if (!strcmp(argv[i + 1], "F")) full_opt=false;
            else PrintErrorAndQuit("-full must be T or F");
            i++;
        }
        else if ( !strcmp(argv[i],"-d") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -d");
            d0_scale = atof(argv[i + 1]); d_opt = true; i++;
        }
        else if ( !strcmp(argv[i],"-closeK") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -closeK");
            closeK_opt = atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-hinge") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -hinge");
            hinge_opt = atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-v") )
        {
            v_opt = true;
        }
        else if ( !strcmp(argv[i],"-h") )
        {
            h_opt = true;
        }
        else if ( !strcmp(argv[i],"-i") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -i");
            if (i_opt==3)
                PrintErrorAndQuit("ERROR! -i and -I cannot be used together");
            fname_lign = argv[i + 1];      i_opt = 1; i++;
        }
        else if (!strcmp(argv[i], "-I") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -I");
            if (i_opt==1)
                PrintErrorAndQuit("ERROR! -I and -i cannot be used together");
            fname_lign = argv[i + 1];      i_opt = 3; i++;
        }
        else if (!strcmp(argv[i], "-m") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -m");
            fname_matrix = argv[i + 1];    m_opt = true; i++;
        }// get filename for rotation matrix
        else if (!strcmp(argv[i], "-fast"))
        {
            fast_opt = true;
        }
        else if (!strcmp(argv[i], "-se"))
        {
            se_opt = true;
        }
        else if ( !strcmp(argv[i],"-infmt1") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -infmt1");
            infmt1_opt=atoi(argv[i + 1]); i++;
            if (infmt1_opt<-1 || infmt1_opt>3)
                PrintErrorAndQuit("ERROR! -infmt1 can only be -1, 0, 1, 2, or 3");
        }
        else if ( !strcmp(argv[i],"-infmt2") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -infmt2");
            infmt2_opt=atoi(argv[i + 1]); i++;
            if (infmt2_opt<-1 || infmt2_opt>3)
                PrintErrorAndQuit("ERROR! -infmt2 can only be -1, 0, 1, 2, or 3");
        }
        else if ( !strcmp(argv[i],"-ter") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -ter");
            ter_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-split") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -split");
            split_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-atom") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -atom");
            atom_opt=argv[i + 1]; i++;
            if (atom_opt.size()!=4) PrintErrorAndQuit(
                "ERROR! Atom name must have 4 characters, including space.\n"
                "For example, C alpha, C3' and P atoms should be specified by\n"
                "-atom \" CA \", -atom \" P  \" and -atom \" C3'\", respectively.");
        }
        else if ( !strcmp(argv[i],"-mol") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -mol");
            mol_opt=argv[i + 1]; i++;
            if (mol_opt=="prot") mol_opt="protein";
            else if (mol_opt=="DNA") mol_opt="RNA";
            if (mol_opt!="auto" && mol_opt!="protein" && mol_opt!="RNA")
                PrintErrorAndQuit("ERROR! Molecule type must be one of the "
                    "following:\nauto, prot (the same as 'protein'), and "
                    "RNA (the same as 'DNA').");
        }
        else if ( !strcmp(argv[i],"-dir") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -dir");
            dir_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-dir1") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -dir1");
            dir1_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-dir2") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -dir2");
            dir2_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-suffix") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -suffix");
            suffix_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-outfmt") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -outfmt");
            outfmt_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-TMcut") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -TMcut");
            TMcut=atof(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-byresi")  || 
                  !strcmp(argv[i],"-tmscore") ||
                  !strcmp(argv[i],"-TMscore"))
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -byresi");
            byresi_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-seq") )
        {
            byresi_opt=5;
        }
        else if ( !strcmp(argv[i],"-cp") )
        {
            mm_opt=3;
        }
        else if ( !strcmp(argv[i],"-mirror") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -mirror");
            mirror_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-het") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -het");
            het_opt=atoi(argv[i + 1]); i++;
            if (het_opt!=0 && het_opt!=1 && het_opt!=2)
                PrintErrorAndQuit("-het must be 0, 1, or 2");
        }
        else if ( !strcmp(argv[i],"-mm") )
        {
            if (i>=(argc-1)) 
                PrintErrorAndQuit("ERROR! Missing value for -mm");
            mm_opt=atoi(argv[i + 1]); i++;
        }
        else if (xname.size() == 0) xname=argv[i];
        else if (yname.size() == 0) yname=argv[i];
        else PrintErrorAndQuit(string("ERROR! Undefined option ")+argv[i]);
    }

    if(xname.size()==0 || (yname.size()==0 && dir_opt.size()==0) || 
                          (yname.size()    && dir_opt.size()))
    {
        if (h_opt) print_help(h_opt);
        if (v_opt)
        {
            print_version();
            exit(EXIT_FAILURE);
        }
        if (xname.size()==0)
            PrintErrorAndQuit("Please provide input structures");
        else if (yname.size()==0 && dir_opt.size()==0 && mm_opt!=4)
            PrintErrorAndQuit("Please provide structure B");
        else if (yname.size() && dir_opt.size())
            PrintErrorAndQuit("Please provide only one file name if -dir is set");
    }

    if (suffix_opt.size() && dir_opt.size()+dir1_opt.size()+dir2_opt.size()==0)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
    if ((dir_opt.size() || dir1_opt.size() || dir2_opt.size()))
    {
        if (mm_opt!=2 && mm_opt!=4)
        {
            if (o_opt)
                PrintErrorAndQuit("-o cannot be set with -dir, -dir1 or -dir2");
            if (m_opt && fname_matrix!="-")
                PrintErrorAndQuit("-m can only be - or unset when using -dir, -dir1 or -dir2");
        }
        else if (dir_opt.size() && (dir1_opt.size() || dir2_opt.size()))
            PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");
    }
    if (o_opt && (infmt1_opt!=-1 && infmt1_opt!=0 && infmt1_opt!=3))
        PrintErrorAndQuit("-o can only be used with -infmt1 -1, 0 or 3");

    if (mol_opt=="protein" && atom_opt=="auto")
        atom_opt=" CA ";
    else if (mol_opt=="RNA" && atom_opt=="auto")
        atom_opt=" C3'";

    if (d_opt && d0_scale<=0)
        PrintErrorAndQuit("Wrong value for option -d! It should be >0");
    if (outfmt_opt>=2 && (a_opt || u_opt || d_opt))
        PrintErrorAndQuit("-outfmt 2 cannot be used with -a, -u, -L, -d");
    if (byresi_opt!=0)
    {
        if (i_opt)
            PrintErrorAndQuit("-byresi >=1 cannot be used with -i or -I");
        if (byresi_opt<0 || byresi_opt>6)
            PrintErrorAndQuit("-byresi can only be 0 to 6");
        if ((byresi_opt==2 || byresi_opt==3 || byresi_opt==6) && ter_opt>=2)
            PrintErrorAndQuit("-byresi 2 and 6 must be used with -ter <=1");
    }
    //if (split_opt==1 && ter_opt!=0)
        //PrintErrorAndQuit("-split 1 should be used with -ter 0");
    //else if (split_opt==2 && ter_opt!=0 && ter_opt!=1)
        //PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (split_opt<0 || split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");

    if (mm_opt==3)
    {
        cp_opt=true;
        mm_opt=0;
    }
    if (cp_opt && i_opt)
        PrintErrorAndQuit("-mm 3 cannot be used with -i or -I");

    if (mirror_opt && het_opt!=1)
        cerr<<"WARNING! -mirror was not used with -het 1. "
            <<"D amino acids may not be correctly aligned."<<endl;

    if (mm_opt)
    {
        if (i_opt) PrintErrorAndQuit("-mm cannot be used with -i or -I");
        if (u_opt) PrintErrorAndQuit("-mm cannot be used with -u or -L");
        //if (cp_opt) PrintErrorAndQuit("-mm cannot be used with -cp");
        if (dir_opt.size() && (mm_opt==1||mm_opt==2)) PrintErrorAndQuit("-mm 1 or 2 cannot be used with -dir");
        if (byresi_opt) PrintErrorAndQuit("-mm cannot be used with -byresi");
        if (ter_opt>=2 && (mm_opt==1 || mm_opt==2)) PrintErrorAndQuit("-mm 1 or 2 must be used with -ter 0 or -ter 1");
        if (mm_opt==4 && (yname.size() || dir2_opt.size()))
            cerr<<"WARNING! structure_2 is ignored for -mm 4"<<endl;
    }
    else if (full_opt) PrintErrorAndQuit("-full can only be used with -mm");

    if (o_opt && ter_opt<=1 && split_opt==2)
    {
        if (mm_opt && o_opt==2) cerr<<"WARNING! -mm may generate incorrect" 
            <<" RasMol output due to limitations in PDB file format. "
            <<"When -mm is used, -o is recommended over -rasmol"<<endl;
        else if (mm_opt==0) cerr<<"WARNING! Only the superposition of the"
            <<"last aligned chain pair will be generated"<<endl;
    }

    if (closeK_opt<0)
    {
        if (mm_opt==5) closeK_opt=5;
        else closeK_opt=0;
    }

    if (mm_opt==7 && hinge_opt>=10)
        PrintErrorAndQuit("ERROR! -hinge must be <10");


    /* read initial alignment file from 'align.txt' */
    if (i_opt) read_user_alignment(sequence, fname_lign, i_opt);

    if (byresi_opt==6) mm_opt=1;
    else if (byresi_opt) i_opt=3;

    if (m_opt && fname_matrix == "") // Output rotation matrix: matrix.txt
        PrintErrorAndQuit("ERROR! Please provide a file name for option -m!");

    /* parse file list */
    if (dir1_opt.size()+dir_opt.size()==0) chain1_list.push_back(xname);
    else file2chainlist(chain1_list, xname, dir_opt+dir1_opt, suffix_opt);

    int i; 
    if (dir_opt.size())
        for (i=0;i<chain1_list.size();i++)
            chain2_list.push_back(chain1_list[i]);
    else if (dir2_opt.size()==0) chain2_list.push_back(yname);
    else file2chainlist(chain2_list, yname, dir2_opt, suffix_opt);

    if (outfmt_opt==2)
    {
        if (mm_opt==2) cout<<"#Query\tTemplate\tTM"<<endl;
        else cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
            <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<endl;
    }

    /* real alignment. entry functions are MMalign_main and 
     * TMalign_main */
    if (mm_opt==0) TMalign(xname, yname, fname_super, fname_lign, fname_matrix,
        sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt,
        u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
        split_opt, outfmt_opt, fast_opt, cp_opt, mirror_opt, het_opt,
        atom_opt, mol_opt, dir_opt, dir1_opt, dir2_opt, byresi_opt,
        chain1_list, chain2_list, se_opt);
    else if (mm_opt==1) MMalign(xname, yname, fname_super, fname_lign,
        fname_matrix, sequence, d0_scale, m_opt, o_opt,
        a_opt, d_opt, full_opt, TMcut, infmt1_opt, infmt2_opt,
        ter_opt, split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
        atom_opt, mol_opt, dir1_opt, dir2_opt, chain1_list, chain2_list,
        byresi_opt);
    else if (mm_opt==2) MMdock(xname, yname, fname_super, 
        fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, o_opt, a_opt,
        u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
        split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
        atom_opt, mol_opt, dir1_opt, dir2_opt, chain1_list, chain2_list);
    else if (mm_opt==3) ; // should be changed to mm_opt=0, cp_opt=true
    else if (mm_opt==4) mTMalign(xname, yname, fname_super, fname_matrix,
        sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt,
        u_opt, d_opt, full_opt, TMcut, infmt1_opt, ter_opt,
        split_opt, outfmt_opt, fast_opt, het_opt,
        atom_opt, mol_opt, dir_opt, byresi_opt, chain1_list);
    else if (mm_opt==5 || mm_opt==6) SOIalign(xname, yname, fname_super, fname_lign,
        fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt,
        a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
        split_opt, outfmt_opt, fast_opt, cp_opt, mirror_opt, het_opt,
        atom_opt, mol_opt, dir_opt, dir1_opt, dir2_opt, 
        chain1_list, chain2_list, se_opt, closeK_opt, mm_opt);
    else if (mm_opt==7) flexalign(xname, yname, fname_super, fname_lign, 
        fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt,
        a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
        split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
        atom_opt, mol_opt, dir_opt, dir1_opt, dir2_opt, byresi_opt,
        chain1_list, chain2_list, hinge_opt);
    else cerr<<"WARNING! -mm "<<mm_opt<<" not implemented"<<endl;

    /* clean up */
    vector<string>().swap(chain1_list);
    vector<string>().swap(chain2_list);
    vector<string>().swap(sequence);

    t2 = clock();
    float diff = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    if (outfmt_opt<2) printf("#Total CPU time is %5.2f seconds\n", diff);
    return 0;
}
