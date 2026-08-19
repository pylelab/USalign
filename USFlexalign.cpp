/* command line argument parsing and document of US-align main program */

#include "se.h"
#include "flexalign.h"

using namespace std;

void print_version()
{
    cout << 
"\n"
" ********************************************************************\n"
" * US-Flexalign (Version 20260819)                                  *\n"
" * Flexible Structure Alignment of Proteins and Nucleic Acids       *\n"
" * Reference: Y Zhu, S Yan, Y Zhang, C Zhang (2026)                 *\n"
" * Please email comments and suggestions to cx.zhang2@siat.ac.cn    *\n"
" ********************************************************************"
    << endl;
}

void print_extra_help()
{
    cout << "Additional options:\n"
            "      -v  Print the version of US-alignFlex\n"
            "\n"
            "    -mol  Type of molecule(s) to align.\n"
            "          auto: (default) align both protein and nucleic acids.\n"
            "          prot: only align proteins in a structure.\n"
            "          RNA : only align RNA and DNA in a structure.\n"
            "\n"
            "    -ter  Number of chains to align.\n"
            "          0: align all chains from all models (recommended for aligning\n"
            "             biological assemblies, i.e. biounits)\n"
            "          1: align all chains of the first model (recommended for aligning\n"
            "             asymmetric units)\n"
            "          2: (default) only align the first chain\n"
            "          3: only align the first chain, or the first segment of the\n"
            "             first chain as marked by the 'TER' string in PDB file\n"
            "\n"
            "      -m  Output rotation matrix for superposition, e.g., '-m matrix.txt'\n"
            "          prints the matrix to 'matrix.txt'; '-m -' prints to stdout.\n"
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
            "-chimerax Output superposed structure1 to sup.* for ChimeraX viewing.\n"
            "          $ USalign structure1.pdb structure2.pdb -chimerax sup\n"
            "          $ chimerax --script sup.cxc             # C-alpha trace aligned region\n"
            "          $ chimerax --script sup_all.cxc         # C-alpha trace whole chain\n"
            "          $ chimerax --script sup_atm.cxc         # full-atom aligned region\n"
            "          $ chimerax --script sup_all_atm.cxc     # full-atom whole chain\n"
            "          $ chimerax --script sup_all_atm_lig.cxc # full-atom with all molecules\n"
            "\n"
            "      -a  TM-score normalized by the average length of two structures\n"
            "          T or F, (default F). -a does not change the final alignment.\n"
            "\n"
            "   -fast  Fast but slightly inaccurate alignment\n"
            "\n"
            "    -dir  Perform all-against-all alignment among the list of PDB\n"
            "          chains listed by 'chain_list' under 'chain_folder'.\n"
            "          $ USalign -dir chain_folder/ chain_list\n"
            "\n"
            //"-dirpair  Perform batch alignment for each pair of chains listed by\n"
            //"          'chain_pair_list' under 'chain_folder'. Each line consist of\n"
            //"          two chains, separated by tab or space.\n"
            //"          $ USalign -dirpair chain_folder/ chain_pair_list\n"
            //"\n"
            "   -dir1  Use chain2 to search a list of PDB chains listed by 'chain1_list'\n"
            "          under 'chain1_folder'.\n"
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
            " -infmt1  Input format for structure_1\n"
            " -infmt2  Input format for structure_2\n"
            "          -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
            "           0: PDB format\n"
            "           1: SPICKER format\n"
            //"           2: xyz format\n"
            "           3: PDBx/mmCIF format\n"
            "\n"
            "-chain1   Chains to parse in structure_1\n"
            "-chain2   Chains to parse in structure_2. Use _ for a chain without chain ID.\n"
            "          Multiple chains can be separated by commas, e.g.,\n"
            "          USalign -chain1 C,D,E,F 5jdo.pdb -chain2 A,B,C,D 3wtg.pdb -ter 0\n"
            "\n"
            "-model1   Models to parse in structure_1\n"
            "-model2   Models to parse in structure_2.\n"
            "          Multiple models can be separated by commas, e.g.,\n"
            "          USalign -model1 1,2 1a03.pdb -model2 3,4 1a0n.pdb -ter 0\n"
            "\n"
         << endl;
}

void print_help(bool h_opt = false)
{
    print_version();
    cout << "\n"
            "Usage: USalignFlex PDB1.pdb PDB2.pdb [Options]\n"
            "\n"
            "Options:\n"
            "  -hinge  Maximum number of hinge allowed in flexible alignment.\n"
            "          default: 9\n"
            "\n"
            "    -afp  Enable AFP-enhancement mechanism.\n"
            "\n"
            " -TMpass  Early stopping threshold for AFP-enhancement mechanism.\n"
            "          default: 0.85\n"
            "\n"
            "      -h  Print the full help message, including additional options\n"
            "\n"
            "Example usages ('gunzip' program is needed to read .gz compressed files):\n"
            "    USalignFlex 1oid.pdb 2wde.pdb\n"
         << endl;

    if (h_opt)
        print_extra_help();

    exit(EXIT_SUCCESS);
}

// Unified engine replacing flexalign_greedy and flexalign_usbcat
int flexalign_unified(string &xname, string &yname, const string &fname_super,
                      const string &fname_lign, const string &fname_matrix,
                      vector<string> &sequence, const double Lnorm_ass, const double d0_scale,
                      const bool m_opt, const int i_opt, const int o_opt, const int a_opt,
                      const bool u_opt, const bool d_opt, const double TMcut,
                      const int infmt1_opt, const int infmt2_opt, const int ter_opt,
                      const int split_opt, const int outfmt_opt, const bool fast_opt,
                      const int mirror_opt, const int het_opt, const string &atom_opt,
                      const bool autojustify, const string &mol_opt, const string &dir_opt,
                      const string &dirpair_opt, const string &dir1_opt, const string &dir2_opt,
                      const vector<string> &chain2parse1, const vector<string> &chain2parse2,
                      const vector<string> &model2parse1, const vector<string> &model2parse2,
                      const int byresi_opt, const vector<string> &chain1_list,
                      const vector<string> &chain2_list, const int hinge_opt, const int ss_opt,
                      FlexAlignMode mode = FLEX_STANDARD, bool hinge_set = false, double TMpass = 0.85)
{
    vector<vector<string> > PDB_lines1;
    vector<vector<string> > PDB_lines2;
    vector<int> mol_vec1;
    vector<int> mol_vec2;
    vector<string> chainID_list1;
    vector<string> chainID_list2;
    int i, j, chain_i, chain_j, r, xlen, ylen, xchainnum, ychainnum;
    char *seqx, *seqy, *secx, *secy;
    double **xa, **ya;
    vector<string> resi_vec1;
    vector<string> resi_vec2;
    int read_resi = byresi_opt;
    if (byresi_opt == 0 && o_opt)
        read_resi = 2;

    for (i = 0; i < chain1_list.size(); i++)
    {
        xname = chain1_list[i];
        xchainnum = get_PDB_lines(xname, PDB_lines1, chainID_list1,
                                  mol_vec1, ter_opt, infmt1_opt, atom_opt, autojustify,
                                  split_opt, het_opt, chain2parse1, model2parse1);
        if (!xchainnum)
        {
            cerr << "Warning! Cannot parse file: " << xname << ". Chain number 0." << endl;
            continue;
        }
        for (chain_i = 0; chain_i < xchainnum; chain_i++)
        {
            xlen = PDB_lines1[chain_i].size();
            if (mol_opt == "RNA")
                mol_vec1[chain_i] = 1;
            else if (mol_opt == "protein")
                mol_vec1[chain_i] = -1;
            if (xlen < 3)
                continue;

            NewArray(&xa, xlen, 3);
            seqx = new char[xlen + 1];
            secx = new char[xlen + 1];
            read_PDB(PDB_lines1[chain_i], xa, seqx, resi_vec1, read_resi);
            if (mirror_opt)
                for (r = 0; r < xlen; r++)
                    xa[r][2] = -xa[r][2];
            (mol_vec1[chain_i] > 0) ? make_sec(seqx, xa, xlen, secx, atom_opt) : make_sec(xa, xlen, secx);

            for (j = (dir_opt.size() > 0) * (i + 1); j < chain2_list.size(); j++)
            {
                if (dirpair_opt.size() && i != j)
                    continue;
                if (PDB_lines2.size() == 0)
                {
                    yname = chain2_list[j];
                    ychainnum = get_PDB_lines(yname, PDB_lines2, chainID_list2,
                                              mol_vec2, ter_opt, infmt2_opt, atom_opt, autojustify,
                                              split_opt, het_opt, chain2parse2, model2parse2);
                    if (!ychainnum)
                        continue;
                }
                for (chain_j = 0; chain_j < ychainnum; chain_j++)
                {
                    ylen = PDB_lines2[chain_j].size();
                    if (mol_opt == "RNA")
                        mol_vec2[chain_j] = 1;
                    else if (mol_opt == "protein")
                        mol_vec2[chain_j] = -1;
                    if (ylen < 3)
                        continue;

                    NewArray(&ya, ylen, 3);
                    seqy = new char[ylen + 1];
                    secy = new char[ylen + 1];
                    read_PDB(PDB_lines2[chain_j], ya, seqy, resi_vec2, read_resi);
                    (mol_vec2[chain_j] > 0) ? make_sec(seqy, ya, ylen, secy, atom_opt) : make_sec(ya, ylen, secy);

                    if (byresi_opt)
                        extract_aln_from_resi(sequence, seqx, seqy, resi_vec1, resi_vec2, byresi_opt);

                    // --- CORE DISPATCH LOGIC START ---
                    if (mode == FLEX_USBCAT)
                    {
                        FlexAlignResult usbcat_res;
                        bool force_fast_opt = (getmin(xlen, ylen) > 1500) ? true : fast_opt;

                        usbcat_res.hingeNum = flexalign_usbcat_main(
                            xa, ya, seqx, seqy, secx, secy,
                            usbcat_res.t0, usbcat_res.u0, usbcat_res.tu_vec,
                            usbcat_res.TM1, usbcat_res.TM2, usbcat_res.TM3, usbcat_res.TM4, usbcat_res.TM5,
                            usbcat_res.d0_0, usbcat_res.TM_0,
                            usbcat_res.d0A, usbcat_res.d0B, usbcat_res.d0u, usbcat_res.d0a, usbcat_res.d0_out,
                            usbcat_res.seqM, usbcat_res.seqxA, usbcat_res.seqyA, usbcat_res.do_vec,
                            usbcat_res.rmsd0, usbcat_res.L_ali, usbcat_res.Liden,
                            usbcat_res.TM_ali, usbcat_res.rmsd_ali, usbcat_res.n_ali, usbcat_res.n_ali8,
                            xlen, ylen, sequence, Lnorm_ass, d0_scale,
                            i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                            mol_vec1[chain_i] + mol_vec2[chain_j], hinge_opt, ss_opt, 0, hinge_set, TMpass);

                        if (outfmt_opt == 0)
                            print_version();
                        output_flexalign_results(
                            xname.substr(dir1_opt.size() + dir_opt.size() + dirpair_opt.size()),
                            yname.substr(dir2_opt.size() + dir_opt.size() + dirpair_opt.size()),
                            chainID_list1[chain_i], chainID_list2[chain_j],
                            xlen, ylen, usbcat_res.t0, usbcat_res.u0, usbcat_res.tu_vec,
                            usbcat_res.TM1, usbcat_res.TM2, usbcat_res.TM3, usbcat_res.TM4, usbcat_res.TM5,
                            usbcat_res.rmsd0, usbcat_res.d0_out, usbcat_res.seqM.c_str(),
                            usbcat_res.seqxA.c_str(), usbcat_res.seqyA.c_str(), usbcat_res.Liden,
                            usbcat_res.n_ali8, usbcat_res.L_ali, usbcat_res.TM_ali, usbcat_res.rmsd_ali,
                            usbcat_res.TM_0, usbcat_res.d0_0,
                            usbcat_res.d0A, usbcat_res.d0B, Lnorm_ass, d0_scale, usbcat_res.d0a, usbcat_res.d0u,
                            (m_opt ? fname_matrix : "").c_str(),
                            outfmt_opt, ter_opt, false, split_opt, o_opt,
                            fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
                            resi_vec1, resi_vec2);
                    }
                    else
                    {
                        // === Standard & Best specific logic ===
                        FlexAlignResult best_res;
                        double global_max_TM = -1.0;

                        int start_ss = (mode == FLEX_BEST) ? 0 : ss_opt;
                        int end_ss = (mode == FLEX_BEST) ? 1 : ss_opt;

                        bool force_fast_opt = (getmin(xlen, ylen) > 1500) ? true : fast_opt;

                        for (int cur_ss_opt = start_ss; cur_ss_opt <= end_ss; cur_ss_opt++)
                        {
                            FlexAlignResult cur_res;
                            execute_flexalign_with_fallback(
                                xa, ya, seqx, seqy, secx, secy, xlen, ylen, sequence, Lnorm_ass, d0_scale,
                                i_opt, a_opt, u_opt, d_opt, force_fast_opt, mol_vec1[chain_i] + mol_vec2[chain_j],
                                hinge_opt, cur_ss_opt, cur_res);

                            double cur_max_TM = (cur_res.TM1 > cur_res.TM2) ? cur_res.TM1 : cur_res.TM2;
                            if (cur_max_TM > global_max_TM)
                            {
                                global_max_TM = cur_max_TM;
                                best_res = cur_res;
                            }
                        }

                        if (outfmt_opt == 0)
                            print_version();
                        output_flexalign_results(
                            xname.substr(dir1_opt.size() + dir_opt.size() + dirpair_opt.size()),
                            yname.substr(dir2_opt.size() + dir_opt.size() + dirpair_opt.size()),
                            chainID_list1[chain_i], chainID_list2[chain_j],
                            xlen, ylen, best_res.t0, best_res.u0, best_res.tu_vec, best_res.TM1, best_res.TM2, best_res.TM3, best_res.TM4, best_res.TM5,
                            best_res.rmsd0, best_res.d0_out, best_res.seqM.c_str(),
                            best_res.seqxA.c_str(), best_res.seqyA.c_str(), best_res.Liden,
                            best_res.n_ali8, best_res.L_ali, best_res.TM_ali, best_res.rmsd_ali, best_res.TM_0, best_res.d0_0,
                            best_res.d0A, best_res.d0B, Lnorm_ass, d0_scale, best_res.d0a, best_res.d0u,
                            (m_opt ? fname_matrix : "").c_str(),
                            outfmt_opt, ter_opt, false, split_opt, o_opt,
                            fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
                            resi_vec1, resi_vec2);
                    }
                    // --- CORE DISPATCH LOGIC END ---

                    // Cleanup memory
                    DeleteArray(&ya, ylen);
                    delete[] seqy;
                    delete[] secy;
                    resi_vec2.clear();
                }
                if (chain2_list.size() > 1)
                {
                    yname.clear();
                    for (chain_j = 0; chain_j < ychainnum; chain_j++)
                        PDB_lines2[chain_j].clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            }
            PDB_lines1[chain_i].clear();
            DeleteArray(&xa, xlen);
            delete[] seqx;
            delete[] secx;
            resi_vec1.clear();
        }
        xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    }
    if (chain2_list.size() == 1)
    {
        yname.clear();
        for (chain_j = 0; chain_j < ychainnum; chain_j++)
            PDB_lines2[chain_j].clear();
        PDB_lines2.clear();
        resi_vec2.clear();
        chainID_list2.clear();
        mol_vec2.clear();
    }
    return 0;
}

// =======================================================================
// Direct Drop-in Wrappers (No changes needed in main() bindings)
// =======================================================================

int flexalign_greedy(string &xname, string &yname, const string &fname_super, const string &fname_lign, const string &fname_matrix, vector<string> &sequence, const double Lnorm_ass, const double d0_scale, const bool m_opt, const int i_opt, const int o_opt, const int a_opt, const bool u_opt, const bool d_opt, const double TMcut, const int infmt1_opt, const int infmt2_opt, const int ter_opt, const int split_opt, const int outfmt_opt, const bool fast_opt, const int mirror_opt, const int het_opt, const string &atom_opt, const bool autojustify, const string &mol_opt, const string &dir_opt, const string &dirpair_opt, const string &dir1_opt, const string &dir2_opt, const vector<string> &chain2parse1, const vector<string> &chain2parse2, const vector<string> &model2parse1, const vector<string> &model2parse2, const int byresi_opt, const vector<string> &chain1_list, const vector<string> &chain2_list, const int hinge_opt)
{
    return flexalign_unified(xname, yname, fname_super, fname_lign, fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt, split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt, atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt, dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2, byresi_opt, chain1_list, chain2_list, hinge_opt, 0 /* ss_opt is ignored in BEST mode */, FLEX_BEST);
}

int flexalign_usbcat(string &xname, string &yname, const string &fname_super, const string &fname_lign, const string &fname_matrix, vector<string> &sequence, const double Lnorm_ass, const double d0_scale, const bool m_opt, const int i_opt, const int o_opt, const int a_opt, const bool u_opt, const bool d_opt, const double TMcut, const int infmt1_opt, const int infmt2_opt, const int ter_opt, const int split_opt, const int outfmt_opt, const bool fast_opt, const int mirror_opt, const int het_opt, const string &atom_opt, const bool autojustify, const string &mol_opt, const string &dir_opt, const string &dirpair_opt, const string &dir1_opt, const string &dir2_opt, const vector<string> &chain2parse1, const vector<string> &chain2parse2, const vector<string> &model2parse1, const vector<string> &model2parse2, const int byresi_opt, const vector<string> &chain1_list, const vector<string> &chain2_list, const int hinge_opt, bool hinge_set = false, double TMpass = 0.85)
{
    return flexalign_unified(xname, yname, fname_super, fname_lign, fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt, split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt, atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt, dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2, byresi_opt, chain1_list, chain2_list, hinge_opt, 0 /* ss_opt ignore */, FLEX_USBCAT, hinge_set, TMpass);
}

int main(int argc, char *argv[])
{
    if (argc < 2)
        print_help();

    clock_t t1, t2;
    t1 = clock();

    /**********************/
    /*    get argument    */
    /**********************/
    string xname = "";
    string yname = "";
    string fname_super = "";  // file name for superposed structure
    string fname_lign = "";   // file name for user alignment
    string fname_matrix = ""; // file name for output matrix
    vector<string> sequence;  // get value from alignment file
    double Lnorm_ass, d0_scale;

    bool h_opt = false;  // print full help message
    bool v_opt = false;  // print version
    bool m_opt = false;  // flag for -m, output rotation matrix
    int i_opt = 0;       // 1 for -i, 3 for -I
    int o_opt = 0;       // 1 for -o, 2 for -rasmol, 3 for -chimerax
    int a_opt = 0;       // flag for -a, do not normalized by average length
    bool u_opt = false;  // flag for -u, normalized by user specified length
    bool d_opt = false;  // flag for -d, user specified d0

    bool full_opt = false; // do not show chain level alignment
    double TMcut = -1;
    int infmt1_opt = -1;   // PDB or PDBx/mmCIF format for chain_1
    int infmt2_opt = -1;   // PDB or PDBx/mmCIF format for chain_2
    int ter_opt = -1;      // default change to 2 (END, or different chainID)
    int split_opt = -1;    // default change to 2 (split each chains)
    int outfmt_opt = 0;    // set -outfmt to full output
    bool fast_opt = false; // flags for -fast, fTM-align algorithm
    int hinge_opt = 9;     // maximum number of hinge allowed for flexible
    bool hinge_set = false;
    double TMpass_opt = 0.85;
    int mirror_opt = 0;         // do not align mirror
    int het_opt = 0;            // do not read HETATM residues
    int mm_opt = 0;             // do not perform MM-align
    bool usbcat_opt = false;    // flag for -afp, only valid with -mm 7
    string atom_opt = "auto";   // use C alpha atom for protein and C3' for RNA
    string mol_opt = "auto";    // auto-detect the molecule type as protein/RNA
    string suffix_opt = "";     // set -suffix to empty
    string dir_opt = "";        // set -dir to empty
    string dirpair_opt = "";    // set -dirpair to empty
    string dir1_opt = "";       // set -dir1 to empty
    string dir2_opt = "";       // set -dir2 to empty
    int byresi_opt = 0;         // set -byresi to 0
    vector<string> chain1_list; // only when -dir1 is set
    vector<string> chain2_list; // only when -dir2 is set
    vector<string> chain2parse1;
    vector<string> chain2parse2;
    vector<string> model2parse1;
    vector<string> model2parse2;
    vector<pair<string, string> > chain_pair_list; // only when -dirpair is set

    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-o"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -o");
            if (o_opt == 2)
                cerr << "Warning! -rasmol is already set. Ignore -o" << endl;
            else if (o_opt == 3)
                cerr << "Warning! -chimerax is already set. Ignore -o" << endl;
            else
            {
                fname_super = argv[i + 1];
                o_opt = 1;
            }
            i++;
        }
        else if (!strcmp(argv[i], "-rasmol"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -rasmol");
            if (o_opt == 1)
                cerr << "Warning! -o is already set. Ignore -rasmol" << endl;
            else if (o_opt == 3)
                cerr << "Warning! -chimerax is already set. Ignore -rasmol" << endl;
            else
            {
                fname_super = argv[i + 1];
                o_opt = 2;
            }
            i++;
        }
        else if (!strcmp(argv[i], "-chimerax"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -chimerax");
            if (o_opt == 1)
                cerr << "Warning! -o is already set. Ignore -chimerax" << endl;
            else if (o_opt == 2)
                cerr << "Warning! -rasmol is already set. Ignore -chimerax" << endl;
            else
            {
                fname_super = argv[i + 1];
                o_opt = 3;
            }
            i++;
        }
        else if (!strcmp(argv[i], "-u") || !strcmp(argv[i], "-L"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -u or -L");
            Lnorm_ass = atof(argv[i + 1]);
            u_opt = true;
            i++;
            if (Lnorm_ass <= 0)
                PrintErrorAndQuit(
                    "ERROR! The value for -u or -L should be >0");
        }
        else if (!strcmp(argv[i], "-a"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -a");
            if (!strcmp(argv[i + 1], "T"))
                a_opt = true;
            else if (!strcmp(argv[i + 1], "F"))
                a_opt = false;
            else
            {
                a_opt = atoi(argv[i + 1]);
                if (a_opt != -2 && a_opt != -1 && a_opt != 1)
                    PrintErrorAndQuit("-a must be -2, -1, 1, T or F");
            }
            i++;
        }
        else if (!strcmp(argv[i], "-d"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -d");
            d0_scale = atof(argv[i + 1]);
            d_opt = true;
            i++;
        }
        else if (!strcmp(argv[i], "-hinge"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -hinge");
            hinge_set = true;
            hinge_opt = atoi(argv[i + 1]);
            i++;
        }
        else if (!strcmp(argv[i], "-v"))
        {
            v_opt = true;
        }
        else if (!strcmp(argv[i], "-h"))
        {
            h_opt = true;
        }
        else if (!strcmp(argv[i], "-chain1"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -chain1");
            split(argv[i + 1], chain2parse1, ',');
            i++;
        }
        else if (!strcmp(argv[i], "-chain2"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -chain2");
            split(argv[i + 1], chain2parse2, ',');
            i++;
        }
        else if (!strcmp(argv[i], "-model1"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -model1");
            split(argv[i + 1], model2parse1, ',');
            i++;
        }
        else if (!strcmp(argv[i], "-model2"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -model2");
            split(argv[i + 1], model2parse2, ',');
            i++;
        }
        else if (!strcmp(argv[i], "-m"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -m");
            fname_matrix = argv[i + 1];
            m_opt = true;
            i++;
        } // get filename for rotation matrix
        else if (!strcmp(argv[i], "-fast"))
        {
            fast_opt = true;
        }
        else if (!strcmp(argv[i], "-infmt1"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -infmt1");
            infmt1_opt = atoi(argv[i + 1]);
            i++;
            if (infmt1_opt < -1 || infmt1_opt > 3)
                PrintErrorAndQuit("ERROR! -infmt1 can only be -1, 0, 1, 2, or 3");
        }
        else if (!strcmp(argv[i], "-infmt2"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -infmt2");
            infmt2_opt = atoi(argv[i + 1]);
            i++;
            if (infmt2_opt < -1 || infmt2_opt > 3)
                PrintErrorAndQuit("ERROR! -infmt2 can only be -1, 0, 1, 2, or 3");
        }
        else if (!strcmp(argv[i], "-ter"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -ter");
            ter_opt = atoi(argv[i + 1]);
            i++;
        }
        else if (!strcmp(argv[i], "-split"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -split");
            split_opt = atoi(argv[i + 1]);
            i++;
        }
        else if (!strcmp(argv[i], "-atom"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -atom");
            atom_opt = argv[i + 1];
            i++;
        }
        else if (!strcmp(argv[i], "-mol"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -mol");
            mol_opt = argv[i + 1];
            i++;
            if (mol_opt == "prot")
                mol_opt = "protein";
            else if (mol_opt == "DNA")
                mol_opt = "RNA";
            if (mol_opt != "auto" && mol_opt != "protein" && mol_opt != "RNA")
                PrintErrorAndQuit("ERROR! Molecule type must be one of the "
                                  "following:\nauto, prot (the same as 'protein'), and "
                                  "RNA (the same as 'DNA').");
        }
        else if (!strcmp(argv[i], "-dir"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -dir");
            dir_opt = argv[i + 1];
            i++;
        }
        else if (!strcmp(argv[i], "-dirpair"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -dirpair");
            dirpair_opt = argv[i + 1];
            i++;
        }
        else if (!strcmp(argv[i], "-dir1"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -dir1");
            dir1_opt = argv[i + 1];
            i++;
        }
        else if (!strcmp(argv[i], "-dir2"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -dir2");
            dir2_opt = argv[i + 1];
            i++;
        }
        else if (!strcmp(argv[i], "-suffix"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -suffix");
            suffix_opt = argv[i + 1];
            i++;
        }
        else if (!strcmp(argv[i], "-outfmt"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -outfmt");
            outfmt_opt = atoi(argv[i + 1]);
            i++;
        }
        else if (!strcmp(argv[i], "-mirror"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -mirror");
            mirror_opt = atoi(argv[i + 1]);
            i++;
        }
        else if (!strcmp(argv[i], "-het"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -het");
            het_opt = atoi(argv[i + 1]);
            i++;
            if (het_opt != 0 && het_opt != 1 && het_opt != 2)
                PrintErrorAndQuit("-het must be 0, 1, or 2");
        }
        else if (!strcmp(argv[i], "-afp"))
        {
            usbcat_opt = true;
        }
        else if (!strcmp(argv[i], "-TMpass")) // Parse the -TMpass argument
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -TMpass");
            TMpass_opt = atof(argv[i + 1]);
            i++;
        }
        else if (xname.size() == 0)
            xname = argv[i];
        else if (yname.size() == 0)
            yname = argv[i];
        else
            PrintErrorAndQuit(string("ERROR! Undefined option ") + argv[i]);
    }

    if (xname.size() == 0 || (yname.size() && dir_opt.size()) ||
        (yname.size() && dirpair_opt.size()) ||
        (yname.size() == 0 && dir_opt.size() == 0 && dirpair_opt.size() == 0))
    {
        if (h_opt)
            print_help(h_opt);
        if (v_opt)
        {
            print_version();
            exit(EXIT_FAILURE);
        }
        if (xname.size() == 0)
            PrintErrorAndQuit("Please provide input structures");
        else if (yname.size() == 0 && dir_opt.size() == 0 && dirpair_opt.size() == 0 && mm_opt != 4)
            PrintErrorAndQuit("Please provide structure B");
        else if (yname.size() && dir_opt.size() + dirpair_opt.size())
            PrintErrorAndQuit("Please provide only one file name if -dir is set");
    }

    if (suffix_opt.size() && dir_opt.size() + dirpair_opt.size() + dir1_opt.size() + dir2_opt.size() == 0)
        PrintErrorAndQuit("-suffix is only valid if -dir, -dir1 or -dir2 is set");
    if ((dir_opt.size() || dirpair_opt.size() || dir1_opt.size() || dir2_opt.size()))
    {
        if (mm_opt != 2 && mm_opt != 4)
        {
            if (o_opt)
                PrintErrorAndQuit("-o cannot be set with -dir, -dir1 or -dir2");
            if (m_opt && fname_matrix != "-")
                PrintErrorAndQuit("-m can only be - or unset when using -dir, -dir1 or -dir2");
        }
        else if ((dir_opt.size() || dirpair_opt.size()) && (dir1_opt.size() || dir2_opt.size()))
            PrintErrorAndQuit("-dir cannot be set with -dir1 or -dir2");
        else if (dir_opt.size() && dirpair_opt.size())
            PrintErrorAndQuit("-dir cannot be set with -dirpair");
    }
    if (o_opt && (infmt1_opt != -1 && infmt1_opt != 0 && infmt1_opt != 3))
        PrintErrorAndQuit("-o can only be used with -infmt1 -1, 0 or 3");

    bool autojustify = (atom_opt == "auto" || atom_opt == "PC4'"); // auto re-pad atom name
    if (mol_opt == "protein" && atom_opt == "auto")
        atom_opt = " CA ";
    else if (mol_opt == "RNA" && atom_opt == "auto")
        atom_opt = " C3'";
    if (atom_opt.size() != 4)
    {
        cerr << "ERROR! Atom name must have 4 characters, including space.\n"
                "For example, C alpha, C3' and P atoms should be specified by\n"
                "-atom \" CA \", -atom \" P  \" and -atom \" C3'\", respectively."
             << endl;
        if (atom_opt.size() >= 5 || atom_opt.size() == 0)
            return 1;
        else if (atom_opt.size() == 1)
            atom_opt = " " + atom_opt + "  ";
        else if (atom_opt.size() == 2)
            atom_opt = " " + atom_opt + " ";
        else if (atom_opt.size() == 3)
            atom_opt = " " + atom_opt;
        cerr << "Change -atom to \"" << atom_opt << "\"" << endl;
    }

    if (d_opt && d0_scale <= 0)
        PrintErrorAndQuit("Wrong value for option -d! It should be >0");
    if (outfmt_opt >= 2 && (a_opt || u_opt || d_opt))
        PrintErrorAndQuit("-outfmt 2 cannot be used with -a, -u, -L, -d");
    if (split_opt < 0)
            split_opt = 2;
    else if (split_opt > 2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");

    if (mirror_opt && het_opt != 1)
        cerr << "WARNING! -mirror was not used with -het 1. "
             << "D amino acids may not be correctly aligned." << endl;

    if (ter_opt < 0)
            ter_opt = 2;

    if (o_opt==2 && ter_opt <= 1 && split_opt == 2)
    {
            cerr << "WARNING! -mm may generate incorrect"
                 << " RasMol output due to limitations in PDB file format. "
                 << "When -mm is used, -o is recommended over -rasmol" << endl;
    }

    if (hinge_opt >= 10)
        PrintErrorAndQuit("ERROR! -hinge must be <10");

    /* read initial alignment file from 'align.txt' */
    if (i_opt)
        read_user_alignment(sequence, fname_lign, i_opt);

    if (m_opt && fname_matrix == "") // Output rotation matrix: matrix.txt
        PrintErrorAndQuit("ERROR! Please provide a file name for option -m!");

    /* parse file list */
    int i;
    if (dirpair_opt.size())
        file2chainpairlist(chain1_list, chain2_list, xname, dirpair_opt, suffix_opt);
    else
    {
        if (dir1_opt.size() + dir_opt.size() == 0)
            chain1_list.push_back(xname);
        else
            file2chainlist(chain1_list, xname, dir_opt + dir1_opt, suffix_opt);

        if (dir_opt.size())
            for (i = 0; i < chain1_list.size(); i++)
                chain2_list.push_back(chain1_list[i]);
        else if (dir2_opt.size() == 0)
            chain2_list.push_back(yname);
        else
            file2chainlist(chain2_list, yname, dir2_opt, suffix_opt);
    }

    if (outfmt_opt == 2)
    {
            cout << "#PDBchain1\tPDBchain2\tTM1\tTM2\t"
                 << "RMSD\tID1\tID2\tIDali\tL1\tL2\tLali\tNblk" << endl;
    }

    /* real alignment. entry functions are MMalign_main and
     * TMalign_main */
    {
        if (usbcat_opt)
            flexalign_usbcat(xname, yname, fname_super, fname_lign,
                             fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt,
                             a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
                             split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
                             atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt,
                             dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2,
                             byresi_opt, chain1_list, chain2_list, hinge_opt, hinge_set, TMpass_opt);
        else
            flexalign_greedy(xname, yname, fname_super, fname_lign,
                             fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt,
                             a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
                             split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
                             atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt,
                             dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2,
                             byresi_opt, chain1_list, chain2_list, hinge_opt);
    }

    /* clean up */
    vector<string>().swap(chain1_list);
    vector<string>().swap(chain2_list);
    vector<string>().swap(chain2parse1);
    vector<string>().swap(chain2parse2);
    vector<string>().swap(model2parse1);
    vector<string>().swap(model2parse2);
    vector<string>().swap(sequence);
    vector<pair<string, string> >().swap(chain_pair_list);

    t2 = clock();
    float diff = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;
    if (outfmt_opt < 2)
        printf("#Total CPU time is %5.2f seconds\n", diff);
    return 0;
}
