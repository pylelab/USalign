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
" * US-align (Version 20260527)                                      *\n"
" * Universal Structure Alignment of Proteins and Nucleic Acids      *\n"
" * Reference: C Zhang, L Freddolino, Y Zhang. (2026) Nat Protoc     *\n"
" *            C Zhang, M Shine, AM Pyle, Y Zhang. (2022) Nat Methods*\n"
" *            C Zhang, AM Pyle (2022) iScience.                     *\n"
" * Please email comments and suggestions to zhang@zhanggroup.org    *\n"
" ********************************************************************"
    << endl;
}

void print_extra_help()
{
    cout << "Additional options:\n"
            "      -v  Print the version of US-align\n"
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
            "              (default for -TMscore 2)\n"
            "           1: treat each MODEL as a separate chain\n"
            "           2: (default for other cases) treat each chain as a separate chain\n"
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
            "\n"
            " -hinge   Maximum number of hinge allowed in flexible alignment. default: 9\n"
            "\n"
            "   -se    Do not perform superposition. Useful for extracting alignment from\n"
            "          superposed structure pairs\n"
            "\n"
            " -infmt1  Input format for structure_1\n"
            " -infmt2  Input format for structure_2\n"
            "          -1: (default) automatically detect PDB or PDBx/mmCIF format\n"
            "           0: PDB format\n"
            "           1: SPICKER format\n"
            //"           2: xyz format\n"
            "           3: PDBx/mmCIF format\n"
            "\n"
            "-chainmap (only useful for -mm 1) use the final chain mapping 'chainmap.txt'\n"
            "          specified by user. 'chainmap.txt' is a tab-seperated text with two\n"
            "          columns, one for each complex\n"
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
            "Advanced usage 1 (generate an image for a pair of superposed structures):\n"
            "    USalign 1cpc.pdb 1mba.pdb -o sup\n"
            "    pymol -c -d @sup_all_atm.pml -g sup_all_atm.png\n"
            "\n"
            "Advanced usage 2 (a quick search of query.pdb against I-TASSER PDB library):\n"
            "    wget https://zhanggroup.org/library/PDB.tar.bz2\n"
            "    tar -xjvf PDB.tar.bz2\n"
            "    USalign query.pdb -dir2 PDB/ PDB/list -suffix .pdb -outfmt 2 -fast\n"
         << endl;
}

void print_help(bool h_opt = false)
{
    print_version();
    cout << "\n"
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
            "          4: MSTA, i.e., alignment of multiple monomeric chains into a\n"
            "             consensus alignment\n"
            "             $ USalign -dir chains/ list -suffix .pdb -mm 4\n"
            "          5: fully non-sequential (fNS) alignment\n"
            "          6: semi-non-sequential (sNS) alignment\n"
            "          To use -mm 1 or -mm 2, '-ter' option must be 0 or 1.\n"
            "\n"
            "    -ter  Number of chains to align.\n"
            "          0: align all chains from all models (recommended for aligning\n"
            "             biological assemblies, i.e. biounits)\n"
            "          1: align all chains of the first model (recommended for aligning\n"
            "             asymmetric units, default for -mm 1,2 and -TMscore 2,6,7)\n"
            "          2: (default for other cases) only align the first chain\n"
            "          3: only align the first chain, or the first segment of the\n"
            "             first chain as marked by the 'TER' string in PDB file\n"
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
            "          7: sequence dependent alignment of two complex structures:\n"
            "             perform global sequence alignment of each chain pair, derive\n"
            "             optimal chain mapping, and then superpose two complex\n"
            "             structures by TM-score\n"
            "\n"
            "      -I  Use the final alignment specified by FASTA file 'align.txt'\n"
            "\n"
            "      -i  Use alignment specified by 'align.txt' as an initial alignment\n"
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
            "     -do  Output distance of aligned residue pairs\n"
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
         << endl;

    // if (h_opt)
    print_extra_help();

    exit(EXIT_SUCCESS);
}

/* TMalign, RNAalign, CPalign, TMscore */
int TMalign(string &xname, string &yname, const string &fname_super,
            const string &fname_lign, const string &fname_matrix,
            vector<string> &sequence, const double Lnorm_ass, const double d0_scale,
            const bool m_opt, const int i_opt, const int o_opt, const int a_opt,
            const bool u_opt, const bool d_opt, const double TMcut,
            const int infmt1_opt, const int infmt2_opt, const int ter_opt,
            const int split_opt, const int outfmt_opt, const bool fast_opt,
            const int cp_opt, const int mirror_opt, const int het_opt,
            const string &atom_opt, const bool autojustify, const string &mol_opt,
            const string &dir_opt, const string &dirpair_opt, const string &dir1_opt,
            const string &dir2_opt, const vector<string> &chain2parse1,
            const vector<string> &chain2parse2, const vector<string> &model2parse1,
            const vector<string> &model2parse2, const int byresi_opt,
            const vector<string> &chain1_list, const vector<string> &chain2_list,
            const bool se_opt, const bool do_opt)
{
    /* declare previously global variables */
    vector<vector<string>> PDB_lines1; // text of chain1
    vector<vector<string>> PDB_lines2; // text of chain2
    vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    vector<string> chainID_list1;      // list of chainID1
    vector<string> chainID_list2;      // list of chainID2
    int i, j;                          // file index
    int chain_i, chain_j;              // chain index
    int r;                             // residue index
    int xlen, ylen;                    // chain length
    int xchainnum, ychainnum;          // number of chains in a PDB file
    char *seqx, *seqy;                 // for the protein sequence
    char *secx, *secy;                 // for the secondary structure
    double **xa, **ya;                 // for input vectors xa[0...xlen-1][0..2] and
                                       // ya[0...ylen-1][0..2], in general,
                                       // ya is regarded as native structure
                                       // --> superpose xa onto ya
    vector<string> resi_vec1;          // residue index for chain1
    vector<string> resi_vec2;          // residue index for chain2
    int read_resi = byresi_opt;        // whether to read residue index
    if (byresi_opt == 0 && o_opt)
        read_resi = 2;

    /* loop over file names */
    for (i = 0; i < chain1_list.size(); i++)
    {
        /* parse chain 1 */
        xname = chain1_list[i];
        xchainnum = get_PDB_lines(xname, PDB_lines1, chainID_list1, mol_vec1,
                                  ter_opt, infmt1_opt, atom_opt, autojustify, split_opt, het_opt,
                                  chain2parse1, model2parse1);
        if (!xchainnum)
        {
            cerr << "Warning! Cannot parse file: " << xname
                 << ". Chain number 0." << endl;
            continue;
        }
        for (chain_i = 0; chain_i < xchainnum; chain_i++)
        {
            xlen = PDB_lines1[chain_i].size();
            if (mol_opt == "RNA")
                mol_vec1[chain_i] = 1;
            else if (mol_opt == "protein")
                mol_vec1[chain_i] = -1;
            if (!xlen)
            {
                cerr << "Warning! Cannot parse file: " << xname
                     << ". Chain length 0." << endl;
                continue;
            }
            else if (xlen < 3)
            {
                cerr << "Sequence is too short <3!: " << xname << endl;
                continue;
            }
            NewArray(&xa, xlen, 3);
            seqx = new char[xlen + 1];
            secx = new char[xlen + 1];
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx,
                            resi_vec1, read_resi);
            if (mirror_opt)
                for (r = 0; r < xlen; r++)
                    xa[r][2] = -xa[r][2];
            if (mol_vec1[chain_i] > 0)
                make_sec(seqx, xa, xlen, secx, atom_opt);
            else
                make_sec(xa, xlen, secx); // secondary structure assignment

            for (j = (dir_opt.size() > 0) * (i + 1); j < chain2_list.size(); j++)
            {
                if (dirpair_opt.size() && j != i)
                    continue;
                /* parse chain 2 */
                if (PDB_lines2.size() == 0)
                {
                    yname = chain2_list[j];
                    ychainnum = get_PDB_lines(yname, PDB_lines2, chainID_list2,
                                              mol_vec2, ter_opt, infmt2_opt, atom_opt, autojustify,
                                              split_opt, het_opt, chain2parse2, model2parse2);
                    if (!ychainnum)
                    {
                        cerr << "Warning! Cannot parse file: " << yname
                             << ". Chain number 0." << endl;
                        continue;
                    }
                }
                for (chain_j = 0; chain_j < ychainnum; chain_j++)
                {
                    ylen = PDB_lines2[chain_j].size();
                    if (mol_opt == "RNA")
                        mol_vec2[chain_j] = 1;
                    else if (mol_opt == "protein")
                        mol_vec2[chain_j] = -1;
                    if (!ylen)
                    {
                        cerr << "Warning! Cannot parse file: " << yname
                             << ". Chain length 0." << endl;
                        continue;
                    }
                    else if (ylen < 3)
                    {
                        cerr << "Sequence is too short <3!: " << yname << endl;
                        continue;
                    }
                    NewArray(&ya, ylen, 3);
                    seqy = new char[ylen + 1];
                    secy = new char[ylen + 1];
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy,
                                    resi_vec2, read_resi);
                    if (mol_vec2[chain_j] > 0)
                        make_sec(seqy, ya, ylen, secy, atom_opt);
                    else
                        make_sec(ya, ylen, secy);

                    if (byresi_opt)
                        extract_aln_from_resi(sequence,
                                              seqx, seqy, resi_vec1, resi_vec2, byresi_opt);

                    /* declare variable specific to this pair of TMalign */
                    double t0[3], u0[3][3];
                    double TM1, TM2;
                    double TM3, TM4, TM5; // for a_opt, u_opt, d_opt
                    double d0_0, TM_0;
                    double d0A, d0B, d0u, d0a;
                    double d0_out = 5.0;
                    string seqM, seqxA, seqyA; // for output alignment
                    double rmsd0 = 0.0;
                    int L_ali; // Aligned length in standard_TMscore
                    double Liden = 0;
                    double TM_ali, rmsd_ali; // TMscore and rmsd in standard_TMscore
                    int n_ali = 0;
                    int n_ali8 = 0;
                    bool force_fast_opt = (getmin(xlen, ylen) > 1500) ? true : fast_opt;
                    vector<double> do_vec;

                    /* entry function for structure alignment */
                    if (cp_opt)
                        CPalign_main(
                            xa, ya, seqx, seqy, secx, secy,
                            t0, u0, TM1, TM2, TM3, TM4, TM5,
                            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                            seqM, seqxA, seqyA, do_vec,
                            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                            xlen, ylen, sequence, Lnorm_ass, d0_scale,
                            i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                            mol_vec1[chain_i] + mol_vec2[chain_j], TMcut);
                    else if (se_opt)
                    {
                        int *invmap = new int[ylen + 1];
                        u0[0][0] = u0[1][1] = u0[2][2] = 1;
                        u0[0][1] = u0[0][2] =
                            u0[1][0] = u0[1][2] =
                                u0[2][0] = u0[2][1] =
                                    t0[0] = t0[1] = t0[2] = 0;
                        se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                                seqM, seqxA, seqyA, do_vec,
                                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                                xlen, ylen, sequence, Lnorm_ass, d0_scale,
                                i_opt, a_opt, u_opt, d_opt,
                                mol_vec1[chain_i] + mol_vec2[chain_j],
                                outfmt_opt, invmap);
                        if (outfmt_opt >= 2)
                        {
                            Liden = L_ali = 0;
                            int r1, r2;
                            for (r2 = 0; r2 < ylen; r2++)
                            {
                                r1 = invmap[r2];
                                if (r1 < 0)
                                    continue;
                                L_ali += 1;
                                Liden += (seqx[r1] == seqy[r2]);
                            }
                        }
                        delete[] invmap;
                    }
                    else
                        TMalign_main(
                            xa, ya, seqx, seqy, secx, secy,
                            t0, u0, TM1, TM2, TM3, TM4, TM5,
                            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                            seqM, seqxA, seqyA, do_vec,
                            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                            xlen, ylen, sequence, Lnorm_ass, d0_scale,
                            i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                            mol_vec1[chain_i] + mol_vec2[chain_j], TMcut, 0);

                    /* print result */
                    if (outfmt_opt == 0)
                        print_version();
                    int left_num = 0;
                    int right_num = 0;
                    int left_aln_num = 0;
                    int right_aln_num = 0;
                    bool after_cp = false;
                    if (cp_opt)
                        after_cp = output_cp(
                            xname.substr(dir1_opt.size() + dir_opt.size()),
                            yname.substr(dir2_opt.size() + dir_opt.size()),
                            seqxA, seqyA, outfmt_opt, left_num, right_num,
                            left_aln_num, right_aln_num);
                    output_results(
                        xname.substr(dir1_opt.size() + dir_opt.size() + dirpair_opt.size()),
                        yname.substr(dir2_opt.size() + dir_opt.size() + dirpair_opt.size()),
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5,
                        rmsd0, d0_out, seqM.c_str(),
                        seqxA.c_str(), seqyA.c_str(), Liden,
                        n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
                        d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u,
                        (m_opt ? fname_matrix : "").c_str(),
                        outfmt_opt, ter_opt, false, split_opt, o_opt,
                        fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
                        resi_vec1, resi_vec2);
                    if (do_opt || (cp_opt && outfmt_opt <= 0))
                    {
                        cout << "###############\t###############\t#########" << endl;
                        cout << "#Aligned atom 1\tAligned atom 2 \tDistance#" << endl;
                        size_t r1 = right_num;
                        size_t r2 = 0;
                        size_t r;
                        int postcp = 0;
                        for (r = 0; r < seqxA.size(); r++)
                        {
                            r1 += seqxA[r] != '-';
                            r2 += seqyA[r] != '-';
                            if (seqxA[r] == '*')
                            {
                                cout << "###### Circular\tPermutation ###\t#########\n";
                                r1 = 0;
                                postcp = 1;
                            }
                            else if (seqxA[r] != '-' && seqyA[r] != '-')
                            {
                                cout << PDB_lines1[chain_i][r1 - 1].substr(12, 15) << '\t'
                                     << PDB_lines2[chain_j][r2 - 1].substr(12, 15) << '\t'
                                     << setw(9) << setiosflags(ios::fixed) << setprecision(3)
                                     << do_vec[r - postcp] << '\n';
                            }
                        }
                        cout << "###############\t###############\t#########" << endl;
                    }

                    /* Done! Free memory */
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();
                    DeleteArray(&ya, ylen);
                    delete[] seqy;
                    delete[] secy;
                    resi_vec2.clear();
                    do_vec.clear();
                } // chain_j
                if (chain2_list.size() > 1)
                {
                    yname.clear();
                    for (chain_j = 0; chain_j < ychainnum; chain_j++)
                        PDB_lines2[chain_j].clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
            PDB_lines1[chain_i].clear();
            DeleteArray(&xa, xlen);
            delete[] seqx;
            delete[] secx;
            resi_vec1.clear();
        } // chain_i
        xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
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

/* MMalign if more than two chains. TMalign if only one chain */
int MMalign(const string &xname, const string &yname,
            const string &fname_super, const string &fname_lign,
            const string &fname_matrix, vector<string> &sequence,
            const double d0_scale, const bool m_opt, const int o_opt,
            const int a_opt, const bool d_opt, const bool full_opt,
            const double TMcut, const int infmt1_opt, const int infmt2_opt,
            const int ter_opt, const int split_opt, const int outfmt_opt,
            bool fast_opt, const int mirror_opt, const int het_opt,
            const string &atom_opt, const bool autojustify, const string &mol_opt,
            const string &dir1_opt, const string &dir2_opt,
            const vector<string> &chain2parse1, const vector<string> &chain2parse2,
            const vector<string> &model2parse1, const vector<string> &model2parse2,
            const vector<string> &chain1_list, const vector<string> &chain2_list,
            const int byresi_opt, const string &chainmapfile, const bool se_opt)
{
    /* declare previously global variables */
    vector<vector<vector<double>>> xa_vec; // structure of complex1
    vector<vector<vector<double>>> ya_vec; // structure of complex2
    vector<vector<char>> seqx_vec;         // sequence of complex1
    vector<vector<char>> seqy_vec;         // sequence of complex2
    vector<vector<char>> secx_vec;         // secondary structure of complex1
    vector<vector<char>> secy_vec;         // secondary structure of complex2
    vector<int> mol_vec1;                  // molecule type of complex1, RNA if >0
    vector<int> mol_vec2;                  // molecule type of complex2, RNA if >0
    vector<string> chainID_list1;          // list of chainID1
    vector<string> chainID_list2;          // list of chainID2
    vector<int> xlen_vec;                  // length of complex1
    vector<int> ylen_vec;                  // length of complex2
    int i, j;                              // chain index
    int xlen, ylen;                        // chain length
    double **xa, **ya;                     // structure of single chain
    char *seqx, *seqy;                     // for the protein sequence
    char *secx, *secy;                     // for the secondary structure
    int xlen_aa, ylen_aa;                  // total length of protein
    int xlen_na, ylen_na;                  // total length of RNA/DNA
    vector<string> resi_vec1;              // residue index for chain1
    vector<string> resi_vec2;              // residue index for chain2

    /* parse complex */
    parse_chain_list(chain1_list, xa_vec, seqx_vec, secx_vec, mol_vec1,
                     xlen_vec, chainID_list1, ter_opt, split_opt, mol_opt, infmt1_opt,
                     atom_opt, autojustify, mirror_opt, het_opt, xlen_aa, xlen_na, o_opt,
                     resi_vec1, chain2parse1, model2parse1);
    if (xa_vec.size() == 0)
        PrintErrorAndQuit("ERROR! 0 chain in complex 1");
    parse_chain_list(chain2_list, ya_vec, seqy_vec, secy_vec, mol_vec2,
                     ylen_vec, chainID_list2, ter_opt, split_opt, mol_opt, infmt2_opt,
                     atom_opt, autojustify, 0, het_opt, ylen_aa, ylen_na, o_opt,
                     resi_vec2, chain2parse2, model2parse2);
    if (ya_vec.size() == 0)
        PrintErrorAndQuit("ERROR! 0 chain in complex 2");
    int len_aa = getmin(xlen_aa, ylen_aa);
    int len_na = getmin(xlen_na, ylen_na);
    if (a_opt)
    {
        len_aa = (xlen_aa + ylen_aa) / 2;
        len_na = (xlen_na + ylen_na) / 2;
    }
    int i_opt = 0;
    if (byresi_opt)
        i_opt = 3;

    map<int, int> chainmap;
    if (chainmapfile.size())
    {
        string line;
        int chainidx1, chainidx2;
        vector<string> line_vec;
        ifstream fin;
        bool fromStdin = (chainmapfile == "-");
        if (!fromStdin)
            fin.open(chainmapfile.c_str());
        while (fromStdin ? cin.good() : fin.good())
        {
            if (fromStdin)
                getline(cin, line);
            else
                getline(fin, line);
            if (line.size() == 0 || line[0] == '#')
                continue;
            split(line, line_vec, '\t');
            if (line_vec.size() == 2)
            {
                chainidx1 = -1;
                chainidx2 = -1;

                for (i = 0; i < chainID_list1.size(); i++)
                {
                    if (line_vec[0] == chainID_list1[i] ||
                        ":" + line_vec[0] == chainID_list1[i] ||
                        ":1," + line_vec[0] == chainID_list1[i])
                    {
                        chainidx1 = i;
                        break;
                    }
                }
                for (i = 0; i < chainID_list2.size(); i++)
                {
                    if (line_vec[1] == chainID_list2[i] ||
                        ":" + line_vec[1] == chainID_list2[i] ||
                        ":1," + line_vec[1] == chainID_list2[i])
                    {
                        chainidx2 = i;
                        break;
                    }
                }
                if (chainidx1 >= 0 && chainidx2 >= 0)
                {
                    if (chainmap.count(chainidx1))
                        cerr << "ERROR! " << line_vec[0] << " already mapped" << endl;
                    chainmap[chainidx1] = chainidx2;
                }
                else
                    cerr << "ERROR! Cannot map " << line << endl;
            }
            else
                cerr << "ERROR! Cannot map " << line << endl;
            for (i = 0; i < line_vec.size(); i++)
                line_vec[i].clear();
            line_vec.clear();
        }
        if (!fromStdin)
            fin.close();
        if (chainmap.size() == 0)
            cerr << "ERROR! cannot map any chain pair from " << chainmapfile << endl;
    }

    /* perform monomer alignment if there is only one chain */
    if (xa_vec.size() == 1 && ya_vec.size() == 1)
    {
        xlen = xlen_vec[0];
        ylen = ylen_vec[0];
        seqx = new char[xlen + 1];
        seqy = new char[ylen + 1];
        secx = new char[xlen + 1];
        secy = new char[ylen + 1];
        NewArray(&xa, xlen, 3);
        NewArray(&ya, ylen, 3);
        copy_chain_data(xa_vec[0], seqx_vec[0], secx_vec[0], xlen, xa, seqx, secx);
        copy_chain_data(ya_vec[0], seqy_vec[0], secy_vec[0], ylen, ya, seqy, secy);

        /* declare variable specific to this pair of TMalign */
        double t0[3], u0[3][3];
        double TM1, TM2;
        double TM3, TM4, TM5; // for a_opt, u_opt, d_opt
        double d0_0, TM_0;
        double d0A, d0B, d0u, d0a;
        double d0_out = 5.0;
        string seqM, seqxA, seqyA; // for output alignment
        double rmsd0 = 0.0;
        int L_ali; // Aligned length in standard_TMscore
        double Liden = 0;
        double TM_ali, rmsd_ali; // TMscore and rmsd in standard_TMscore
        int n_ali = 0;
        int n_ali8 = 0;
        vector<double> do_vec;

        if (byresi_opt)
            extract_aln_from_resi(sequence,
                                  seqx, seqy, resi_vec1, resi_vec2, byresi_opt);

        /* entry function for structure alignment */
        if (se_opt)
        {
            int *invmap = new int[ylen + 1];
            u0[0][0] = u0[1][1] = u0[2][2] = 1;
            u0[0][1] = u0[0][2] =
                u0[1][0] = u0[1][2] =
                    u0[2][0] = u0[2][1] =
                        t0[0] = t0[1] = t0[2] = 0;
            se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                    seqM, seqxA, seqyA, do_vec,
                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, 0, d0_scale,
                    i_opt, a_opt, false, d_opt,
                    mol_vec1[0] + mol_vec2[0], outfmt_opt, invmap);
            if (outfmt_opt >= 2)
            {
                Liden = L_ali = 0;
                int r1, r2;
                for (r2 = 0; r2 < ylen; r2++)
                {
                    r1 = invmap[r2];
                    if (r1 < 0)
                        continue;
                    L_ali += 1;
                    Liden += (seqx[r1] == seqy[r2]);
                }
            }
            delete[] invmap;
        }
        else
            TMalign_main(xa, ya, seqx, seqy, secx, secy,
                         t0, u0, TM1, TM2, TM3, TM4, TM5,
                         d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                         seqM, seqxA, seqyA, do_vec,
                         rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                         xlen, ylen, sequence, 0, d0_scale,
                         i_opt, a_opt, false, d_opt, fast_opt,
                         mol_vec1[0] + mol_vec2[0], TMcut, 0);

        /* print result */
        output_results(
            xname.substr(dir1_opt.size()),
            yname.substr(dir2_opt.size()),
            chainID_list1[0], chainID_list2[0],
            xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
            seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden,
            n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, d0B,
            0, d0_scale, d0a, d0u, (m_opt ? fname_matrix : "").c_str(),
            outfmt_opt, ter_opt, true, split_opt, o_opt, fname_super,
            0, a_opt, false, d_opt, mirror_opt, resi_vec1, resi_vec2);

        /* clean up */
        seqM.clear();
        seqxA.clear();
        seqyA.clear();
        delete[] seqx;
        delete[] seqy;
        delete[] secx;
        delete[] secy;
        DeleteArray(&xa, xlen);
        DeleteArray(&ya, ylen);
        do_vec.clear();

        vector<vector<vector<double>>>().swap(xa_vec); // structure of complex1
        vector<vector<vector<double>>>().swap(ya_vec); // structure of complex2
        vector<vector<char>>().swap(seqx_vec);         // sequence of complex1
        vector<vector<char>>().swap(seqy_vec);         // sequence of complex2
        vector<vector<char>>().swap(secx_vec);         // secondary structure of complex1
        vector<vector<char>>().swap(secy_vec);         // secondary structure of complex2
        mol_vec1.clear();                              // molecule type of complex1, RNA if >0
        mol_vec2.clear();                              // molecule type of complex2, RNA if >0
        chainID_list1.clear();                         // list of chainID1
        chainID_list2.clear();                         // list of chainID2
        xlen_vec.clear();                              // length of complex1
        ylen_vec.clear();                              // length of complex2
        return 0;
    }

    /* declare TM-score tables */
    int chain1_num = xa_vec.size();
    int chain2_num = ya_vec.size();
    int chain_num = MAX(chain1_num, chain2_num);
    vector<string> tmp_str_vec(chain2_num, "");
    double **TMave_mat;
    double **ut_mat; // rotation matrices for all-against-all alignment
    int ui, uj, ut_idx;
    NewArray(&TMave_mat, chain_num, chain_num);
    NewArray(&ut_mat, chain1_num * chain2_num, 4 * 3);
    vector<vector<string>> seqxA_mat(chain1_num, tmp_str_vec);
    vector<vector<string>> seqM_mat(chain1_num, tmp_str_vec);
    vector<vector<string>> seqyA_mat(chain1_num, tmp_str_vec);

    double maxTMmono = -1;
    int maxTMmono_i, maxTMmono_j;

    /* get all-against-all alignment */
    if (len_aa + len_na > 500)
        fast_opt = true;
    for (i = 0; i < chain1_num; i++)
    {
        xlen = xlen_vec[i];
        if (xlen < 3)
        {
            for (j = 0; j < chain2_num; j++)
                TMave_mat[i][j] = TMave_mat[j][i] = -1;
            continue;
        }
        seqx = new char[xlen + 1];
        secx = new char[xlen + 1];
        NewArray(&xa, xlen, 3);
        copy_chain_data(xa_vec[i], seqx_vec[i], secx_vec[i],
                        xlen, xa, seqx, secx);

        for (j = 0; j < chain2_num; j++)
        {
            ut_idx = i * chain2_num + j;
            for (ui = 0; ui < 4; ui++)
                for (uj = 0; uj < 3; uj++)
                    ut_mat[ut_idx][ui * 3 + uj] = 0;
            ut_mat[ut_idx][0] = 1;
            ut_mat[ut_idx][4] = 1;
            ut_mat[ut_idx][8] = 1;

            if (mol_vec1[i] * mol_vec2[j] < 0) // no protein-RNA alignment
            {
                TMave_mat[i][j] = TMave_mat[j][i] = -1;
                continue;
            }
            if (chainmap.size() && (!chainmap.count(i) || chainmap[i] != j))
            {
                TMave_mat[i][j] = TMave_mat[j][i] = -1;
                continue;
            }

            ylen = ylen_vec[j];
            if (ylen < 3)
            {
                TMave_mat[i][j] = TMave_mat[j][i] = -1;
                continue;
            }
            seqy = new char[ylen + 1];
            secy = new char[ylen + 1];
            NewArray(&ya, ylen, 3);
            copy_chain_data(ya_vec[j], seqy_vec[j], secy_vec[j],
                            ylen, ya, seqy, secy);

            /* declare variable specific to this pair of TMalign */
            double t0[3], u0[3][3];
            double TM1, TM2;
            double TM3, TM4, TM5; // for a_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out = 5.0;
            string seqM, seqxA, seqyA; // for output alignment
            double rmsd0 = 0.0;
            int L_ali; // Aligned length in standard_TMscore
            double Liden = 0;
            double TM_ali, rmsd_ali; // TMscore and rmsd in standard_TMscore
            int n_ali = 0;
            int n_ali8 = 0;
            vector<double> do_vec;

            int Lnorm_tmp = len_aa;
            if (mol_vec1[i] + mol_vec2[j] > 0)
                Lnorm_tmp = len_na;

            if (byresi_opt)
            {
                int total_aln = extract_aln_from_resi(sequence, seqx, seqy,
                                                      resi_vec1, resi_vec2, xlen_vec, ylen_vec, i, j, byresi_opt);
                seqxA_mat[i][j] = sequence[0];
                seqyA_mat[i][j] = sequence[1];
                if (total_aln > xlen + ylen - 3)
                {
                    for (ui = 0; ui < 3; ui++)
                        for (uj = 0; uj < 3; uj++)
                            ut_mat[ut_idx][ui * 3 + uj] = (ui == uj) ? 1 : 0;
                    for (uj = 0; uj < 3; uj++)
                        ut_mat[ut_idx][9 + uj] = 0;
                    TMave_mat[i][j] = TMave_mat[j][i] = 0;
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();

                    delete[] seqy;
                    delete[] secy;
                    DeleteArray(&ya, ylen);
                    continue;
                }
            }

            /* entry function for structure alignment */
            if (se_opt)
            {
                int *invmap = new int[ylen + 1];
                u0[0][0] = u0[1][1] = u0[2][2] = 1;
                u0[0][1] = u0[0][2] =
                    u0[1][0] = u0[1][2] =
                        u0[2][0] = u0[2][1] =
                            t0[0] = t0[1] = t0[2] = 0;
                se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA, do_vec,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                        i_opt, false, true, false,
                        mol_vec1[i] + mol_vec2[j], outfmt_opt, invmap);
                if (outfmt_opt >= 2)
                {
                    Liden = L_ali = 0;
                    int r1, r2;
                    for (r2 = 0; r2 < ylen; r2++)
                    {
                        r1 = invmap[r2];
                        if (r1 < 0)
                            continue;
                        L_ali += 1;
                        Liden += (seqx[r1] == seqy[r2]);
                    }
                }
                delete[] invmap;
            }
            else
                TMalign_main(xa, ya, seqx, seqy, secx, secy,
                             t0, u0, TM1, TM2, TM3, TM4, TM5,
                             d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                             seqM, seqxA, seqyA, do_vec,
                             rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                             xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                             i_opt, false, true, false, fast_opt,
                             mol_vec1[i] + mol_vec2[j], TMcut, 0);

            /* store result */
            for (ui = 0; ui < 3; ui++)
                for (uj = 0; uj < 3; uj++)
                    ut_mat[ut_idx][ui * 3 + uj] = u0[ui][uj];
            for (uj = 0; uj < 3; uj++)
                ut_mat[ut_idx][9 + uj] = t0[uj];
            seqxA_mat[i][j] = seqxA;
            seqyA_mat[i][j] = seqyA;
            TMave_mat[i][j] = TMave_mat[j][i] = TM4 * Lnorm_tmp;
            if (TMave_mat[i][j] > maxTMmono)
            {
                maxTMmono = TMave_mat[i][j];
                maxTMmono_i = i;
                maxTMmono_j = j;
            }

            /* clean up */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();

            delete[] seqy;
            delete[] secy;
            DeleteArray(&ya, ylen);
            do_vec.clear();
        }

        delete[] seqx;
        delete[] secx;
        DeleteArray(&xa, xlen);
    }

    /* calculate initial chain-chain assignment */
    int *assign1_list; // value is index of assigned chain2
    int *assign2_list; // value is index of assigned chain1
    assign1_list = new int[chain1_num];
    assign2_list = new int[chain2_num];
    double total_score = enhanced_greedy_search(TMave_mat, assign1_list,
                                                assign2_list, chain1_num, chain2_num);
    if (total_score <= 0)
        PrintErrorAndQuit("ERROR! No assignable chain");

    /* refine alignment for large oligomers */
    int aln_chain_num = count_assign_pair(assign1_list, chain1_num);
    bool is_oligomer = (aln_chain_num >= 3);
    if (aln_chain_num == 2 && chainmap.size() == 0 && !se_opt) // dimer alignment
    {
        int na_chain_num1, na_chain_num2, aa_chain_num1, aa_chain_num2;
        count_na_aa_chain_num(na_chain_num1, aa_chain_num1, mol_vec1);
        count_na_aa_chain_num(na_chain_num2, aa_chain_num2, mol_vec2);

        /* align protein-RNA hybrid dimer to another hybrid dimer */
        if (na_chain_num1 == 1 && na_chain_num2 == 1 &&
            aa_chain_num1 == 1 && aa_chain_num2 == 1)
            is_oligomer = false;
        /* align pure protein dimer or pure RNA dimer */
        else if ((getmin(na_chain_num1, na_chain_num2) == 0 &&
                  aa_chain_num1 == 2 && aa_chain_num2 == 2) ||
                 (getmin(aa_chain_num1, aa_chain_num2) == 0 &&
                  na_chain_num1 == 2 && na_chain_num2 == 2))
        {
            adjust_dimer_assignment(xa_vec, ya_vec, xlen_vec, ylen_vec, mol_vec1,
                                    mol_vec2, assign1_list, assign2_list, seqxA_mat, seqyA_mat);
            is_oligomer = false; // cannot refiner further
        }
        else
            is_oligomer = true; /* align oligomers to dimer */
    }

    if ((aln_chain_num >= 3 || is_oligomer) && chainmap.size() == 0 && !se_opt) // oligomer alignment
    {
        /* extract centroid coordinates */
        double **xcentroids;
        double **ycentroids;
        NewArray(&xcentroids, chain1_num, 3);
        NewArray(&ycentroids, chain2_num, 3);
        double d0MM = getmin(
            calculate_centroids(xa_vec, chain1_num, xcentroids),
            calculate_centroids(ya_vec, chain2_num, ycentroids));

        /* refine enhanced greedy search with centroid superposition */
        // double het_deg=check_heterooligomer(TMave_mat, chain1_num, chain2_num);
        homo_refined_greedy_search(TMave_mat, assign1_list,
                                   assign2_list, chain1_num, chain2_num, xcentroids,
                                   ycentroids, d0MM, len_aa + len_na, ut_mat);

        if (chain1_num <= chain2_num)
        {
            hetero_refined_greedy_search(TMave_mat, assign1_list,
                                         assign2_list, chain1_num, chain2_num, xcentroids,
                                         ycentroids, d0MM, len_aa + len_na);
        }
        else
        {
            hetero_refined_greedy_search(TMave_mat, assign2_list,
                                         assign1_list, chain2_num, chain1_num, ycentroids,
                                         xcentroids, d0MM, len_aa + len_na);
        }

        /* clean up */
        DeleteArray(&xcentroids, chain1_num);
        DeleteArray(&ycentroids, chain2_num);
    }

    /* store initial assignment */
    int init_pair_num = count_assign_pair(assign1_list, chain1_num);
    int *assign1_init, *assign2_init;
    assign1_init = new int[chain1_num];
    assign2_init = new int[chain2_num];
    double **TMave_init;
    NewArray(&TMave_init, chain1_num, chain2_num);
    vector<vector<string>> seqxA_init(chain1_num, tmp_str_vec);
    vector<vector<string>> seqyA_init(chain1_num, tmp_str_vec);
    vector<string> sequence_init;
    copy_chain_assign_data(chain1_num, chain2_num, sequence_init,
                           seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat,
                           seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init);

    /* perform iterative alignment */
    double max_total_score = 0; // ignore old total_score because previous
                                // score was from monomeric chain superpositions
    int max_iter = 5 - (int)((len_aa + len_na) / 200);
    if (max_iter < 2)
        max_iter = 2;
    // if (byresi_opt==0)
    if (!se_opt)
        MMalign_iter(max_total_score, max_iter, xa_vec, ya_vec,
                     seqx_vec, seqy_vec, secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec,
                     ylen_vec, xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num,
                     chain2_num, TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list,
                     sequence, d0_scale, fast_opt, chainmap, byresi_opt);

    if (byresi_opt && aln_chain_num >= 4 && is_oligomer && chainmap.size() == 0 && !se_opt) // oligomer alignment
    {
        MMalign_final(xname.substr(dir1_opt.size()), yname.substr(dir2_opt.size()),
                      chainID_list1, chainID_list2,
                      fname_super, fname_lign, fname_matrix,
                      xa_vec, ya_vec, seqx_vec, seqy_vec,
                      secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
                      xa, ya, seqx, seqy, secx, secy, len_aa, len_na,
                      chain1_num, chain2_num, TMave_mat,
                      seqxA_mat, seqM_mat, seqyA_mat, assign1_list, assign2_list, sequence,
                      d0_scale, 1, 0, 5, ter_opt, split_opt,
                      0, 0, true, true, mirror_opt, resi_vec1, resi_vec2);

        /* extract centroid coordinates */
        double **xcentroids;
        double **ycentroids;
        NewArray(&xcentroids, chain1_num, 3);
        NewArray(&ycentroids, chain2_num, 3);
        double d0MM = getmin(
            calculate_centroids(xa_vec, chain1_num, xcentroids),
            calculate_centroids(ya_vec, chain2_num, ycentroids));

        /* refine enhanced greedy search with centroid superposition */
        // double het_deg=check_heterooligomer(TMave_mat, chain1_num, chain2_num);
        homo_refined_greedy_search(TMave_mat, assign1_list,
                                   assign2_list, chain1_num, chain2_num, xcentroids,
                                   ycentroids, d0MM, len_aa + len_na, ut_mat);

        hetero_refined_greedy_search(TMave_mat, assign1_list,
                                     assign2_list, chain1_num, chain2_num, xcentroids,
                                     ycentroids, d0MM, len_aa + len_na);

        /* clean up */
        DeleteArray(&xcentroids, chain1_num);
        DeleteArray(&ycentroids, chain2_num);
    }

    /* sometime MMalign_iter is even worse than monomer alignment */
    if (byresi_opt == 0 && max_total_score < maxTMmono)
    {
        copy_chain_assign_data(chain1_num, chain2_num, sequence,
                               seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
                               seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
        for (i = 0; i < chain1_num; i++)
        {
            if (i != maxTMmono_i)
                assign1_list[i] = -1;
            else
                assign1_list[i] = maxTMmono_j;
        }
        for (j = 0; j < chain2_num; j++)
        {
            if (j != maxTMmono_j)
                assign2_list[j] = -1;
            else
                assign2_list[j] = maxTMmono_i;
        }
        sequence[0] = seqxA_mat[maxTMmono_i][maxTMmono_j];
        sequence[1] = seqyA_mat[maxTMmono_i][maxTMmono_j];
        max_total_score = maxTMmono;
        MMalign_iter(max_total_score, max_iter, xa_vec, ya_vec, seqx_vec, seqy_vec,
                     secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
                     xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num, chain2_num,
                     TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list, sequence,
                     d0_scale, fast_opt, chainmap);
    }

    /* perform cross chain alignment
     * in some cases, this leads to dramatic improvement, esp for homodimer */
    int iter_pair_num = count_assign_pair(assign1_list, chain1_num);
    if (iter_pair_num >= init_pair_num)
        copy_chain_assign_data(
            chain1_num, chain2_num, sequence_init,
            seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat,
            seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init);
    double max_total_score_cross = max_total_score;
    if (byresi_opt == 0 && len_aa + len_na < 10000)
    {
        MMalign_dimer(max_total_score_cross, xa_vec, ya_vec, seqx_vec, seqy_vec,
                      secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
                      xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num, chain2_num,
                      TMave_init, seqxA_init, seqyA_init, assign1_init, assign2_init,
                      sequence_init, d0_scale, fast_opt);
        if (max_total_score_cross > max_total_score)
        {
            max_total_score = max_total_score_cross;
            copy_chain_assign_data(chain1_num, chain2_num, sequence,
                                   seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
                                   seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
        }
    }

    /* final alignment */
    if (outfmt_opt == 0)
        print_version();
    if (se_opt)
        MMalign_se_final(xname.substr(dir1_opt.size()), yname.substr(dir2_opt.size()),
                         chainID_list1, chainID_list2,
                         fname_super, fname_lign, fname_matrix,
                         xa_vec, ya_vec, seqx_vec, seqy_vec,
                         secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
                         xa, ya, seqx, seqy, secx, secy, len_aa, len_na,
                         chain1_num, chain2_num, TMave_mat,
                         seqxA_mat, seqM_mat, seqyA_mat, assign1_list, assign2_list, sequence,
                         d0_scale, m_opt, o_opt, outfmt_opt, ter_opt, split_opt,
                         a_opt, d_opt, fast_opt, full_opt, mirror_opt, resi_vec1, resi_vec2);
    else
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
    delete[] assign1_list;
    delete[] assign2_list;
    DeleteArray(&TMave_mat, chain_num);
    DeleteArray(&ut_mat, chain1_num * chain2_num);
    vector<vector<string>>().swap(seqxA_mat);
    vector<vector<string>>().swap(seqM_mat);
    vector<vector<string>>().swap(seqyA_mat);
    vector<string>().swap(tmp_str_vec);

    delete[] assign1_init;
    delete[] assign2_init;
    DeleteArray(&TMave_init, chain1_num);
    vector<vector<string>>().swap(seqxA_init);
    vector<vector<string>>().swap(seqyA_init);

    vector<vector<vector<double>>>().swap(xa_vec); // structure of complex1
    vector<vector<vector<double>>>().swap(ya_vec); // structure of complex2
    vector<vector<char>>().swap(seqx_vec);         // sequence of complex1
    vector<vector<char>>().swap(seqy_vec);         // sequence of complex2
    vector<vector<char>>().swap(secx_vec);         // secondary structure of complex1
    vector<vector<char>>().swap(secy_vec);         // secondary structure of complex2
    mol_vec1.clear();                              // molecule type of complex1, RNA if >0
    mol_vec2.clear();                              // molecule type of complex2, RNA if >0
    vector<string>().swap(chainID_list1);          // list of chainID1
    vector<string>().swap(chainID_list2);          // list of chainID2
    xlen_vec.clear();                              // length of complex1
    ylen_vec.clear();                              // length of complex2
    vector<string>().swap(resi_vec1);              // residue index for chain1
    vector<string>().swap(resi_vec2);              // residue index for chain2
    map<int, int>().swap(chainmap);
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
           const string &atom_opt, const bool autojustify, const string &mol_opt,
           const string &dir1_opt, const string &dir2_opt,
           const vector<string> &chain2parse1, const vector<string> &chain2parse2,
           const vector<string> &model2parse1, const vector<string> &model2parse2,
           const vector<string> &chain1_list, const vector<string> &chain2_list,
           const bool do_opt)
{
    /* declare previously global variables */
    vector<vector<vector<double>>> xa_vec; // structure of complex1
    vector<vector<vector<double>>> ya_vec; // structure of complex2
    vector<vector<char>> seqx_vec;         // sequence of complex1
    vector<vector<char>> seqy_vec;         // sequence of complex2
    vector<vector<char>> secx_vec;         // secondary structure of complex1
    vector<vector<char>> secy_vec;         // secondary structure of complex2
    vector<int> mol_vec1;                  // molecule type of complex1, RNA if >0
    vector<int> mol_vec2;                  // molecule type of complex2, RNA if >0
    vector<string> chainID_list1;          // list of chainID1
    vector<string> chainID_list2;          // list of chainID2
    vector<int> xlen_vec;                  // length of complex1
    vector<int> ylen_vec;                  // length of complex2
    int i, j;                              // chain index
    int xlen, ylen;                        // chain length
    double **xa, **ya;                     // structure of single chain
    char *seqx, *seqy;                     // for the protein sequence
    char *secx, *secy;                     // for the secondary structure
    int xlen_aa, ylen_aa;                  // total length of protein
    int xlen_na, ylen_na;                  // total length of RNA/DNA
    vector<string> resi_vec1;              // residue index for chain1
    vector<string> resi_vec2;              // residue index for chain2

    /* parse complex */
    parse_chain_list(chain1_list, xa_vec, seqx_vec, secx_vec, mol_vec1,
                     xlen_vec, chainID_list1, ter_opt, split_opt, mol_opt, infmt1_opt,
                     atom_opt, autojustify, mirror_opt, het_opt, xlen_aa, xlen_na, o_opt,
                     resi_vec1, chain2parse1, model2parse1);
    if (xa_vec.size() == 0)
        PrintErrorAndQuit("ERROR! 0 individual chain");
    parse_chain_list(chain2_list, ya_vec, seqy_vec, secy_vec, mol_vec2,
                     ylen_vec, chainID_list2, ter_opt, split_opt, mol_opt, infmt2_opt,
                     atom_opt, autojustify, 0, het_opt, ylen_aa, ylen_na, o_opt, resi_vec2,
                     chain2parse2, model2parse2);
    if (xa_vec.size() > ya_vec.size())
        PrintErrorAndQuit(
            "ERROR! more individual chains to align than number of chains in complex template");
    int len_aa = getmin(xlen_aa, ylen_aa);
    int len_na = getmin(xlen_na, ylen_na);
    if (a_opt)
    {
        len_aa = (xlen_aa + ylen_aa) / 2;
        len_na = (xlen_na + ylen_na) / 2;
    }

    /* perform monomer alignment if there is only one chain */
    if (xa_vec.size() == 1 && ya_vec.size() == 1)
    {
        xlen = xlen_vec[0];
        ylen = ylen_vec[0];
        seqx = new char[xlen + 1];
        seqy = new char[ylen + 1];
        secx = new char[xlen + 1];
        secy = new char[ylen + 1];
        NewArray(&xa, xlen, 3);
        NewArray(&ya, ylen, 3);
        copy_chain_data(xa_vec[0], seqx_vec[0], secx_vec[0], xlen, xa, seqx, secx);
        copy_chain_data(ya_vec[0], seqy_vec[0], secy_vec[0], ylen, ya, seqy, secy);

        /* declare variable specific to this pair of TMalign */
        double t0[3], u0[3][3];
        double TM1, TM2;
        double TM3, TM4, TM5; // for a_opt, u_opt, d_opt
        double d0_0, TM_0;
        double d0A, d0B, d0u, d0a;
        double d0_out = 5.0;
        string seqM, seqxA, seqyA; // for output alignment
        double rmsd0 = 0.0;
        int L_ali; // Aligned length in standard_TMscore
        double Liden = 0;
        double TM_ali, rmsd_ali; // TMscore and rmsd in standard_TMscore
        int n_ali = 0;
        int n_ali8 = 0;
        vector<double> do_vec;

        /* entry function for structure alignment */
        TMalign_main(xa, ya, seqx, seqy, secx, secy,
                     t0, u0, TM1, TM2, TM3, TM4, TM5,
                     d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                     seqM, seqxA, seqyA, do_vec,
                     rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                     xlen, ylen, sequence, Lnorm_ass, d0_scale,
                     0, a_opt, u_opt, d_opt, fast_opt,
                     mol_vec1[0] + mol_vec2[0], TMcut, 0);

        /* print result */
        output_results(
            xname.substr(dir1_opt.size()),
            yname.substr(dir2_opt.size()),
            chainID_list1[0], chainID_list2[0],
            xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
            seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden,
            n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, d0B,
            Lnorm_ass, d0_scale, d0a, d0u, (m_opt ? fname_matrix : "").c_str(),
            (outfmt_opt == 2 ? outfmt_opt : 3), ter_opt, true, split_opt, o_opt, fname_super,
            0, a_opt, false, d_opt, mirror_opt, resi_vec1, resi_vec2);
        if (outfmt_opt == 2)
            printf("%s%s\t%s%s\t%.4f\n",
                   xname.substr(dir1_opt.size()).c_str(), chainID_list1[0].c_str(),
                   yname.substr(dir2_opt.size()).c_str(), chainID_list2[0].c_str(),
                   sqrt((TM1 * TM1 + TM2 * TM2) / 2));

        /* clean up */
        seqM.clear();
        seqxA.clear();
        seqyA.clear();
        delete[] seqx;
        delete[] seqy;
        delete[] secx;
        delete[] secy;
        DeleteArray(&xa, xlen);
        DeleteArray(&ya, ylen);
        do_vec.clear();

        vector<vector<vector<double>>>().swap(xa_vec); // structure of complex1
        vector<vector<vector<double>>>().swap(ya_vec); // structure of complex2
        vector<vector<char>>().swap(seqx_vec);         // sequence of complex1
        vector<vector<char>>().swap(seqy_vec);         // sequence of complex2
        vector<vector<char>>().swap(secx_vec);         // secondary structure of complex1
        vector<vector<char>>().swap(secy_vec);         // secondary structure of complex2
        mol_vec1.clear();                              // molecule type of complex1, RNA if >0
        mol_vec2.clear();                              // molecule type of complex2, RNA if >0
        chainID_list1.clear();                         // list of chainID1
        chainID_list2.clear();                         // list of chainID2
        xlen_vec.clear();                              // length of complex1
        ylen_vec.clear();                              // length of complex2
        return 0;
    }

    /* declare TM-score tables */
    int chain1_num = xa_vec.size();
    int chain2_num = ya_vec.size();
    vector<string> tmp_str_vec(chain2_num, "");
    double **TMave_mat;
    NewArray(&TMave_mat, chain1_num, chain2_num);
    vector<vector<string>> seqxA_mat(chain1_num, tmp_str_vec);
    vector<vector<string>> seqM_mat(chain1_num, tmp_str_vec);
    vector<vector<string>> seqyA_mat(chain1_num, tmp_str_vec);

    /* trimComplex */
    vector<vector<vector<double>>> ya_trim_vec; // structure of complex2
    vector<vector<char>> seqy_trim_vec;         // sequence of complex2
    vector<vector<char>> secy_trim_vec;         // secondary structure of complex2
    vector<int> ylen_trim_vec;                  // length of complex2
    int Lchain_aa_max1 = 0;
    int Lchain_na_max1 = 0;
    for (i = 0; i < chain1_num; i++)
    {
        xlen = xlen_vec[i];
        if (mol_vec1[i] > 0 && xlen > Lchain_na_max1)
            Lchain_na_max1 = xlen;
        else if (mol_vec1[i] <= 0 && xlen > Lchain_aa_max1)
            Lchain_aa_max1 = xlen;
    }
    int trim_chain_count = trimComplex(ya_trim_vec, seqy_trim_vec,
                                       secy_trim_vec, ylen_trim_vec, ya_vec, seqy_vec, secy_vec, ylen_vec,
                                       mol_vec2, Lchain_aa_max1, Lchain_na_max1);
    int ylen_trim;    // chain length
    double **ya_trim; // structure of single chain
    char *seqy_trim;  // for the protein sequence
    char *secy_trim;  // for the secondary structure
    double **xt;

    /* get all-against-all alignment */
    if (len_aa + len_na > 500)
        fast_opt = true;
    for (i = 0; i < chain1_num; i++)
    {
        xlen = xlen_vec[i];
        if (xlen < 3)
        {
            for (j = 0; j < chain2_num; j++)
                TMave_mat[i][j] = -1;
            continue;
        }
        seqx = new char[xlen + 1];
        secx = new char[xlen + 1];
        NewArray(&xa, xlen, 3);
        copy_chain_data(xa_vec[i], seqx_vec[i], secx_vec[i],
                        xlen, xa, seqx, secx);

        for (j = 0; j < chain2_num; j++)
        {
            if (mol_vec1[i] * mol_vec2[j] < 0) // no protein-RNA alignment
            {
                TMave_mat[i][j] = -1;
                continue;
            }

            ylen = ylen_vec[j];
            if (ylen < 3)
            {
                TMave_mat[i][j] = -1;
                continue;
            }
            seqy = new char[ylen + 1];
            secy = new char[ylen + 1];
            NewArray(&ya, ylen, 3);
            copy_chain_data(ya_vec[j], seqy_vec[j], secy_vec[j],
                            ylen, ya, seqy, secy);

            /* declare variable specific to this pair of TMalign */
            double t0[3], u0[3][3];
            double TM1, TM2;
            double TM3, TM4, TM5; // for a_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out = 5.0;
            string seqM, seqxA, seqyA; // for output alignment
            double rmsd0 = 0.0;
            int L_ali; // Aligned length in standard_TMscore
            double Liden = 0;
            double TM_ali, rmsd_ali; // TMscore and rmsd in standard_TMscore
            int n_ali = 0;
            int n_ali8 = 0;
            vector<double> do_vec;

            int Lnorm_tmp = len_aa;
            if (mol_vec1[i] + mol_vec2[j] > 0)
                Lnorm_tmp = len_na;

            /* entry function for structure alignment */
            if (trim_chain_count && ylen_trim_vec[j] < ylen)
            {
                ylen_trim = ylen_trim_vec[j];
                seqy_trim = new char[ylen_trim + 1];
                secy_trim = new char[ylen_trim + 1];
                NewArray(&ya_trim, ylen_trim, 3);
                copy_chain_data(ya_trim_vec[j], seqy_trim_vec[j], secy_trim_vec[j],
                                ylen_trim, ya_trim, seqy_trim, secy_trim);
                TMalign_main(xa, ya_trim, seqx, seqy_trim, secx, secy_trim,
                             t0, u0, TM1, TM2, TM3, TM4, TM5,
                             d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                             seqM, seqxA, seqyA, do_vec,
                             rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                             xlen, ylen_trim, sequence, Lnorm_tmp, d0_scale,
                             0, false, true, false, fast_opt,
                             mol_vec1[i] + mol_vec2[j], TMcut, 0);
                seqxA.clear();
                seqyA.clear();
                delete[] seqy_trim;
                delete[] secy_trim;
                DeleteArray(&ya_trim, ylen_trim);

                NewArray(&xt, xlen, 3);
                do_rotation(xa, xt, xlen, t0, u0);
                int *invmap = new int[ylen + 1];
                se_main(xt, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
                        do_vec, rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                        0, false, 2, false, mol_vec1[i] + mol_vec2[j], 1, invmap);
                delete[] invmap;

                if (sequence.size() < 2)
                    sequence.push_back("");
                if (sequence.size() < 2)
                    sequence.push_back("");
                sequence[0] = seqxA;
                sequence[1] = seqyA;
                TMalign_main(xt, ya, seqx, seqy, secx, secy,
                             t0, u0, TM1, TM2, TM3, TM4, TM5,
                             d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                             seqM, seqxA, seqyA, do_vec,
                             rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                             xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                             2, false, true, false, fast_opt,
                             mol_vec1[i] + mol_vec2[j], TMcut, 0);
                DeleteArray(&xt, xlen);
            }
            else
            {
                TMalign_main(xa, ya, seqx, seqy, secx, secy,
                             t0, u0, TM1, TM2, TM3, TM4, TM5,
                             d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                             seqM, seqxA, seqyA, do_vec,
                             rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                             xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                             0, false, true, false, fast_opt,
                             mol_vec1[i] + mol_vec2[j], TMcut, 0);
            }

            /* store result */
            seqxA_mat[i][j] = seqxA;
            seqyA_mat[i][j] = seqyA;
            TMave_mat[i][j] = TM4 * Lnorm_tmp;

            /* clean up */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();
            do_vec.clear();

            delete[] seqy;
            delete[] secy;
            DeleteArray(&ya, ylen);
        }

        delete[] seqx;
        delete[] secx;
        DeleteArray(&xa, xlen);
    }
    vector<vector<vector<double>>>().swap(ya_trim_vec);
    vector<vector<char>>().swap(seqy_trim_vec);
    vector<vector<char>>().swap(secy_trim_vec);
    vector<int>().swap(ylen_trim_vec);

    /* calculate initial chain-chain assignment */
    int *assign1_list; // value is index of assigned chain2
    int *assign2_list; // value is index of assigned chain1
    assign1_list = new int[chain1_num];
    assign2_list = new int[chain2_num];
    enhanced_greedy_search(TMave_mat, assign1_list,
                           assign2_list, chain1_num, chain2_num);

    /* final alignment */
    if (outfmt_opt == 0)
        print_version();
    double **ut_mat; // rotation matrices for all-against-all alignment
    NewArray(&ut_mat, chain1_num, 4 * 3);
    int ui, uj;
    vector<string> xname_vec;
    vector<string> yname_vec;
    vector<double> TM_vec;
    for (i = 0; i < chain1_num; i++)
    {
        j = assign1_list[i];
        xname_vec.push_back(xname + chainID_list1[i]);
        if (j < 0)
        {
            cerr << "Warning! " << chainID_list1[i] << " cannot be alighed" << endl;
            for (ui = 0; ui < 3; ui++)
            {
                for (uj = 0; uj < 4; uj++)
                    ut_mat[i][ui * 3 + uj] = 0;
                ut_mat[i][ui * 3 + ui] = 1;
            }
            yname_vec.push_back(yname);
            continue;
        }
        yname_vec.push_back(yname + chainID_list2[j]);

        xlen = xlen_vec[i];
        seqx = new char[xlen + 1];
        secx = new char[xlen + 1];
        NewArray(&xa, xlen, 3);
        copy_chain_data(xa_vec[i], seqx_vec[i], secx_vec[i], xlen, xa, seqx, secx);

        ylen = ylen_vec[j];
        seqy = new char[ylen + 1];
        secy = new char[ylen + 1];
        NewArray(&ya, ylen, 3);
        copy_chain_data(ya_vec[j], seqy_vec[j], secy_vec[j], ylen, ya, seqy, secy);

        /* declare variable specific to this pair of TMalign */
        double t0[3], u0[3][3];
        double TM1, TM2;
        double TM3, TM4, TM5; // for a_opt, u_opt, d_opt
        double d0_0, TM_0;
        double d0A, d0B, d0u, d0a;
        double d0_out = 5.0;
        string seqM, seqxA, seqyA; // for output alignment
        double rmsd0 = 0.0;
        int L_ali; // Aligned length in standard_TMscore
        double Liden = 0;
        double TM_ali, rmsd_ali; // TMscore and rmsd in standard_TMscore
        int n_ali = 0;
        int n_ali8 = 0;
        vector<double> do_vec;

        int c;
        for (c = 0; c < sequence.size(); c++)
            sequence[c].clear();
        sequence.clear();
        sequence.push_back(seqxA_mat[i][j]);
        sequence.push_back(seqyA_mat[i][j]);

        /* entry function for structure alignment */
        TMalign_main(xa, ya, seqx, seqy, secx, secy,
                     t0, u0, TM1, TM2, TM3, TM4, TM5,
                     d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                     seqM, seqxA, seqyA, do_vec,
                     rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                     xlen, ylen, sequence, Lnorm_ass, d0_scale,
                     3, a_opt, u_opt, d_opt, fast_opt,
                     mol_vec1[i] + mol_vec2[j], -1, 0);

        for (ui = 0; ui < 3; ui++)
            for (uj = 0; uj < 3; uj++)
                ut_mat[i][ui * 3 + uj] = u0[ui][uj];
        for (uj = 0; uj < 3; uj++)
            ut_mat[i][9 + uj] = t0[uj];

        TM_vec.push_back(TM1);
        TM_vec.push_back(TM2);

        if (outfmt_opt < 2)
            output_results(
                xname.c_str(), yname.c_str(),
                chainID_list1[i], chainID_list2[j],
                xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5,
                rmsd0, d0_out, seqM.c_str(),
                seqxA.c_str(), seqyA.c_str(), Liden,
                n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
                d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u,
                "", outfmt_opt, ter_opt, false, split_opt,
                false, "", // o_opt, fname_super+chainID_list1[i],
                false, a_opt, u_opt, d_opt, mirror_opt,
                resi_vec1, resi_vec2);

        /* clean up */
        seqM.clear();
        seqxA.clear();
        seqyA.clear();

        delete[] seqy;
        delete[] secy;
        DeleteArray(&ya, ylen);

        delete[] seqx;
        delete[] secx;
        DeleteArray(&xa, xlen);
        do_vec.clear();
    }
    if (outfmt_opt == 2)
    {
        double TM = 0;
        for (i = 0; i < TM_vec.size(); i++)
            TM += TM_vec[i] * TM_vec[i];
        TM = sqrt(TM / TM_vec.size());
        string query_name = xname;
        string template_name = yname;
        for (i = 0; i < chain1_num; i++)
        {
            j = assign1_list[i];
            if (j < 0)
                continue;
            query_name += chainID_list1[i];
            template_name += chainID_list2[j];
        }
        printf("%s\t%s\t%.4f\n", query_name.c_str(), template_name.c_str(), TM);
        query_name.clear();
        template_name.clear();
    }

    if (m_opt)
        output_dock_rotation_matrix(fname_matrix.c_str(),
                                    xname_vec, yname_vec, ut_mat, assign1_list);

    if (o_opt)
        output_dock(chain1_list, ter_opt, split_opt, infmt1_opt,
                    atom_opt, mirror_opt, ut_mat, fname_super);

    /* clean up everything */
    vector<double>().swap(TM_vec);
    vector<string>().swap(xname_vec);
    vector<string>().swap(yname_vec);
    delete[] assign1_list;
    delete[] assign2_list;
    DeleteArray(&TMave_mat, chain1_num);
    DeleteArray(&ut_mat, chain1_num);
    vector<vector<string>>().swap(seqxA_mat);
    vector<vector<string>>().swap(seqM_mat);
    vector<vector<string>>().swap(seqyA_mat);
    vector<string>().swap(tmp_str_vec);

    vector<vector<vector<double>>>().swap(xa_vec); // structure of complex1
    vector<vector<vector<double>>>().swap(ya_vec); // structure of complex2
    vector<vector<char>>().swap(seqx_vec);         // sequence of complex1
    vector<vector<char>>().swap(seqy_vec);         // sequence of complex2
    vector<vector<char>>().swap(secx_vec);         // secondary structure of complex1
    vector<vector<char>>().swap(secy_vec);         // secondary structure of complex2
    mol_vec1.clear();                              // molecule type of complex1, RNA if >0
    mol_vec2.clear();                              // molecule type of complex2, RNA if >0
    vector<string>().swap(chainID_list1);          // list of chainID1
    vector<string>().swap(chainID_list2);          // list of chainID2
    xlen_vec.clear();                              // length of complex1
    ylen_vec.clear();                              // length of complex2
    return 1;
}

int mTMalign(string &xname, string &yname, const string &fname_super,
             const string &fname_matrix,
             vector<string> &sequence, double Lnorm_ass, const double d0_scale,
             const bool m_opt, const int i_opt, const int o_opt, const int a_opt,
             bool u_opt, const bool d_opt, const bool full_opt, const double TMcut,
             const int infmt_opt, const int ter_opt,
             const int split_opt, const int outfmt_opt, bool fast_opt,
             const int het_opt, const string &atom_opt, const bool autojustify,
             const string &mol_opt, const string &dir_opt, const int byresi_opt,
             const vector<string> &chain_list, const vector<string> &chain2parse,
             const vector<string> &model2parse, const bool se_opt)
{
    /* declare previously global variables */
    vector<vector<vector<double>>> a_vec;  // atomic structure
    vector<vector<vector<double>>> ua_vec; // unchanged atomic structure
    vector<vector<char>> seq_vec;          // sequence of complex
    vector<vector<char>> sec_vec;          // secondary structure of complex
    vector<int> mol_vec;                   // molecule type of complex1, RNA if >0
    vector<string> chainID_list;           // list of chainID
    vector<int> len_vec;                   // length of complex
    int i, j;                              // chain index
    int xlen, ylen;                        // chain length
    double **xa, **ya;                     // structure of single chain
    char *seqx, *seqy;                     // for the protein sequence
    char *secx, *secy;                     // for the secondary structure
    int len_aa, len_na;                    // total length of protein and RNA/DNA
    vector<string> resi_vec;               // residue index for chain

    /* parse chain list */
    parse_chain_list(chain_list, a_vec, seq_vec, sec_vec, mol_vec,
                     len_vec, chainID_list, ter_opt, split_opt, mol_opt, infmt_opt,
                     atom_opt, autojustify, false, het_opt, len_aa, len_na, o_opt,
                     resi_vec, chain2parse, model2parse);
    int chain_num = a_vec.size();
    if (chain_num <= 1)
        PrintErrorAndQuit("ERROR! <2 chains for multiple alignment");
    if (m_opt || o_opt)
        for (i = 0; i < chain_num; i++)
            ua_vec.push_back(a_vec[i]);
    int mol_type = 0;
    int total_len = 0;
    xlen = 0;
    for (i = 0; i < chain_num; i++)
    {
        if (len_vec[i] > xlen)
            xlen = len_vec[i];
        total_len += len_vec[i];
        mol_type += mol_vec[i];
    }
    if (!u_opt)
        Lnorm_ass = total_len / chain_num;
    u_opt = true;
    total_len -= xlen;
    if (total_len > 750)
        fast_opt = true;

    /* get all-against-all alignment */
    double **TMave_mat;
    NewArray(&TMave_mat, chain_num, chain_num);
    vector<string> tmp_str_vec(chain_num, "");
    vector<vector<string>> seqxA_mat(chain_num, tmp_str_vec);
    vector<vector<string>> seqyA_mat(chain_num, tmp_str_vec);
    for (i = 0; i < chain_num; i++)
        for (j = 0; j < chain_num; j++)
            TMave_mat[i][j] = 0;
    for (i = 0; i < chain_num; i++)
    {
        xlen = len_vec[i];
        if (xlen < 3)
            continue;
        seqx = new char[xlen + 1];
        secx = new char[xlen + 1];
        NewArray(&xa, xlen, 3);
        copy_chain_data(a_vec[i], seq_vec[i], sec_vec[i], xlen, xa, seqx, secx);
        seqxA_mat[i][i] = seqyA_mat[i][i] = (string)(seqx);
        for (j = i + 1; j < chain_num; j++)
        {
            ylen = len_vec[j];
            if (ylen < 3)
                continue;
            seqy = new char[ylen + 1];
            secy = new char[ylen + 1];
            NewArray(&ya, ylen, 3);
            copy_chain_data(a_vec[j], seq_vec[j], sec_vec[j], ylen, ya, seqy, secy);

            /* declare variable specific to this pair of TMalign */
            double t0[3], u0[3][3];
            double TM1, TM2;
            double TM3, TM4, TM5; // for a_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out = 5.0;
            string seqM, seqxA, seqyA; // for output alignment
            double rmsd0 = 0.0;
            int L_ali; // Aligned length in standard_TMscore
            double Liden = 0;
            double TM_ali, rmsd_ali; // TMscore and rmsd in standard_TMscore
            int n_ali = 0;
            int n_ali8 = 0;
            vector<double> do_vec;

            /* entry function for structure alignment */
            if (se_opt)
            {
                int *invmap = new int[ylen + 1];
                u0[0][0] = u0[1][1] = u0[2][2] = 1;
                u0[0][1] = u0[0][2] =
                    u0[1][0] = u0[1][2] =
                        u0[2][0] = u0[2][1] =
                            t0[0] = t0[1] = t0[2] = 0;
                se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA, do_vec,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        0, false, u_opt, false, mol_type, outfmt_opt, invmap);
                if (outfmt_opt >= 2)
                {
                    Liden = L_ali = 0;
                    int r1, r2;
                    for (r2 = 0; r2 < ylen; r2++)
                    {
                        r1 = invmap[r2];
                        if (r1 < 0)
                            continue;
                        L_ali += 1;
                        Liden += (seqx[r1] == seqy[r2]);
                    }
                }
                delete[] invmap;
            }
            else
                TMalign_main(xa, ya, seqx, seqy, secx, secy,
                             t0, u0, TM1, TM2, TM3, TM4, TM5,
                             d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                             seqM, seqxA, seqyA, do_vec,
                             rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                             xlen, ylen, sequence, Lnorm_ass, d0_scale,
                             0, false, u_opt, false, fast_opt,
                             mol_type, TMcut, 0);

            /* store result */
            TMave_mat[i][j] = TMave_mat[j][i] = TM4;
            seqxA_mat[i][j] = seqyA_mat[j][i] = seqxA;
            seqyA_mat[i][j] = seqxA_mat[j][i] = seqyA;
            // cout<<chain_list[i]<<':'<<chainID_list[i]
            //<<chain_list[j]<<':'<<chainID_list[j]<<"\tTM4="<<TM4<<endl;
            if (full_opt)
                output_results(
                    chain_list[i], chain_list[j], chainID_list[i], chainID_list[j],
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

            delete[] seqy;
            delete[] secy;
            DeleteArray(&ya, ylen);
            do_vec.clear();
        }

        delete[] seqx;
        delete[] secx;
        DeleteArray(&xa, xlen);
    }

    /* representative related variables */
    int r;
    int repr_idx = 0;
    vector<string> xname_vec;
    for (i = 0; i < chain_num; i++)
        xname_vec.push_back(
            chain_list[i].substr(dir_opt.size()) + chainID_list[i]);
    vector<string> yname_vec;
    double *TMave_list;
    TMave_list = new double[chain_num];
    int *assign_list;
    assign_list = new int[chain_num];
    vector<string> msa(ylen, ""); // row is position along msa; column is sequence

    int compare_num;
    double TM1_total, TM2_total;
    double TM3_total, TM4_total, TM5_total; // for a_opt, u_opt, d_opt
    double d0_0_total, TM_0_total;
    double d0A_total, d0B_total, d0u_total, d0a_total;
    double d0_out_total;
    double rmsd0_total;
    int L_ali_total; // Aligned length in standard_TMscore
    double Liden_total;
    double TM_ali_total, rmsd_ali_total; // TMscore and rmsd in standard_TMscore
    int n_ali_total;
    int n_ali8_total;
    int xlen_total, ylen_total;
    double TM4_total_max = 0;

    int max_iter = 5 - (int)(total_len / 200);
    if (max_iter < 2)
        max_iter = 2;
    int iter = 0;
    vector<double> TM_vec(chain_num, 0);
    vector<double> d0_vec(chain_num, 0);
    vector<double> seqID_vec(chain_num, 0);
    vector<vector<double>> TM_mat(chain_num, TM_vec);
    vector<vector<double>> d0_mat(chain_num, d0_vec);
    vector<vector<double>> seqID_mat(chain_num, seqID_vec);
    for (iter = 0; iter < max_iter; iter++)
    {
        /* select representative */
        for (j = 0; j < chain_num; j++)
            TMave_list[j] = 0;
        for (i = 0; i < chain_num; i++)
        {
            for (j = 0; j < chain_num; j++)
            {
                // cout<<'\t'<<setprecision(4)<<TMave_mat[i][j];
                TMave_list[j] += TMave_mat[i][j];
            }
            // cout<<'\t'<<chain_list[i]<<':'<<chainID_list[i]<<endl;
        }
        repr_idx = 0;
        double repr_TM = 0;
        for (j = 0; j < chain_num; j++)
        {
            // cout<<chain_list[j]<<'\t'<<len_vec[j]<<'\t'<<TMave_list[j]<<endl;
            if (TMave_list[j] < repr_TM)
                continue;
            repr_TM = TMave_list[j];
            repr_idx = j;
        }
        // cout<<"repr="<<repr_idx<<"; "<<chain_list[repr_idx]<<"; TM="<<repr_TM<<endl;

        /* superpose */
        yname = chain_list[repr_idx].substr(dir_opt.size()) + chainID_list[repr_idx];
        double **xt;
        vector<pair<double, int>> TM_pair_vec; // TM vs chain

        for (i = 0; i < chain_num; i++)
            assign_list[i] = -1;
        assign_list[repr_idx] = repr_idx;
        // ylen = len_vec[repr_idx];
        // seqy = new char[ylen+1];
        // secy = new char[ylen+1];
        // NewArray(&ya, ylen, 3);
        // copy_chain_data(a_vec[repr_idx],seq_vec[repr_idx],sec_vec[repr_idx], ylen,ya,seqy,secy);
        for (r = 0; r < sequence.size(); r++)
            sequence[r].clear();
        sequence.clear();
        sequence.push_back("");
        sequence.push_back("");
        for (i = 0; i < chain_num; i++)
        {
            yname_vec.push_back(yname);
            xlen = len_vec[i];
            if (i == repr_idx || xlen < 3)
                continue;
            TM_pair_vec.push_back(make_pair(-TMave_mat[i][repr_idx], i));
        }
        sort(TM_pair_vec.begin(), TM_pair_vec.end());

        int tm_idx;
        if (outfmt_opt < 0)
            cout << "#PDBchain1\tPDBchain2\tTM1\tTM2\t"
                 << "RMSD\tID1\tID2\tIDali\tL1\tL2\tLali" << endl;
        for (tm_idx = 0; tm_idx < TM_pair_vec.size(); tm_idx++)
        {
            i = TM_pair_vec[tm_idx].second;
            xlen = len_vec[i];
            seqx = new char[xlen + 1];
            secx = new char[xlen + 1];
            NewArray(&xa, xlen, 3);
            copy_chain_data(a_vec[i], seq_vec[i], sec_vec[i], xlen, xa, seqx, secx);

            double maxTM = TMave_mat[i][repr_idx];
            int maxj = repr_idx;
            for (j = 0; j < chain_num; j++)
            {
                if (i == j || assign_list[j] < 0 || TMave_mat[i][j] <= maxTM)
                    continue;
                maxj = j;
                maxTM = TMave_mat[i][j];
            }
            j = maxj;
            assign_list[i] = j;
            ylen = len_vec[j];
            seqy = new char[ylen + 1];
            secy = new char[ylen + 1];
            NewArray(&ya, ylen, 3);
            copy_chain_data(a_vec[j], seq_vec[j], sec_vec[j], ylen, ya, seqy, secy);

            sequence[0] = seqxA_mat[i][j];
            sequence[1] = seqyA_mat[i][j];
            // cout<<"tm_idx="<<tm_idx<<"\ti="<<i<<"\tj="<<j<<endl;
            // cout<<"superpose "<<xname_vec[i]<<" to "<<xname_vec[j]<<endl;

            /* declare variable specific to this pair of TMalign */
            double t0[3], u0[3][3];
            double TM1, TM2;
            double TM3, TM4, TM5; // for a_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out = 5.0;
            string seqM, seqxA, seqyA; // for output alignment
            double rmsd0 = 0.0;
            int L_ali; // Aligned length in standard_TMscore
            double Liden = 0;
            double TM_ali, rmsd_ali; // TMscore and rmsd in standard_TMscore
            int n_ali = 0;
            int n_ali8 = 0;
            vector<double> do_vec;

            /* entry function for structure alignment */
            if (se_opt)
            {
                int *invmap = new int[ylen + 1];
                u0[0][0] = u0[1][1] = u0[2][2] = 1;
                u0[0][1] = u0[0][2] =
                    u0[1][0] = u0[1][2] =
                        u0[2][0] = u0[2][1] =
                            t0[0] = t0[1] = t0[2] = 0;
                se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                        seqM, seqxA, seqyA, do_vec,
                        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        2, a_opt, u_opt, d_opt, mol_type, outfmt_opt, invmap);
                if (outfmt_opt >= 2)
                {
                    Liden = L_ali = 0;
                    int r1, r2;
                    for (r2 = 0; r2 < ylen; r2++)
                    {
                        r1 = invmap[r2];
                        if (r1 < 0)
                            continue;
                        L_ali += 1;
                        Liden += (seqx[r1] == seqy[r2]);
                    }
                }
                delete[] invmap;
            }
            else
                TMalign_main(xa, ya, seqx, seqy, secx, secy,
                             t0, u0, TM1, TM2, TM3, TM4, TM5,
                             d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                             seqM, seqxA, seqyA, do_vec,
                             rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                             xlen, ylen, sequence, Lnorm_ass, d0_scale,
                             2, a_opt, u_opt, d_opt, fast_opt, mol_type, -1, 0);

            if (outfmt_opt < 0)
                output_results(
                    xname_vec[i].c_str(), xname_vec[j].c_str(), "", "",
                    xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5,
                    rmsd0, d0_out, seqM.c_str(),
                    seqxA.c_str(), seqyA.c_str(), Liden,
                    n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
                    d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u,
                    "", 2, // outfmt_opt,
                    ter_opt, false, split_opt,
                    false, "", // o_opt, fname_super+chainID_list1[i],
                    false, a_opt, u_opt, d_opt, false,
                    resi_vec, resi_vec);

            NewArray(&xt, xlen, 3);
            do_rotation(xa, xt, xlen, t0, u0);
            for (r = 0; r < xlen; r++)
            {
                a_vec[i][r][0] = xt[r][0];
                a_vec[i][r][1] = xt[r][1];
                a_vec[i][r][2] = xt[r][2];
            }
            DeleteArray(&xt, xlen);

            /* clean up */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();
            sequence[0].clear();
            sequence[1].clear();

            delete[] seqx;
            delete[] secx;
            DeleteArray(&xa, xlen);

            delete[] seqy;
            delete[] secy;
            DeleteArray(&ya, ylen);
            do_vec.clear();
        }
        ylen = len_vec[repr_idx];
        seqy = new char[ylen + 1];
        secy = new char[ylen + 1];
        NewArray(&ya, ylen, 3);
        copy_chain_data(a_vec[repr_idx], seq_vec[repr_idx], sec_vec[repr_idx], ylen, ya, seqy, secy);

        /* recover alignment */
        int ylen_ext = ylen; // chain length
        double **ya_ext;     // structure of single chain
        char *seqy_ext;      // for the protein sequence
        char *secy_ext;      // for the secondary structure
        for (r = 0; r < msa.size(); r++)
            msa[r].clear();
        msa.clear();
        msa.assign(ylen, "");   // row is position along msa; column is sequence
        vector<string> msa_ext; // row is position along msa; column is sequence
        for (r = 0; r < ylen; r++)
            msa[r] = seqy[r];
        // for (r=0;r<msa.size();r++) cout<<"["<<r<<"]\t"<<msa[r]<<endl;
        // cout<<"start recover"<<endl;
        assign_list[repr_idx] = 0;
        for (tm_idx = 0; tm_idx < TM_pair_vec.size(); tm_idx++)
        {
            i = TM_pair_vec[tm_idx].second;
            assign_list[i] = tm_idx + 1;

            xlen = len_vec[i];
            seqx = new char[xlen + 1];
            secx = new char[xlen + 1];
            NewArray(&xa, xlen, 3);
            copy_chain_data(a_vec[i], seq_vec[i], sec_vec[i], xlen, xa, seqx, secx);

            /* declare variable specific to this pair of TMalign */
            double TM1, TM2;
            double TM3, TM4, TM5; // for a_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out = 5.0;
            string seqM, seqxA, seqyA; // for output alignment
            double rmsd0 = 0.0;
            int L_ali; // Aligned length in standard_TMscore
            double Liden = 0;
            double TM_ali, rmsd_ali; // TMscore and rmsd in standard_TMscore
            int n_ali = 0;
            int n_ali8 = 0;
            int *invmap = new int[ylen + 1];
            vector<double> do_vec;

            se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
                    do_vec, rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                    xlen, ylen, sequence, Lnorm_ass, d0_scale,
                    0, a_opt, u_opt, d_opt, mol_type, 1, invmap);

            int rx = 0, ry = 0;
            ylen_ext = seqxA.size();
            NewArray(&ya_ext, ylen_ext, 3);    // structure of single chain
            seqy_ext = new char[ylen_ext + 1]; // for the protein sequence
            secy_ext = new char[ylen_ext + 1]; // for the secondary structure
            string tmp_gap = "";
            for (r = 0; r < msa[0].size(); r++)
                tmp_gap += '-';
            for (r = msa_ext.size(); r < ylen_ext; r++)
                msa_ext.push_back("");
            // cout<<"x:"<<xname_vec[i]<<'\n'<<seqxA<<endl;
            // cout<<"y:"<<xname_vec[repr_idx]<<'\n'<<seqyA<<endl;
            for (r = 0; r < ylen_ext; r++)
            {
                if (seqyA[r] == '-')
                {
                    msa_ext[r] = tmp_gap + seqxA[r];
                    ya_ext[r][0] = xa[rx][0];
                    ya_ext[r][1] = xa[rx][1];
                    ya_ext[r][2] = xa[rx][2];
                    seqy_ext[r] = seqx[rx];
                    secy_ext[r] = secx[rx];
                }
                else
                {
                    msa_ext[r] = msa[ry] + seqxA[r];
                    ya_ext[r][0] = ya[ry][0];
                    ya_ext[r][1] = ya[ry][1];
                    ya_ext[r][2] = ya[ry][2];
                    seqy_ext[r] = seqy[ry];
                    secy_ext[r] = secy[ry];
                }
                rx += (seqxA[r] != '-');
                ry += (seqyA[r] != '-');
            }

            /* copy ya_ext to ya */
            delete[] seqy;
            delete[] secy;
            DeleteArray(&ya, ylen);

            ylen = ylen_ext;
            NewArray(&ya, ylen, 3);
            seqy = new char[ylen + 1];
            secy = new char[ylen + 1];
            for (r = 0; r < ylen; r++)
            {
                ya[r][0] = ya_ext[r][0];
                ya[r][1] = ya_ext[r][1];
                ya[r][2] = ya_ext[r][2];
                seqy[r] = seqy_ext[r];
                secy[r] = secy_ext[r];
            }
            for (r = 0; r < ylen; r++)
            {
                if (r < msa.size())
                    msa[r] = msa_ext[r];
                else
                    msa.push_back(msa_ext[r]);
            }
            // for (r=0;r<ylen_ext;r++) cout<<"["<<r<<"]\t"<<msa_ext[r]<<'\t'<<seqy[r]<<'\t'
            //<<ya[r][0]<<'\t'<<ya[r][1]<<'\t'<<ya[r][2]<<'\t'<<secy[r]<<endl;

            /* clean up */
            tmp_gap.clear();
            delete[] invmap;
            seqM.clear();
            seqxA.clear();
            seqyA.clear();

            delete[] seqx;
            delete[] secx;
            DeleteArray(&xa, xlen);

            delete[] seqy_ext;
            delete[] secy_ext;
            DeleteArray(&ya_ext, ylen_ext);
            do_vec.clear();
        }
        vector<string>().swap(msa_ext);
        vector<pair<double, int>>().swap(TM_pair_vec);
        for (i = 0; i < chain_num; i++)
        {
            tm_idx = assign_list[i];
            if (tm_idx < 0)
                continue;
            seqyA_mat[i][i] = "";
            for (r = 0; r < ylen; r++)
                seqyA_mat[i][i] += msa[r][tm_idx];
            seqxA_mat[i][i] = seqyA_mat[i][i];
            // cout<<xname_vec[i]<<'\t'<<seqxA_mat[i][i]<<endl;
        }
        for (i = 0; i < chain_num; i++)
        {
            if (assign_list[i] < 0)
                continue;
            string seqxA = seqxA_mat[i][i];
            for (j = 0; j < chain_num; j++)
            {
                if (i == j || assign_list[j] < 0)
                    continue;
                string seqyA = seqyA_mat[j][j];
                seqxA_mat[i][j] = seqyA_mat[i][j] = "";
                for (r = 0; r < ylen; r++)
                {
                    if (seqxA[r] == '-' && seqyA[r] == '-')
                        continue;
                    seqxA_mat[i][j] += seqxA[r];
                    seqyA_mat[i][j] += seqyA[r];
                }
                seqyA.clear();
            }
            seqxA.clear();
        }

        /* recover statistics such as TM-score */
        compare_num = 0;
        TM1_total = 0, TM2_total = 0;
        TM3_total = 0, TM4_total = 0, TM5_total = 0;
        d0_0_total = 0, TM_0_total = 0;
        d0A_total = 0, d0B_total = 0, d0u_total = 0, d0a_total = 0;
        d0_out_total = 0;
        rmsd0_total = 0.0;
        L_ali_total = 0;
        Liden_total = 0;
        TM_ali_total = 0, rmsd_ali_total = 0;
        n_ali_total = 0;
        n_ali8_total = 0;
        xlen_total = 0, ylen_total = 0;
        for (i = 0; i < chain_num; i++)
        {
            xlen = len_vec[i];
            if (xlen < 3)
                continue;
            seqx = new char[xlen + 1];
            secx = new char[xlen + 1];
            NewArray(&xa, xlen, 3);
            copy_chain_data(a_vec[i], seq_vec[i], sec_vec[i], xlen, xa, seqx, secx);
            for (j = i + 1; j < chain_num; j++)
            {
                ylen = len_vec[j];
                if (ylen < 3)
                    continue;
                compare_num++;
                seqy = new char[ylen + 1];
                secy = new char[ylen + 1];
                NewArray(&ya, ylen, 3);
                copy_chain_data(a_vec[j], seq_vec[j], sec_vec[j], ylen, ya, seqy, secy);
                sequence[0] = seqxA_mat[i][j];
                sequence[1] = seqyA_mat[i][j];

                /* declare variable specific to this pair of TMalign */
                double TM1, TM2;
                double TM3, TM4, TM5; // for a_opt, u_opt, d_opt
                double d0_0, TM_0;
                double d0A, d0B, d0u, d0a;
                double d0_out = 5.0;
                string seqM, seqxA, seqyA; // for output alignment
                double rmsd0 = 0.0;
                int L_ali = 0; // Aligned length in standard_TMscore
                double Liden = 0;
                double TM_ali, rmsd_ali; // TMscore and rmsd in standard_TMscore
                int n_ali = 0;
                int n_ali8 = 0;
                int *invmap = new int[ylen + 1];
                vector<double> do_vec;

                se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
                        do_vec, rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                        xlen, ylen, sequence, Lnorm_ass, d0_scale,
                        true, a_opt, u_opt, d_opt, mol_type, 1, invmap);

                if (xlen <= ylen)
                {
                    xlen_total += xlen;
                    ylen_total += ylen;
                    TM1_total += TM1;
                    TM2_total += TM2;
                    d0A_total += d0A;
                    d0B_total += d0B;
                }
                else
                {
                    xlen_total += ylen;
                    ylen_total += xlen;
                    TM1_total += TM2;
                    TM2_total += TM1;
                    d0A_total += d0B;
                    d0B_total += d0A;
                }
                TM_mat[i][j] = TM2;
                TM_mat[j][i] = TM1;
                d0_mat[i][j] = d0B;
                d0_mat[j][i] = d0A;
                seqID_mat[i][j] = 1. * Liden / xlen;
                seqID_mat[j][i] = 1. * Liden / ylen;

                TM3_total += TM3;
                TM4_total += TM4;
                TM5_total += TM5;
                d0_0_total += d0_0;
                TM_0_total += TM_0;
                d0u_total += d0u;
                d0_out_total += d0_out;
                rmsd0_total += rmsd0;
                L_ali_total += L_ali; // Aligned length in standard_TMscore
                Liden_total += Liden;
                TM_ali_total += TM_ali;
                rmsd_ali_total += rmsd_ali; // TMscore and rmsd in standard_TMscore
                n_ali_total += n_ali;
                n_ali8_total += n_ali8;

                /* clean up */
                delete[] invmap;
                seqM.clear();
                seqxA.clear();
                seqyA.clear();

                delete[] seqy;
                delete[] secy;
                DeleteArray(&ya, ylen);
                do_vec.clear();
            }
            delete[] seqx;
            delete[] secx;
            DeleteArray(&xa, xlen);
        }
        if (TM4_total <= TM4_total_max)
            break;
        TM4_total_max = TM4_total;
    }
    for (i = 0; i < chain_num; i++)
    {
        for (j = 0; j < chain_num; j++)
        {
            if (i == j)
                continue;
            TM_vec[i] += TM_mat[i][j];
            d0_vec[i] += d0_mat[i][j];
            seqID_vec[i] += seqID_mat[i][j];
        }
        TM_vec[i] /= (chain_num - 1);
        d0_vec[i] /= (chain_num - 1);
        seqID_vec[i] /= (chain_num - 1);
    }
    xlen_total /= compare_num;
    ylen_total /= compare_num;
    TM1_total /= compare_num;
    TM2_total /= compare_num;
    d0A_total /= compare_num;
    d0B_total /= compare_num;
    TM3_total /= compare_num;
    TM4_total /= compare_num;
    TM5_total /= compare_num;
    d0_0_total /= compare_num;
    TM_0_total /= compare_num;
    d0u_total /= compare_num;
    d0_out_total /= compare_num;
    rmsd0_total /= compare_num;
    L_ali_total /= compare_num;
    Liden_total /= compare_num;
    TM_ali_total /= compare_num;
    rmsd_ali_total /= compare_num;
    n_ali_total /= compare_num;
    n_ali8_total /= compare_num;
    xname = "shorter";
    yname = "longer";
    string seqM = "";
    string seqxA = "";
    string seqyA = "";
    double t0[3];
    double u0[3][3];
    stringstream buf;
    for (i = 0; i < chain_num; i++)
    {
        if (assign_list[i] < 0)
            continue;
        buf << ">" << xname_vec[i] << "\tL=" << len_vec[i]
            << "\td0=" << setiosflags(ios::fixed) << setprecision(2) << d0_vec[i]
            << "\tseqID=" << setiosflags(ios::fixed) << setprecision(3) << seqID_vec[i]
            << "\tTM-score=" << setiosflags(ios::fixed) << setprecision(5) << TM_vec[i];
        if (i == repr_idx)
            buf << "\t*";
        buf << '\n'
            << seqxA_mat[i][i] << endl;
    }
    seqM = buf.str();
    seqM = seqM.substr(0, seqM.size() - 1);
    buf.str(string());
    // MergeAlign(seqxA_mat,seqyA_mat,repr_idx,xname_vec,chain_num,seqM);
    if (outfmt_opt == 0)
        print_version();
    output_mTMalign_results(xname, yname, "", "",
                            xlen_total, ylen_total, t0, u0, TM1_total, TM2_total,
                            TM3_total, TM4_total, TM5_total, rmsd0_total, d0_out_total,
                            seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden_total,
                            n_ali8_total, L_ali_total, TM_ali_total, rmsd_ali_total,
                            TM_0_total, d0_0_total, d0A_total, d0B_total,
                            Lnorm_ass, d0_scale, d0a_total, d0u_total,
                            "", outfmt_opt, ter_opt, 0, split_opt, false,
                            "", false, a_opt, u_opt, d_opt, false,
                            resi_vec, resi_vec);

    if (m_opt || o_opt)
    {
        double **ut_mat; // rotation matrices for all-against-all alignment
        int ui, uj;
        double t[3], u[3][3];
        double rmsd;
        NewArray(&ut_mat, chain_num, 4 * 3);
        for (i = 0; i < chain_num; i++)
        {
            xlen = ylen = a_vec[i].size();
            NewArray(&xa, xlen, 3);
            NewArray(&ya, xlen, 3);
            for (r = 0; r < xlen; r++)
            {
                xa[r][0] = ua_vec[i][r][0];
                xa[r][1] = ua_vec[i][r][1];
                xa[r][2] = ua_vec[i][r][2];
                ya[r][0] = a_vec[i][r][0];
                ya[r][1] = a_vec[i][r][1];
                ya[r][2] = a_vec[i][r][2];
            }
            Kabsch(xa, ya, xlen, 1, &rmsd, t, u);
            for (ui = 0; ui < 3; ui++)
                for (uj = 0; uj < 3; uj++)
                    ut_mat[i][ui * 3 + uj] = u[ui][uj];
            for (uj = 0; uj < 3; uj++)
                ut_mat[i][9 + uj] = t[uj];
            DeleteArray(&xa, xlen);
            DeleteArray(&ya, xlen);
        }
        vector<vector<vector<double>>>().swap(ua_vec);

        if (m_opt)
        {
            assign_list[repr_idx] = -1;
            output_dock_rotation_matrix(fname_matrix.c_str(),
                                        xname_vec, yname_vec, ut_mat, assign_list);
        }

        // if (o_opt) output_dock(chain_list, ter_opt, split_opt,
        // infmt_opt, atom_opt, false, ut_mat, fname_super);
        if (o_opt)
            output_mTMalign_pymol(chain_list,
                                  infmt_opt, ut_mat, fname_super, o_opt);

        DeleteArray(&ut_mat, chain_num);
    }

    /* clean up */
    vector<string>().swap(msa);
    vector<string>().swap(tmp_str_vec);
    vector<vector<string>>().swap(seqxA_mat);
    vector<vector<string>>().swap(seqyA_mat);
    vector<string>().swap(xname_vec);
    vector<string>().swap(yname_vec);
    delete[] TMave_list;
    DeleteArray(&TMave_mat, chain_num);
    vector<vector<vector<double>>>().swap(a_vec); // structure of complex
    vector<vector<char>>().swap(seq_vec);         // sequence of complex
    vector<vector<char>>().swap(sec_vec);         // secondary structure of complex
    vector<int>().swap(mol_vec);                  // molecule type of complex1, RNA if >0
    vector<string>().swap(chainID_list);          // list of chainID
    vector<int>().swap(len_vec);                  // length of complex
    vector<double>().swap(TM_vec);
    vector<double>().swap(d0_vec);
    vector<double>().swap(seqID_vec);
    vector<vector<double>>().swap(TM_mat);
    vector<vector<double>>().swap(d0_mat);
    vector<vector<double>>().swap(seqID_mat);
    return 1;
}

/* sequence order independent alignment */
int SOIalign(string &xname, string &yname, const string &fname_super,
             const string &fname_lign, const string &fname_matrix,
             vector<string> &sequence, const double Lnorm_ass, const double d0_scale,
             const bool m_opt, const int i_opt, const int o_opt, const int a_opt,
             const bool u_opt, const bool d_opt, const double TMcut,
             const int infmt1_opt, const int infmt2_opt, const int ter_opt,
             const int split_opt, const int outfmt_opt, const bool fast_opt,
             const int cp_opt, const int mirror_opt, const int het_opt,
             const string &atom_opt, const bool autojustify, const string &mol_opt,
             const string &dir_opt, const string &dirpair_opt, const string &dir1_opt,
             const string &dir2_opt, const vector<string> &chain2parse1,
             const vector<string> &chain2parse2, const vector<string> &model2parse1,
             const vector<string> &model2parse2, const vector<string> &chain1_list,
             const vector<string> &chain2_list, const bool se_opt,
             const int closeK_opt, const int mm_opt)
{
    /* declare previously global variables */
    vector<vector<string>> PDB_lines1; // text of chain1
    vector<vector<string>> PDB_lines2; // text of chain2
    vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    vector<string> chainID_list1;      // list of chainID1
    vector<string> chainID_list2;      // list of chainID2
    int i, j;                          // file index
    int chain_i, chain_j;              // chain index
    int r;                             // residue index
    int xlen, ylen;                    // chain length
    int xchainnum, ychainnum;          // number of chains in a PDB file
    char *seqx, *seqy;                 // for the protein sequence
    char *secx, *secy;                 // for the secondary structure
    int **secx_bond;                   // boundary of secondary structure
    int **secy_bond;                   // boundary of secondary structure
    double **xa, **ya;                 // for input vectors xa[0...xlen-1][0..2] and
                                       // ya[0...ylen-1][0..2], in general,
                                       // ya is regarded as native structure
                                       // --> superpose xa onto ya
    double **xk, **yk;                 // k closest residues
    vector<string> resi_vec1;          // residue index for chain1
    vector<string> resi_vec2;          // residue index for chain2
    int read_resi = 0;                 // whether to read residue index
    if (o_opt)
        read_resi = 2;

    /* loop over file names */
    for (i = 0; i < chain1_list.size(); i++)
    {
        /* parse chain 1 */
        xname = chain1_list[i];
        xchainnum = get_PDB_lines(xname, PDB_lines1, chainID_list1, mol_vec1,
                                  ter_opt, infmt1_opt, atom_opt, autojustify, split_opt, het_opt,
                                  chain2parse1, model2parse1);
        if (!xchainnum)
        {
            cerr << "Warning! Cannot parse file: " << xname
                 << ". Chain number 0." << endl;
            continue;
        }
        for (chain_i = 0; chain_i < xchainnum; chain_i++)
        {
            xlen = PDB_lines1[chain_i].size();
            if (mol_opt == "RNA")
                mol_vec1[chain_i] = 1;
            else if (mol_opt == "protein")
                mol_vec1[chain_i] = -1;
            if (!xlen)
            {
                cerr << "Warning! Cannot parse file: " << xname
                     << ". Chain length 0." << endl;
                continue;
            }
            else if (xlen < 3)
            {
                cerr << "Sequence is too short <3!: " << xname << endl;
                continue;
            }
            NewArray(&xa, xlen, 3);
            if (closeK_opt >= 3)
                NewArray(&xk, xlen * closeK_opt, 3);
            seqx = new char[xlen + 1];
            secx = new char[xlen + 1];
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx,
                            resi_vec1, read_resi);
            if (mirror_opt)
                for (r = 0; r < xlen; r++)
                    xa[r][2] = -xa[r][2];
            if (mol_vec1[chain_i] > 0)
                make_sec(seqx, xa, xlen, secx, atom_opt);
            else
                make_sec(xa, xlen, secx); // secondary structure assignment
            if (closeK_opt >= 3)
                getCloseK(xa, xlen, closeK_opt, xk);
            if (mm_opt == 6)
            {
                NewArray(&secx_bond, xlen, 2);
                assign_sec_bond(secx_bond, secx, xlen);
            }

            for (j = (dir_opt.size() > 0) * (i + 1); j < chain2_list.size(); j++)
            {
                if (dirpair_opt.size() && i != j)
                    continue;
                /* parse chain 2 */
                if (PDB_lines2.size() == 0)
                {
                    yname = chain2_list[j];
                    ychainnum = get_PDB_lines(yname, PDB_lines2, chainID_list2,
                                              mol_vec2, ter_opt, infmt2_opt, atom_opt, autojustify,
                                              split_opt, het_opt, chain2parse2, model2parse2);
                    if (!ychainnum)
                    {
                        cerr << "Warning! Cannot parse file: " << yname
                             << ". Chain number 0." << endl;
                        continue;
                    }
                }
                for (chain_j = 0; chain_j < ychainnum; chain_j++)
                {
                    ylen = PDB_lines2[chain_j].size();
                    if (mol_opt == "RNA")
                        mol_vec2[chain_j] = 1;
                    else if (mol_opt == "protein")
                        mol_vec2[chain_j] = -1;
                    if (!ylen)
                    {
                        cerr << "Warning! Cannot parse file: " << yname
                             << ". Chain length 0." << endl;
                        continue;
                    }
                    else if (ylen < 3)
                    {
                        cerr << "Sequence is too short <3!: " << yname << endl;
                        continue;
                    }
                    NewArray(&ya, ylen, 3);
                    if (closeK_opt >= 3)
                        NewArray(&yk, ylen * closeK_opt, 3);
                    seqy = new char[ylen + 1];
                    secy = new char[ylen + 1];
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy,
                                    resi_vec2, read_resi);
                    if (mol_vec2[chain_j] > 0)
                        make_sec(seqy, ya, ylen, secy, atom_opt);
                    else
                        make_sec(ya, ylen, secy);
                    if (closeK_opt >= 3)
                        getCloseK(ya, ylen, closeK_opt, yk);
                    if (mm_opt == 6)
                    {
                        NewArray(&secy_bond, ylen, 2);
                        assign_sec_bond(secy_bond, secy, ylen);
                    }

                    /* declare variable specific to this pair of TMalign */
                    double t0[3], u0[3][3];
                    double TM1, TM2;
                    double TM3, TM4, TM5; // for a_opt, u_opt, d_opt
                    double d0_0, TM_0;
                    double d0A, d0B, d0u, d0a;
                    double d0_out = 5.0;
                    string seqM, seqxA, seqyA; // for output alignment
                    double rmsd0 = 0.0;
                    int L_ali; // Aligned length in standard_TMscore
                    double Liden = 0;
                    double TM_ali, rmsd_ali; // TMscore and rmsd in standard_TMscore
                    int n_ali = 0;
                    int n_ali8 = 0;
                    bool force_fast_opt = (getmin(xlen, ylen) > 1500) ? true : fast_opt;
                    int *invmap = new int[ylen + 1];
                    double *dist_list = new double[ylen + 1];

                    /* entry function for structure alignment */
                    if (se_opt)
                    {
                        u0[0][0] = u0[1][1] = u0[2][2] = 1;
                        u0[0][1] = u0[0][2] =
                            u0[1][0] = u0[1][2] =
                                u0[2][0] = u0[2][1] =
                                    t0[0] = t0[1] = t0[2] = 0;
                        soi_se_main(xa, ya, seqx, seqy, TM1, TM2, TM3, TM4, TM5,
                                    d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                                    seqM, seqxA, seqyA,
                                    rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                                    xlen, ylen, Lnorm_ass, d0_scale,
                                    i_opt, a_opt, u_opt, d_opt,
                                    mol_vec1[chain_i] + mol_vec2[chain_j],
                                    outfmt_opt, invmap, dist_list,
                                    secx_bond, secy_bond, mm_opt);
                        if (outfmt_opt >= 2)
                        {
                            Liden = L_ali = 0;
                            int r1, r2;
                            for (r2 = 0; r2 < ylen; r2++)
                            {
                                r1 = invmap[r2];
                                if (r1 < 0)
                                    continue;
                                L_ali += 1;
                                Liden += (seqx[r1] == seqy[r2]);
                            }
                        }
                    }
                    else
                        SOIalign_main(xa, ya, xk, yk, closeK_opt,
                                      seqx, seqy, secx, secy,
                                      t0, u0, TM1, TM2, TM3, TM4, TM5,
                                      d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                                      seqM, seqxA, seqyA, invmap,
                                      rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                                      xlen, ylen, sequence, Lnorm_ass, d0_scale,
                                      i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                                      mol_vec1[chain_i] + mol_vec2[chain_j], dist_list,
                                      secx_bond, secy_bond, mm_opt);

                    /* print result */
                    if (outfmt_opt == 0)
                        print_version();
                    output_results(
                        xname.substr(dir1_opt.size() + dir_opt.size() + dirpair_opt.size()),
                        yname.substr(dir2_opt.size() + dir_opt.size() + dirpair_opt.size()),
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5,
                        rmsd0, d0_out, seqM.c_str(),
                        seqxA.c_str(), seqyA.c_str(), Liden,
                        n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0,
                        d0A, d0B, Lnorm_ass, d0_scale, d0a, d0u,
                        (m_opt ? fname_matrix : "").c_str(),
                        outfmt_opt, ter_opt, false, split_opt, o_opt,
                        fname_super, i_opt, a_opt, u_opt, d_opt, mirror_opt,
                        resi_vec1, resi_vec2);
                    if (outfmt_opt <= 0)
                    {
                        cout << "###############\t###############\t#########" << endl;
                        cout << "#Aligned atom 1\tAligned atom 2 \tDistance#" << endl;
                        int r1, r2;
                        for (r2 = 0; r2 < ylen; r2++)
                        {
                            r1 = invmap[r2];
                            if (r1 < 0)
                                continue;
                            cout << PDB_lines1[chain_i][r1].substr(12, 15) << '\t'
                                 << PDB_lines2[chain_j][r2].substr(12, 15) << '\t'
                                 << setw(9) << setiosflags(ios::fixed) << setprecision(3)
                                 << dist_list[r2] << '\n';
                        }
                        cout << "###############\t###############\t#########" << endl;
                    }

                    /* Done! Free memory */
                    delete[] invmap;
                    delete[] dist_list;
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();
                    DeleteArray(&ya, ylen);
                    if (closeK_opt >= 3)
                        DeleteArray(&yk, ylen * closeK_opt);
                    delete[] seqy;
                    delete[] secy;
                    resi_vec2.clear();
                    if (mm_opt == 6)
                        DeleteArray(&secy_bond, ylen);
                } // chain_j
                if (chain2_list.size() > 1)
                {
                    yname.clear();
                    for (chain_j = 0; chain_j < ychainnum; chain_j++)
                        PDB_lines2[chain_j].clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
            PDB_lines1[chain_i].clear();
            DeleteArray(&xa, xlen);
            if (closeK_opt >= 3)
                DeleteArray(&xk, xlen * closeK_opt);
            delete[] seqx;
            delete[] secx;
            resi_vec1.clear();
            if (mm_opt == 6)
                DeleteArray(&secx_bond, xlen);
        } // chain_i
        xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
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
// Data structures and Helpers for flexalign unified pipeline
// =======================================================================

// Data structure to hold outputs of flexalign_main to avoid parameter clutter
struct FlexAlignResult
{
    double t0[3];
    double u0[3][3];
    vector<vector<double>> tu_vec;
    double TM1, TM2, TM3, TM4, TM5;
    double d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out;
    string seqM, seqxA, seqyA;
    vector<double> do_vec;
    double rmsd0, Liden, TM_ali, rmsd_ali;
    int L_ali, n_ali, n_ali8, hingeNum;

    FlexAlignResult() : TM1(-1.0), TM2(-1.0), TM3(-1.0), TM4(-1.0), TM5(-1.0),
                        d0_0(0.0), TM_0(0.0), d0A(0.0), d0B(0.0), d0u(0.0), d0a(0.0), d0_out(5.0),
                        rmsd0(0.0), Liden(0.0), TM_ali(0.0), rmsd_ali(0.0),
                        L_ali(0), n_ali(0), n_ali8(0), hingeNum(0)
    {
        for (int i = 0; i < 3; i++)
        {
            t0[i] = 0.0;
            for (int j = 0; j < 3; j++)
                u0[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
};

enum FlexAlignMode
{
    FLEX_STANDARD = 0,
    FLEX_BEST = 1,
    FLEX_FATCAT = 2
};

// Encapsulates the execution of flexalign_main and its fallback refinement logic
void execute_flexalign_with_fallback(
    double **xa, double **ya, char *seqx, char *seqy, char *secx, char *secy,
    int xlen, int ylen, vector<string> &sequence, const double Lnorm_ass, const double d0_scale,
    const int i_opt, const int a_opt, const bool u_opt, const bool d_opt, const bool force_fast_opt,
    const int mol_type, const int hinge_opt, const int ss_opt, FlexAlignResult &res)
{
    res.hingeNum = flexalign_main(
        xa, ya, seqx, seqy, secx, secy,
        res.t0, res.u0, res.tu_vec, res.TM1, res.TM2, res.TM3, res.TM4, res.TM5,
        res.d0_0, res.TM_0, res.d0A, res.d0B, res.d0u, res.d0a, res.d0_out,
        res.seqM, res.seqxA, res.seqyA, res.do_vec,
        res.rmsd0, res.L_ali, res.Liden, res.TM_ali, res.rmsd_ali, res.n_ali, res.n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt, force_fast_opt,
        mol_type, hinge_opt, ss_opt);

    // Fallback compensation when too few hinges are found
    if (hinge_opt && res.hingeNum <= 1 && res.n_ali8 < 0.6 * getmin(xlen, ylen))
    {
        FlexAlignResult res_h;
        res_h.tu_vec.push_back(res.tu_vec[0]);
        tu2t_u(res.tu_vec[0], res_h.t0, res_h.u0);

        res_h.hingeNum = flexalign_main(
            xa, ya, seqx, seqy, secx, secy,
            res_h.t0, res_h.u0, res_h.tu_vec,
            res_h.TM1, res_h.TM2, res_h.TM3, res_h.TM4, res_h.TM5,
            res_h.d0_0, res_h.TM_0, res.d0A, res.d0B, res.d0u, res.d0a, res_h.d0_out,
            res_h.seqM, res_h.seqxA, res_h.seqyA, res_h.do_vec,
            res_h.rmsd0, res_h.L_ali, res_h.Liden, res_h.TM_ali, res_h.rmsd_ali,
            res_h.n_ali, res_h.n_ali8,
            xlen, ylen, sequence, Lnorm_ass, d0_scale, i_opt,
            a_opt, u_opt, d_opt, force_fast_opt,
            mol_type, hinge_opt, ss_opt);

        double TM = (res.TM1 > res.TM2) ? res.TM1 : res.TM2;
        double TM_h = (res_h.TM1 > res_h.TM2) ? res_h.TM1 : res_h.TM2;
        if (TM_h > TM)
        {
            res = res_h; // Safely overwrite with the better refined results
        }
    }
}

// ==========================================
// FATCAT Core Algorithm (flexalign_fatcat_main)
// ==========================================
struct FATCAT_AFP
{
    int i, j, len;
    double score;
    double R[3][3];
    double t[3];
};

int flexalign_fatcat_main(double **xa, double **ya,
                          const char *seqx, const char *seqy, const char *secx, const char *secy,
                          double t0[3], double u0[3][3], std::vector<std::vector<double>> &tu_vec,
                          double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
                          double &d0_0, double &TM_0,
                          double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
                          std::string &seqM, std::string &seqxA, std::string &seqyA, std::vector<double> &do_vec,
                          double &rmsd0, int &L_ali, double &Liden,
                          double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
                          const int xlen, const int ylen,
                          const std::vector<std::string> sequence, const double Lnorm_ass,
                          const double d0_scale, const int i_opt, const int a_opt,
                          const bool u_opt, const bool d_opt, const bool fast_opt,
                          const int mol_type, const int hinge_opt, const int ss_opt,
                          int sparse_val = 0)
{
    // ==========================================
    // TRUE -mm 9 BASELINE (Defender)
    // Run full sequence without generate_bounds slicing!
    // This perfectly simulates FLEX_BEST (-mm 9) behavior.
    // ==========================================
    double best_global_max_TM = -1.0;
    std::vector<std::vector<double>> best_tu_vec;
    double best_t0[3], best_u0[3][3];
    double best_TM1 = 0.0, best_TM2 = 0.0, best_TM3 = 0.0, best_TM4 = 0.0, best_TM5 = 0.0;
    double best_rmsd0 = 0.0, best_Liden = 0.0, best_TM_ali = 0.0, best_rmsd_ali = 0.0;
    int best_L_ali = 0, best_n_ali = 0, best_n_ali8 = 0;
    std::string best_seqM = "", best_seqxA = "", best_seqyA = "";
    std::vector<double> best_do_vec;
    double best_d0A = 0.0, best_d0B = 0.0, best_d0a = 0.0, best_d0u = 0.0;

    bool force_fast_opt_global = (std::min(xlen, ylen) > 1500) ? true : fast_opt;
    std::vector<std::string> local_sequence = sequence;

    for (int cur_ss_opt = 0; cur_ss_opt <= 1; cur_ss_opt++)
    {
        FlexAlignResult base_res;
        // Pass full unbroken sequences directly to flexalign (identical to -mm 9)
        execute_flexalign_with_fallback(
            xa, ya, (char*)seqx, (char*)seqy, (char*)secx, (char*)secy,
            xlen, ylen, local_sequence, Lnorm_ass, d0_scale,
            i_opt, a_opt, u_opt, d_opt, force_fast_opt_global,
            mol_type, 9, cur_ss_opt, base_res); // -mm 9 explicitly uses 9 hinges

        double cur_max_TM = (base_res.TM1 > base_res.TM2) ? base_res.TM1 : base_res.TM2;
        if (cur_max_TM > best_global_max_TM)
        {
            best_global_max_TM = cur_max_TM;
            for (int a = 0; a < 3; a++) {
                best_t0[a] = base_res.t0[a];
                for (int b = 0; b < 3; b++) best_u0[a][b] = base_res.u0[a][b];
            }
            best_tu_vec = base_res.tu_vec;
            best_TM1 = base_res.TM1; best_TM2 = base_res.TM2; best_TM3 = base_res.TM3;
            best_TM4 = base_res.TM4; best_TM5 = base_res.TM5;
            best_rmsd0 = base_res.rmsd0; best_Liden = base_res.Liden;
            best_TM_ali = base_res.TM_ali; best_rmsd_ali = base_res.rmsd_ali;
            best_L_ali = base_res.L_ali; best_n_ali = base_res.n_ali; best_n_ali8 = base_res.n_ali8;
            best_seqM = base_res.seqM; best_seqxA = base_res.seqxA; best_seqyA = base_res.seqyA;
            best_do_vec = base_res.do_vec;
            best_d0A = base_res.d0A; best_d0B = base_res.d0B;
            best_d0a = base_res.d0a; best_d0u = base_res.d0u;
        }
    }

    // Early exit if the true -mm 9 baseline is already excellent
    if (best_global_max_TM >= 0.85)
    {
        TM1 = best_TM1; TM2 = best_TM2; TM3 = best_TM3; TM4 = best_TM4; TM5 = best_TM5;
        rmsd0 = best_rmsd0; Liden = best_Liden; TM_ali = best_TM_ali; rmsd_ali = best_rmsd_ali;
        L_ali = best_L_ali; n_ali = best_n_ali; n_ali8 = best_n_ali8;
        seqM = best_seqM; seqxA = best_seqxA; seqyA = best_seqyA;
        do_vec = best_do_vec; tu_vec = best_tu_vec;
        d0A = best_d0A; d0B = best_d0B; d0a = best_d0a; d0u = best_d0u;
        for (int a = 0; a < 3; a++) {
            t0[a] = best_t0[a];
            for (int b = 0; b < 3; b++) u0[a][b] = best_u0[a][b];
        }
        return tu_vec.size();
    }

    // ==========================================
    // Proceed to FATCAT sliced bounds logic...
    // ==========================================

    // FATCAT base parameters
    int fragLen = 8;
    double resScore = 3.0;
    double gap_ext = -0.5;
    double disCut = 5.0;
    double disSmooth = 4.0;
    double twist_pen = -25.0;
    int max_gap = 40;
    double max_penalty = -5.0;
    int misCut = 2 * fragLen;
    int maxGapFrag = fragLen + max_gap;
    double afp_dis_cut = fragLen * fragLen * (disCut * disCut);
    int max_twists = hinge_opt;

    // OPTIMIZATION 1: Precompute local intra-protein distance matrices
    int max_dist_window = max_gap + 2 * fragLen + 1;
    std::vector<std::vector<double>> disTable1(xlen, std::vector<double>(max_dist_window, 0.0));
    std::vector<std::vector<double>> disTable2(ylen, std::vector<double>(max_dist_window, 0.0));

    for (int i = 0; i < xlen; i++)
    {
        for (int j = i; j < std::min(xlen, i + max_dist_window); j++)
            disTable1[i][j - i] = std::sqrt(dist(xa[i], xa[j]));
    }
    for (int i = 0; i < ylen; i++)
    {
        for (int j = i; j < std::min(ylen, i + max_dist_window); j++)
            disTable2[i][j - i] = std::sqrt(dist(ya[i], ya[j]));
    }

    // Wrapper for generating bounds
    auto generate_bounds = [&](double cur_rmsdCut, double cur_badRmsd, double cur_local_badRmsd) -> std::pair<std::vector<int>, std::vector<int>>
    {
        // Step 1: Extract initial AFPs in batches
        std::vector<FATCAT_AFP> initial_afps;
        int step = sparse_val + 1;

        double r1_static[8][3], r2_static[8][3];
        double *r1[8], *r2[8];
        for (int k = 0; k < 8; k++)
        {
            r1[k] = r1_static[k];
            r2[k] = r2_static[k];
        }

        for (int i = 0; i <= xlen - fragLen; i += step)
        {
            for (int j = 0; j <= ylen - fragLen; j += step)
            {
                int d3_term = std::min(i, j) + std::min(xlen - (i + fragLen - 1), ylen - (j + fragLen)) + fragLen;
                if (d3_term < 0.3 * std::min(xlen, ylen))
                    continue;

                double dist1 = disTable1[i][fragLen - 1]; 
                double dist2 = disTable2[j][fragLen - 1];

                if (std::fabs(dist1 - dist2) > 2.0 * cur_rmsdCut)
                    continue;

                for (int k = 0; k < fragLen; k++)
                {
                    r1[k][0] = xa[i + k][0]; r1[k][1] = xa[i + k][1]; r1[k][2] = xa[i + k][2];
                    r2[k][0] = ya[j + k][0]; r2[k][1] = ya[j + k][1]; r2[k][2] = ya[j + k][2];
                }

                double rms_sum_sq, t_tmp[3], u_tmp[3][3];
                Kabsch(r1, r2, fragLen, 0, &rms_sum_sq, t_tmp, u_tmp);
                double rmsd_tmp = std::sqrt(rms_sum_sq / fragLen);

                if (rmsd_tmp < cur_rmsdCut)
                {
                    FATCAT_AFP afp;
                    afp.i = i; afp.j = j; afp.len = fragLen;
                    afp.score = resScore * fragLen * (1.0 - (rmsd_tmp / cur_badRmsd) * (rmsd_tmp / cur_badRmsd));
                    for (int a = 0; a < 3; a++)
                    {
                        afp.t[a] = t_tmp[a];
                        for (int b = 0; b < 3; b++) afp.R[a][b] = u_tmp[a][b];
                    }
                    initial_afps.push_back(afp);
                }
            }
        }

        // Step 2: Merge diagonal AFPs
        int max_diagonal_idx = xlen + ylen + 1;
        std::vector<std::vector<FATCAT_AFP>> diagonals(max_diagonal_idx);
        for (size_t k = 0; k < initial_afps.size(); k++)
        {
            diagonals[initial_afps[k].i - initial_afps[k].j + ylen].push_back(initial_afps[k]);
        }

        std::vector<FATCAT_AFP> merged_afps;
        int max_merge_len = std::min(xlen, ylen);
        double **r1_merge, **r2_merge;
        NewArray(&r1_merge, max_merge_len, 3);
        NewArray(&r2_merge, max_merge_len, 3);

        for (int d = 0; d < max_diagonal_idx; d++)
        {
            if (diagonals[d].empty())
                continue;
            std::vector<FATCAT_AFP> &group = diagonals[d];

            std::sort(group.begin(), group.end(), [](const FATCAT_AFP &a, const FATCAT_AFP &b) { return a.i < b.i; });

            int n_group = group.size();
            std::vector<bool> invalid(n_group, false);
            for (int idx = 0; idx < n_group; idx++)
            {
                if (invalid[idx])
                    continue;
                FATCAT_AFP curr = group[idx];
                for (int nxt_idx = idx + 1; nxt_idx < n_group; nxt_idx++)
                {
                    FATCAT_AFP nxt = group[nxt_idx];
                    if (nxt.i > curr.i + curr.len)
                        break;

                    if (nxt.i + nxt.len > curr.i + curr.len)
                    {
                        int new_len = (nxt.i + nxt.len) - curr.i;

                        for (int k = 0; k < new_len; k++)
                        {
                            r1_merge[k][0] = xa[curr.i + k][0]; r1_merge[k][1] = xa[curr.i + k][1]; r1_merge[k][2] = xa[curr.i + k][2];
                            r2_merge[k][0] = ya[curr.j + k][0]; r2_merge[k][1] = ya[curr.j + k][1]; r2_merge[k][2] = ya[curr.j + k][2];
                        }

                        double rms_sum_sq, t_tmp[3], u_tmp[3][3];
                        Kabsch(r1_merge, r2_merge, new_len, 0, &rms_sum_sq, t_tmp, u_tmp);
                        double rmsd_tmp = std::sqrt(rms_sum_sq / new_len);

                        if (rmsd_tmp < cur_rmsdCut)
                        {
                            curr.len = new_len;
                            for (int a = 0; a < 3; a++)
                            {
                                curr.t[a] = t_tmp[a];
                                for (int b = 0; b < 3; b++) curr.R[a][b] = u_tmp[a][b];
                            }
                            curr.score = resScore * new_len * (1.0 - (rmsd_tmp / cur_badRmsd) * (rmsd_tmp / cur_badRmsd));
                            invalid[nxt_idx] = true;
                        }
                    }
                }
                merged_afps.push_back(curr);
            }
        }
        DeleteArray(&r1_merge, max_merge_len);
        DeleteArray(&r2_merge, max_merge_len);

        std::sort(merged_afps.begin(), merged_afps.end(), [](const FATCAT_AFP &a, const FATCAT_AFP &b) {
            if (a.i == b.i) return a.j < b.j;
            return a.i < b.i; 
        });

        int n_afps = merged_afps.size();
        std::vector<int> ret_b1, ret_b2;
        if (n_afps == 0)
            return std::make_pair(ret_b1, ret_b2);

        // Step 3 & 4: Dual Dynamic Programming and Domain Splitting
        std::vector<int> afp_aft_index(xlen * ylen, -1);
        std::vector<int> afp_bef_index(xlen * ylen, -1);

        std::vector<std::vector<std::pair<int, int>>> i_to_j(xlen);
        for (int m = 0; m < n_afps; m++)
        {
            i_to_j[merged_afps[m].i].push_back(std::make_pair(merged_afps[m].j, m));
        }

        for (int i_val = 0; i_val < xlen; i_val++)
        {
            if (i_to_j[i_val].empty()) continue;
            for (size_t p = 0; p < i_to_j[i_val].size(); p++)
            {
                int j_val = i_to_j[i_val][p].first;
                afp_aft_index[i_val * ylen + j_val] = i_to_j[i_val][p].second;
                afp_bef_index[i_val * ylen + j_val] = i_to_j[i_val][p].second;
            }
            int curr_bef = -1;
            for (int j_val = 0; j_val < ylen; j_val++)
            {
                if (afp_bef_index[i_val * ylen + j_val] != -1) curr_bef = afp_bef_index[i_val * ylen + j_val];
                else afp_bef_index[i_val * ylen + j_val] = curr_bef;
            }
            int curr_aft = -1;
            for (int j_val = ylen - 1; j_val >= 0; j_val--)
            {
                if (afp_aft_index[i_val * ylen + j_val] != -1) curr_aft = afp_aft_index[i_val * ylen + j_val];
                else afp_aft_index[i_val * ylen + j_val] = curr_aft;
            }
        }

        auto get_dvar = [&](const FATCAT_AFP &prv, const FATCAT_AFP &curr) -> double
        {
            double rms_sq = 0;
            for (int i_idx = 0; i_idx < fragLen; i_idx++)
            {
                for (int j_idx = 0; j_idx < fragLen; j_idx++)
                {
                    double dist1, dist2;
                    int idx1_a = curr.i + i_idx, idx1_b = prv.i + j_idx;
                    if (idx1_a >= idx1_b) dist1 = disTable1[idx1_b][idx1_a - idx1_b];
                    else dist1 = disTable1[idx1_a][idx1_b - idx1_a];

                    int idx2_a = curr.j + i_idx, idx2_b = prv.j + j_idx;
                    if (idx2_a >= idx2_b) dist2 = disTable2[idx2_b][idx2_a - idx2_b];
                    else dist2 = disTable2[idx2_a][idx2_b - idx2_a];

                    rms_sq += (dist1 - dist2) * (dist1 - dist2);
                }
            }
            if (rms_sq > afp_dis_cut) return 1e9;
            return std::sqrt(rms_sq / (fragLen * fragLen));
        };

        auto calc_block_rmsd = [&](const std::vector<FATCAT_AFP> &afp_list) -> double
        {
            std::vector<int> r1, r2;
            for (size_t a = 0; a < afp_list.size(); a++)
            {
                for (int l = 0; l < afp_list[a].len; l++)
                {
                    r1.push_back(afp_list[a].i + l);
                    r2.push_back(afp_list[a].j + l);
                }
            }
            int n = r1.size();
            if (n < 3) return 0.0;
            double **p1; NewArray(&p1, n, 3);
            double **p2; NewArray(&p2, n, 3);
            for (int i = 0; i < n; i++)
            {
                p1[i][0] = xa[r1[i]][0]; p1[i][1] = xa[r1[i]][1]; p1[i][2] = xa[r1[i]][2];
                p2[i][0] = ya[r2[i]][0]; p2[i][1] = ya[r2[i]][1]; p2[i][2] = ya[r2[i]][2];
            }
            double rms_sq_sum, t_tmp[3], u_tmp[3][3];
            Kabsch(p1, p2, n, 0, &rms_sq_sum, t_tmp, u_tmp);
            DeleteArray(&p1, n); DeleteArray(&p2, n);
            return std::sqrt(rms_sq_sum / n);
        };

        struct Region { int s1, e1, s2, e2; };

        std::vector<double> sco(n_afps);
        std::vector<int> twi(n_afps, 0);
        std::vector<int> pre(n_afps, -1);
        for (int m = 0; m < n_afps; m++) sco[m] = merged_afps[m].score;

        for (int m = 0; m < n_afps; m++)
        {
            int curr_i = merged_afps[m].i;
            int curr_j = merged_afps[m].j;
            int a3 = curr_i - fragLen;
            int a2 = std::max(0, a3 - misCut);
            int a1 = std::max(0, curr_i - maxGapFrag);
            int b3 = curr_j - fragLen;
            int b2 = std::max(0, b3 - misCut);
            int b1 = std::max(0, curr_j - maxGapFrag);

            std::vector<int> valid_prevs;
            for (int st = 0; st < 2; st++)
            {
                int a_s, a_e, b_s, b_e;
                if (st == 0)
                {
                    a_s = std::max(a1, 0); a_e = std::min(a3, xlen - 1);
                    b_s = std::max(b2, 0); b_e = std::min(b3, ylen - 1);
                }
                else
                {
                    a_s = std::max(a2, 0); a_e = std::min(a3, xlen - 1);
                    b_s = std::max(b1, 0); b_e = std::min(b2 - 1, ylen - 1);
                }

                if (b_s >= ylen || b_e < 0) continue;
                for (int prev_i = a_s; prev_i <= a_e; prev_i++)
                {
                    int s1 = afp_aft_index[prev_i * ylen + b_s];
                    int s2 = afp_bef_index[prev_i * ylen + b_e];
                    if (s1 != -1 && s2 != -1 && s1 <= s2)
                        for (int s = s1; s <= s2; s++) valid_prevs.push_back(s);
                }
            }

            double curr_sco = merged_afps[m].score;
            for (size_t v = 0; v < valid_prevs.size(); v++)
            {
                int prev = valid_prevs[v];
                int prev_twi = twi[prev];
                if (prev_twi > max_twists) continue;

                int gap_i = curr_i - (merged_afps[prev].i + merged_afps[prev].len);
                int gap_j = curr_j - (merged_afps[prev].j + merged_afps[prev].len);
                int m_gap = std::max(gap_i, gap_j);

                double gp = 0.0;
                int m_mis = 0;
                if (gap_i < 0 || gap_j < 0) m_mis = (gap_i < gap_j) ? -gap_i : -gap_j;
                gp = gap_ext * m_mis;
                if (m_gap > 0) gp += gap_ext * m_gap;
                if (gp < max_penalty) gp = max_penalty;

                double rms_sq = 0;
                for (int k = 0; k < fragLen; k++)
                {
                    for (int l = 0; l < fragLen; l++)
                    {
                        double dist1, dist2;
                        int idx1_a = curr_i + k, idx1_b = merged_afps[prev].i + l;
                        if (idx1_a >= idx1_b) dist1 = disTable1[idx1_b][idx1_a - idx1_b];
                        else dist1 = disTable1[idx1_a][idx1_b - idx1_a];

                        int idx2_a = curr_j + k, idx2_b = merged_afps[prev].j + l;
                        if (idx2_a >= idx2_b) dist2 = disTable2[idx2_b][idx2_a - idx2_b];
                        else dist2 = disTable2[idx2_a][idx2_b - idx2_a];

                        rms_sq += (dist1 - dist2) * (dist1 - dist2);
                    }
                }

                double tp = 0.0;
                int is_twist = 0;
                if (rms_sq >= afp_dis_cut)
                {
                    tp = twist_pen; is_twist = 1;
                }
                else
                {
                    double dvar = std::sqrt(rms_sq / (fragLen * fragLen));
                    if (dvar > disCut - disSmooth) tp = twist_pen * std::sqrt((dvar - disCut + disSmooth) / disSmooth);
                }

                if (prev_twi + is_twist > max_twists) continue;

                double stmp = sco[prev] + curr_sco + tp + gp;
                if (stmp > sco[m])
                {
                    sco[m] = stmp; pre[m] = prev; twi[m] = prev_twi + is_twist;
                }
            }
        }

        int best_m = 0;
        for (int m = 1; m < n_afps; m++)
            if (sco[m] > sco[best_m]) best_m = m;

        std::vector<int> path;
        int curr_m = best_m;
        while (curr_m != -1)
        {
            path.push_back(curr_m); curr_m = pre[curr_m];
        }
        std::reverse(path.begin(), path.end());

        if (path.empty()) return std::make_pair(ret_b1, ret_b2);

        struct Block
        {
            std::vector<FATCAT_AFP> afps;
            std::vector<double> dvars;
        };
        std::vector<Block> blocks;
        Block curr_block;
        curr_block.afps.push_back(merged_afps[path[0]]);
        curr_block.dvars.push_back(0.0);

        for (size_t k = 1; k < path.size(); k++)
        {
            FATCAT_AFP curr = merged_afps[path[k]];
            FATCAT_AFP prv = merged_afps[path[k - 1]];
            double dvar = get_dvar(prv, curr);

            if (dvar >= disCut)
            {
                blocks.push_back(curr_block);
                curr_block.afps.clear(); curr_block.dvars.clear();
                curr_block.afps.push_back(curr); curr_block.dvars.push_back(0.0);
            }
            else
            {
                curr_block.afps.push_back(curr); curr_block.dvars.push_back(dvar);
            }
        }
        if (!curr_block.afps.empty()) blocks.push_back(curr_block);

        bool splitted = true;
        while (splitted && blocks.size() < (size_t)(max_twists + 1))
        {
            splitted = false;
            double max_rmsd = 0.0;
            int target_b = -1;

            for (size_t b = 0; b < blocks.size(); b++)
            {
                if (blocks[b].afps.size() > 2)
                {
                    double cur_rmsd = calc_block_rmsd(blocks[b].afps);
                    if (cur_rmsd > max_rmsd) { max_rmsd = cur_rmsd; target_b = b; }
                }
            }

            if (max_rmsd >= cur_local_badRmsd && target_b != -1)
            {
                double max_t = 0; int cut_idx = 0;
                for (size_t i = 1; i < blocks[target_b].afps.size(); i++)
                {
                    if (blocks[target_b].dvars[i] > max_t) { max_t = blocks[target_b].dvars[i]; cut_idx = i; }
                }

                if (cut_idx > 0)
                {
                    Block right_blk;
                    right_blk.afps.assign(blocks[target_b].afps.begin() + cut_idx, blocks[target_b].afps.end());
                    right_blk.dvars.assign(blocks[target_b].dvars.begin() + cut_idx, blocks[target_b].dvars.end());
                    right_blk.dvars[0] = 0.0;
                    blocks[target_b].afps.erase(blocks[target_b].afps.begin() + cut_idx, blocks[target_b].afps.end());
                    blocks[target_b].dvars.erase(blocks[target_b].dvars.begin() + cut_idx, blocks[target_b].dvars.end());
                    blocks.insert(blocks.begin() + target_b + 1, right_blk);
                    splitted = true;
                }
            }
        }

        for (int b = 0; b < (int)blocks.size(); b++)
        {
            if (blocks[b].afps.size() <= 1)
            {
                int e1 = (b < (int)blocks.size() - 1) ? blocks[b + 1].afps.front().i : xlen;
                int e2 = (b < (int)blocks.size() - 1) ? blocks[b + 1].afps.front().j : ylen;
                int b1 = (b > 0) ? blocks[b - 1].afps.back().i + blocks[b - 1].afps.back().len : 0;
                int b2 = (b > 0) ? blocks[b - 1].afps.back().j + blocks[b - 1].afps.back().len : 0;
                int span = std::min(e1 - b1, e2 - b2);
                if (span < 2 * fragLen) { blocks.erase(blocks.begin() + b); b--; }
            }
        }

        bool merged = true;
        while (merged && blocks.size() > 1)
        {
            merged = false;
            double min_rmsd = 1e9;
            int min_b = -1;
            for (size_t b = 0; b < blocks.size() - 1; b++)
            {
                std::vector<FATCAT_AFP> temp_merged = blocks[b].afps;
                temp_merged.insert(temp_merged.end(), blocks[b + 1].afps.begin(), blocks[b + 1].afps.end());
                double cur_rmsd = calc_block_rmsd(temp_merged);
                if (cur_rmsd < min_rmsd) { min_rmsd = cur_rmsd; min_b = b; }
            }

            if (min_rmsd < cur_local_badRmsd && min_b != -1)
            {
                blocks[min_b].afps.insert(blocks[min_b].afps.end(), blocks[min_b + 1].afps.begin(), blocks[min_b + 1].afps.end());
                blocks.erase(blocks.begin() + min_b + 1);
                merged = true;
            }
        }

        std::vector<Region> real_blocks;
        int last_i = 0, last_j = 0;
        for (size_t b = 0; b < blocks.size(); b++)
        {
            int b_s1 = -1, b_e1 = -1, b_s2 = -1, b_e2 = -1;
            for (size_t a = 0; a < blocks[b].afps.size(); a++)
            {
                FATCAT_AFP afp = blocks[b].afps[a];
                int skip = std::max(std::max(last_i - afp.i, last_j - afp.j), 0);
                if (skip >= afp.len) continue;

                int eff_i = afp.i + skip; int eff_j = afp.j + skip; int eff_L = afp.len - skip;
                if (b_s1 == -1) { b_s1 = eff_i; b_s2 = eff_j; }
                b_e1 = eff_i + eff_L; b_e2 = eff_j + eff_L;
                last_i = b_e1; last_j = b_e2;
            }
            if (b_s1 != -1)
            {
                if (b_e1 - b_s1 >= 4 && b_e2 - b_s2 >= 4) {
                    Region r = {b_s1, b_e1, b_s2, b_e2};
                    real_blocks.push_back(r);
                }
            }
        }

        if (real_blocks.empty()) return std::make_pair(ret_b1, ret_b2);

        ret_b1.push_back(0); ret_b2.push_back(0);
        for (size_t k = 0; k < real_blocks.size() - 1; k++)
        {
            ret_b1.push_back((real_blocks[k].e1 + real_blocks[k + 1].s1) / 2);
            ret_b2.push_back((real_blocks[k].e2 + real_blocks[k + 1].s2) / 2);
        }
        ret_b1.push_back(xlen); ret_b2.push_back(ylen);

        return std::make_pair(ret_b1, ret_b2);
    };

    auto bounds_fatcat = generate_bounds(3.0, 4.0, 4.0);
    auto bounds_strict = generate_bounds(2.0, 3.0, 2.0);

    std::vector<std::pair<std::vector<int>, std::vector<int>>> all_bounds;
    all_bounds.push_back(bounds_fatcat);
    if (bounds_strict.first != bounds_fatcat.first || bounds_strict.second != bounds_fatcat.second)
    {
        all_bounds.push_back(bounds_strict);
    }

    // Loop through both bound sets, updating best_global_max_TM if we beat the -mm 9 defender
    for (size_t b_idx = 0; b_idx < all_bounds.size(); b_idx++)
    {
        std::vector<int>& bounds1 = all_bounds[b_idx].first;
        std::vector<int>& bounds2 = all_bounds[b_idx].second;

        // Skip if only one interval (block) is generated, as the full unbroken sequence
        // has already been processed by the baseline (-mm 9) above.
        if (bounds1.size() <= 2) continue;

        // ================== DEBUG START ==================
        // Output the interval mapping for the current boundary set
        // std::cout << "\n[DEBUG] --- Region Mapping Table ---" << std::endl;
        // std::cout << "[DEBUG] Mode: " << (b_idx == 0 ? "FATCAT Bounds" : "Strict Bounds") << std::endl;
        // std::cout << "[DEBUG] Total Blocks: " << (bounds1.size() - 1) << std::endl;
        
        // for (size_t k = 0; k < bounds1.size() - 1; k++)
        // {
        //     std::cout << "[DEBUG] Block " << (k + 1) << ": "
        //               << "Chain1 [" << bounds1[k] << " -> " << bounds1[k + 1] << "]  <==>  "
        //               << "Chain2 [" << bounds2[k] << " -> " << bounds2[k + 1] << "]" 
        //               << std::endl;
        // }
        // std::cout << "[DEBUG] ----------------------------\n" << std::endl;
        // =================== DEBUG END ===================

        // Step 5: Iteratively align each block
        std::string cur_global_seqM = "", cur_global_seqxA = "", cur_global_seqyA = "";
        cur_global_seqM.reserve(xlen + ylen + max_gap);
        cur_global_seqxA.reserve(xlen + ylen + max_gap);
        cur_global_seqyA.reserve(xlen + ylen + max_gap);

        std::vector<std::vector<double>> cur_tu_vec;
        std::vector<int> cur_global_res_tu(xlen, -1);

        for (size_t k = 0; k < bounds1.size() - 1; k++)
        {
            int x_s = bounds1[k], x_e = bounds1[k + 1];
            int y_s = bounds2[k], y_e = bounds2[k + 1];
            int L1_sub = x_e - x_s;
            int L2_sub = y_e - y_s;

            if (L1_sub < 3 || L2_sub < 3)
            {
                for (int i = 0; i < L1_sub; i++)
                {
                    cur_global_seqxA += seqx[x_s + i]; cur_global_seqyA += '-'; cur_global_seqM += ' ';
                }
                for (int i = 0; i < L2_sub; i++)
                {
                    cur_global_seqxA += '-'; cur_global_seqyA += seqy[y_s + i]; cur_global_seqM += ' ';
                }
                continue;
            }

            double **xa_sub, **ya_sub;
            NewArray(&xa_sub, L1_sub, 3);
            NewArray(&ya_sub, L2_sub, 3);
            char *seqx_sub = new char[L1_sub + 1];
            char *seqy_sub = new char[L2_sub + 1];
            char *secx_sub = new char[L1_sub + 1];
            char *secy_sub = new char[L2_sub + 1];

            for (int i = 0; i < L1_sub; i++)
            {
                xa_sub[i][0] = xa[x_s + i][0]; xa_sub[i][1] = xa[x_s + i][1]; xa_sub[i][2] = xa[x_s + i][2];
                seqx_sub[i] = seqx[x_s + i]; secx_sub[i] = secx[x_s + i];
            }
            seqx_sub[L1_sub] = '\0'; secx_sub[L1_sub] = '\0';

            for (int i = 0; i < L2_sub; i++)
            {
                ya_sub[i][0] = ya[y_s + i][0]; ya_sub[i][1] = ya[y_s + i][1]; ya_sub[i][2] = ya[y_s + i][2];
                seqy_sub[i] = seqy[y_s + i]; secy_sub[i] = secy[y_s + i];
            }
            seqy_sub[L2_sub] = '\0'; secy_sub[L2_sub] = '\0';

            double t0_best[3], u0_best[3][3];
            double TM_best_max = -1.0;
            std::string seqM_best, seqxA_best, seqyA_best;
            std::vector<std::vector<double>> tu_vec_best;

            bool force_fast_opt = (std::min(L1_sub, L2_sub) > 1500) ? true : fast_opt;

            // Set hinge_opt to 0 if the block length is less than 2 * fragLen
            int local_hinge_opt = (std::min(L1_sub, L2_sub) < 2 * fragLen) ? 0 : hinge_opt;

            for (int cur_ss_opt = 0; cur_ss_opt <= 1; cur_ss_opt++)
            {
                FlexAlignResult cur_res;
                execute_flexalign_with_fallback(
                    xa_sub, ya_sub, seqx_sub, seqy_sub, secx_sub, secy_sub,
                    L1_sub, L2_sub, local_sequence, Lnorm_ass, d0_scale,
                    i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                    mol_type, local_hinge_opt, cur_ss_opt, cur_res);

                double cur_max_TM = (cur_res.TM1 > cur_res.TM2) ? cur_res.TM1 : cur_res.TM2;
                if (cur_max_TM > TM_best_max)
                {
                    TM_best_max = cur_max_TM;
                    for (int a = 0; a < 3; a++)
                    {
                        t0_best[a] = cur_res.t0[a];
                        for (int b = 0; b < 3; b++) u0_best[a][b] = cur_res.u0[a][b];
                    }
                    seqM_best = cur_res.seqM; seqxA_best = cur_res.seqxA; seqyA_best = cur_res.seqyA;
                    tu_vec_best = cur_res.tu_vec;
                }
            }

            if (TM_best_max <= 0)
            {
                for (int i = 0; i < L1_sub; i++)
                {
                    cur_global_seqxA += seqx_sub[i]; cur_global_seqyA += '-'; cur_global_seqM += ' ';
                }
                for (int i = 0; i < L2_sub; i++)
                {
                    cur_global_seqxA += '-'; cur_global_seqyA += seqy_sub[i]; cur_global_seqM += ' ';
                }
                DeleteArray(&xa_sub, L1_sub); DeleteArray(&ya_sub, L2_sub);
                delete[] seqx_sub; delete[] seqy_sub; delete[] secx_sub; delete[] secy_sub;
                continue;
            }

            if (tu_vec_best.empty())
            {
                std::vector<double> tu_tmp(12);
                t_u2tu(t0_best, u0_best, tu_tmp);
                tu_vec_best.push_back(tu_tmp);
            }

            int base_tu_idx = cur_tu_vec.size();
            for (size_t m = 0; m < tu_vec_best.size(); m++) cur_tu_vec.push_back(tu_vec_best[m]);

            int rx = x_s;
            int current_global_idx = base_tu_idx;

            for (size_t i = 0; i < seqxA_best.length(); i++)
            {
                char c = seqM_best[i];

                if (c != ' ' && c != '.' && c != ':')
                {
                    int local_hinge_idx = -1;
                    if (c >= '0' && c <= '9') local_hinge_idx = c - '0';
                    else if (c >= 'a' && c <= 'z') local_hinge_idx = c - 'a' + 10;
                    else if (c >= 'A' && c <= 'Z') local_hinge_idx = c - 'A' + 36;
                    if (local_hinge_idx >= 0 && local_hinge_idx < tu_vec_best.size()) current_global_idx = base_tu_idx + local_hinge_idx;
                }

                if (seqxA_best[i] != '-')
                {
                    cur_global_res_tu[rx] = current_global_idx;
                    rx++;
                }

                if (seqxA_best[i] != '-' && seqyA_best[i] != '-')
                {
                    if (c != ' ' && c != '.' && c != ':')
                    {
                        char global_c;
                        if (current_global_idx < 10) global_c = '0' + current_global_idx;
                        else if (current_global_idx < 36) global_c = 'a' + (current_global_idx - 10);
                        else if (current_global_idx < 62) global_c = 'A' + (current_global_idx - 36);
                        else global_c = '*';
                        seqM_best[i] = global_c;
                    }
                    else
                    {
                        seqM_best[i] = c;
                    }
                }
                else { seqM_best[i] = ' '; }
            }

            cur_global_seqM += seqM_best; cur_global_seqxA += seqxA_best; cur_global_seqyA += seqyA_best;

            DeleteArray(&xa_sub, L1_sub); DeleteArray(&ya_sub, L2_sub);
            delete[] seqx_sub; delete[] seqy_sub; delete[] secx_sub; delete[] secy_sub;
        }

        // Step 6: Recalculate global metrics correctly for current DP boundary
        double cur_d0A = 1.24 * std::pow(ylen * 1.0 - 15.0, 1.0 / 3.0) - 1.8;
        if (cur_d0A < 0.5) cur_d0A = 0.5;
        double cur_d0B = 1.24 * std::pow(xlen * 1.0 - 15.0, 1.0 / 3.0) - 1.8;
        if (cur_d0B < 0.5) cur_d0B = 0.5;
        double cur_d0a = 1.24 * std::pow((xlen + ylen) * 0.5 - 15.0, 1.0 / 3.0) - 1.8;
        if (cur_d0a < 0.5) cur_d0a = 0.5;
        
        double cur_d0u = 0.0;
        if (u_opt)
        {
            cur_d0u = 1.24 * std::pow(Lnorm_ass - 15.0, 1.0 / 3.0) - 1.8;
            if (cur_d0u < 0.5) cur_d0u = 0.5;
        }

        double cur_TM1 = 0.0, cur_TM2 = 0.0, cur_TM3 = 0.0, cur_TM4 = 0.0, cur_TM5 = 0.0;
        double cur_rmsd0 = 0.0, cur_Liden = 0.0;
        int cur_n_ali8 = 0, cur_n_ali = 0;
        std::vector<double> cur_do_vec;

        int i_res = 0, j_res = 0;
        for (size_t r = 0; r < cur_global_seqxA.length(); r++)
        {
            bool x_valid = (cur_global_seqxA[r] != '-');
            bool y_valid = (cur_global_seqyA[r] != '-');

            if (x_valid && y_valid)
            {
                int matrix_idx = cur_global_res_tu[i_res];
                if (matrix_idx >= 0 && matrix_idx < cur_tu_vec.size())
                {
                    double t_k[3], u_k[3][3];
                    tu2t_u(cur_tu_vec[matrix_idx], t_k, u_k);

                    double x_rot[3];
                    transform(t_k, u_k, xa[i_res], x_rot);
                    double dist2 = dist(x_rot, ya[j_res]);
                    double d = std::sqrt(dist2);

                    cur_TM2 += 1.0 / (1.0 + dist2 / (cur_d0B * cur_d0B));
                    cur_TM1 += 1.0 / (1.0 + dist2 / (cur_d0A * cur_d0A));
                    if (a_opt) cur_TM3 += 1.0 / (1.0 + dist2 / (cur_d0a * cur_d0a));
                    if (u_opt) cur_TM4 += 1.0 / (1.0 + dist2 / (cur_d0u * cur_d0u));
                    if (d_opt) cur_TM5 += 1.0 / (1.0 + dist2 / (d0_scale * d0_scale));

                    cur_n_ali++;
                    cur_do_vec.push_back(d);

                    if (d <= d0_out)
                    {
                        cur_rmsd0 += dist2;
                        cur_n_ali8++;
                        if (seqx[i_res] == seqy[j_res]) cur_Liden += 1.0;
                    }
                }
                else { cur_do_vec.push_back(-1); }
            }
            else { cur_do_vec.push_back(-1); }

            if (x_valid) i_res++;
            if (y_valid) j_res++;
        }

        // Normalize TM-scores
        cur_TM2 /= xlen;
        cur_TM1 /= ylen;
        if (a_opt) cur_TM3 /= (xlen + ylen) * 0.5;
        if (u_opt) cur_TM4 /= Lnorm_ass;
        if (d_opt) cur_TM5 /= ylen;
        if (cur_n_ali8 > 0) cur_rmsd0 = std::sqrt(cur_rmsd0 / cur_n_ali8);
        else cur_rmsd0 = 0.0;

        // Compare against the -mm 9 defender!
        double cur_global_max_TM = (cur_TM1 > cur_TM2) ? cur_TM1 : cur_TM2;

        if (cur_global_max_TM > best_global_max_TM)
        {
            best_global_max_TM = cur_global_max_TM;
            best_tu_vec = cur_tu_vec;
            best_TM1 = cur_TM1; best_TM2 = cur_TM2; best_TM3 = cur_TM3; best_TM4 = cur_TM4; best_TM5 = cur_TM5;
            best_rmsd0 = cur_rmsd0; best_Liden = cur_Liden; best_TM_ali = cur_TM1; best_rmsd_ali = cur_rmsd0;
            best_L_ali = cur_n_ali; best_n_ali = cur_n_ali; best_n_ali8 = cur_n_ali8;
            best_seqM = cur_global_seqM; best_seqxA = cur_global_seqxA; best_seqyA = cur_global_seqyA;
            best_do_vec = cur_do_vec;
            best_d0A = cur_d0A; best_d0B = cur_d0B; best_d0a = cur_d0a; best_d0u = cur_d0u;

            if (!best_tu_vec.empty()) {
                tu2t_u(best_tu_vec[0], best_t0, best_u0);
            }
        }
    }

    // Safety check
    if (best_global_max_TM < 0) return 0;

    // Output best values back to the reference parameters
    TM1 = best_TM1; TM2 = best_TM2; TM3 = best_TM3; TM4 = best_TM4; TM5 = best_TM5;
    rmsd0 = best_rmsd0; Liden = best_Liden; TM_ali = best_TM_ali; rmsd_ali = best_rmsd_ali;
    L_ali = best_L_ali; n_ali = best_n_ali; n_ali8 = best_n_ali8;
    seqM = best_seqM; seqxA = best_seqxA; seqyA = best_seqyA;
    do_vec = best_do_vec; tu_vec = best_tu_vec;
    d0A = best_d0A; d0B = best_d0B; d0a = best_d0a; d0u = best_d0u;

    for (int a = 0; a < 3; a++) {
        t0[a] = best_t0[a];
        for (int b = 0; b < 3; b++) u0[a][b] = best_u0[a][b];
    }

    return tu_vec.size();
}

// Unified engine replacing flexalign, flexalign_best, and flexalign_fatcat
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
                      FlexAlignMode mode = FLEX_STANDARD)
{
    vector<vector<string>> PDB_lines1;
    vector<vector<string>> PDB_lines2;
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
                    if (mode == FLEX_FATCAT)
                    {
                        FlexAlignResult fatcat_res;
                        bool force_fast_opt = (getmin(xlen, ylen) > 1500) ? true : fast_opt;

                        fatcat_res.hingeNum = flexalign_fatcat_main(
                            xa, ya, seqx, seqy, secx, secy,
                            fatcat_res.t0, fatcat_res.u0, fatcat_res.tu_vec,
                            fatcat_res.TM1, fatcat_res.TM2, fatcat_res.TM3, fatcat_res.TM4, fatcat_res.TM5,
                            fatcat_res.d0_0, fatcat_res.TM_0,
                            fatcat_res.d0A, fatcat_res.d0B, fatcat_res.d0u, fatcat_res.d0a, fatcat_res.d0_out,
                            fatcat_res.seqM, fatcat_res.seqxA, fatcat_res.seqyA, fatcat_res.do_vec,
                            fatcat_res.rmsd0, fatcat_res.L_ali, fatcat_res.Liden,
                            fatcat_res.TM_ali, fatcat_res.rmsd_ali, fatcat_res.n_ali, fatcat_res.n_ali8,
                            xlen, ylen, sequence, Lnorm_ass, d0_scale,
                            i_opt, a_opt, u_opt, d_opt, force_fast_opt,
                            mol_vec1[chain_i] + mol_vec2[chain_j], hinge_opt, ss_opt, 0 /* sparse_val */
                        );

                        if (outfmt_opt == 0)
                            print_version();
                        output_flexalign_results(
                            xname.substr(dir1_opt.size() + dir_opt.size() + dirpair_opt.size()),
                            yname.substr(dir2_opt.size() + dir_opt.size() + dirpair_opt.size()),
                            chainID_list1[chain_i], chainID_list2[chain_j],
                            xlen, ylen, fatcat_res.t0, fatcat_res.u0, fatcat_res.tu_vec,
                            fatcat_res.TM1, fatcat_res.TM2, fatcat_res.TM3, fatcat_res.TM4, fatcat_res.TM5,
                            fatcat_res.rmsd0, fatcat_res.d0_out, fatcat_res.seqM.c_str(),
                            fatcat_res.seqxA.c_str(), fatcat_res.seqyA.c_str(), fatcat_res.Liden,
                            fatcat_res.n_ali8, fatcat_res.L_ali, fatcat_res.TM_ali, fatcat_res.rmsd_ali,
                            fatcat_res.TM_0, fatcat_res.d0_0,
                            fatcat_res.d0A, fatcat_res.d0B, Lnorm_ass, d0_scale, fatcat_res.d0a, fatcat_res.d0u,
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

int flexalign(string &xname, string &yname, const string &fname_super, const string &fname_lign, const string &fname_matrix, vector<string> &sequence, const double Lnorm_ass, const double d0_scale, const bool m_opt, const int i_opt, const int o_opt, const int a_opt, const bool u_opt, const bool d_opt, const double TMcut, const int infmt1_opt, const int infmt2_opt, const int ter_opt, const int split_opt, const int outfmt_opt, const bool fast_opt, const int mirror_opt, const int het_opt, const string &atom_opt, const bool autojustify, const string &mol_opt, const string &dir_opt, const string &dirpair_opt, const string &dir1_opt, const string &dir2_opt, const vector<string> &chain2parse1, const vector<string> &chain2parse2, const vector<string> &model2parse1, const vector<string> &model2parse2, const int byresi_opt, const vector<string> &chain1_list, const vector<string> &chain2_list, const int hinge_opt, const int ss_opt)
{
    return flexalign_unified(xname, yname, fname_super, fname_lign, fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt, split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt, atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt, dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2, byresi_opt, chain1_list, chain2_list, hinge_opt, ss_opt, FLEX_STANDARD);
}

int flexalign_best(string &xname, string &yname, const string &fname_super, const string &fname_lign, const string &fname_matrix, vector<string> &sequence, const double Lnorm_ass, const double d0_scale, const bool m_opt, const int i_opt, const int o_opt, const int a_opt, const bool u_opt, const bool d_opt, const double TMcut, const int infmt1_opt, const int infmt2_opt, const int ter_opt, const int split_opt, const int outfmt_opt, const bool fast_opt, const int mirror_opt, const int het_opt, const string &atom_opt, const bool autojustify, const string &mol_opt, const string &dir_opt, const string &dirpair_opt, const string &dir1_opt, const string &dir2_opt, const vector<string> &chain2parse1, const vector<string> &chain2parse2, const vector<string> &model2parse1, const vector<string> &model2parse2, const int byresi_opt, const vector<string> &chain1_list, const vector<string> &chain2_list, const int hinge_opt)
{
    return flexalign_unified(xname, yname, fname_super, fname_lign, fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt, split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt, atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt, dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2, byresi_opt, chain1_list, chain2_list, hinge_opt, 0 /* ss_opt is ignored in BEST mode */, FLEX_BEST);
}

int flexalign_fatcat(string &xname, string &yname, const string &fname_super, const string &fname_lign, const string &fname_matrix, vector<string> &sequence, const double Lnorm_ass, const double d0_scale, const bool m_opt, const int i_opt, const int o_opt, const int a_opt, const bool u_opt, const bool d_opt, const double TMcut, const int infmt1_opt, const int infmt2_opt, const int ter_opt, const int split_opt, const int outfmt_opt, const bool fast_opt, const int mirror_opt, const int het_opt, const string &atom_opt, const bool autojustify, const string &mol_opt, const string &dir_opt, const string &dirpair_opt, const string &dir1_opt, const string &dir2_opt, const vector<string> &chain2parse1, const vector<string> &chain2parse2, const vector<string> &model2parse1, const vector<string> &model2parse2, const int byresi_opt, const vector<string> &chain1_list, const vector<string> &chain2_list, const int hinge_opt)
{
    return flexalign_unified(xname, yname, fname_super, fname_lign, fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt, split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt, atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt, dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2, byresi_opt, chain1_list, chain2_list, hinge_opt, 0 /* ss_opt ignore */, FLEX_FATCAT);
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
    bool do_opt = false; // flag for -do, output distance of i-th aligned pair

    bool full_opt = false; // do not show chain level alignment
    double TMcut = -1;
    bool se_opt = false;
    int infmt1_opt = -1;        // PDB or PDBx/mmCIF format for chain_1
    int infmt2_opt = -1;        // PDB or PDBx/mmCIF format for chain_2
    int ter_opt = -1;           // default change to 2 (END, or different chainID)
    int split_opt = -1;         // default change to 2 (split each chains)
    int outfmt_opt = 0;         // set -outfmt to full output
    bool fast_opt = false;      // flags for -fast, fTM-align algorithm
    int cp_opt = 0;             // do not check circular permutation
    int closeK_opt = -1;        // number of atoms for SOI initial alignment.
                                // 5 and 0 for -mm 5 and 6
    int hinge_opt = 9;          // maximum number of hinge allowed for flexible
    int mirror_opt = 0;         // do not align mirror
    int het_opt = 0;            // do not read HETATM residues
    int mm_opt = 0;             // do not perform MM-align
    string atom_opt = "auto";   // use C alpha atom for protein and C3' for RNA
    string mol_opt = "auto";    // auto-detect the molecule type as protein/RNA
    string suffix_opt = "";     // set -suffix to empty
    string dir_opt = "";        // set -dir to empty
    string dirpair_opt = "";    // set -dirpair to empty
    string dir1_opt = "";       // set -dir1 to empty
    string dir2_opt = "";       // set -dir2 to empty
    string chainmapfile = "";   // chain mapping between two complexes
    int byresi_opt = 0;         // set -byresi to 0
    vector<string> chain1_list; // only when -dir1 is set
    vector<string> chain2_list; // only when -dir2 is set
    vector<string> chain2parse1;
    vector<string> chain2parse2;
    vector<string> model2parse1;
    vector<string> model2parse2;
    vector<pair<string, string>> chain_pair_list; // only when -dirpair is set

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
        else if (!strcmp(argv[i], "-full"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -full");
            if (!strcmp(argv[i + 1], "T"))
                full_opt = true;
            else if (!strcmp(argv[i + 1], "F"))
                full_opt = false;
            else
                PrintErrorAndQuit("-full must be T or F");
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
        else if (!strcmp(argv[i], "-closeK"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -closeK");
            closeK_opt = atoi(argv[i + 1]);
            i++;
        }
        else if (!strcmp(argv[i], "-hinge"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -hinge");
            hinge_opt = atoi(argv[i + 1]);
            i++;
        }
        else if (!strcmp(argv[i], "-v"))
        {
            v_opt = true;
        }
        else if (!strcmp(argv[i], "-do"))
        {
            do_opt = true;
        }
        else if (!strcmp(argv[i], "-h"))
        {
            h_opt = true;
        }
        else if (!strcmp(argv[i], "-i"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -i");
            if (i_opt == 3)
                PrintErrorAndQuit("ERROR! -i and -I cannot be used together");
            fname_lign = argv[i + 1];
            i_opt = 1;
            i++;
        }
        else if (!strcmp(argv[i], "-I"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -I");
            if (i_opt == 1)
                PrintErrorAndQuit("ERROR! -I and -i cannot be used together");
            fname_lign = argv[i + 1];
            i_opt = 3;
            i++;
        }
        else if (!strcmp(argv[i], "-chainmap"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -chainmap");
            chainmapfile = argv[i + 1];
            i++;
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
        else if (!strcmp(argv[i], "-se"))
        {
            se_opt = true;
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
        else if (!strcmp(argv[i], "-TMcut"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -TMcut");
            TMcut = atof(argv[i + 1]);
            i++;
        }
        else if (!strcmp(argv[i], "-byresi") ||
                 !strcmp(argv[i], "-tmscore") ||
                 !strcmp(argv[i], "-TMscore"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -byresi");
            byresi_opt = atoi(argv[i + 1]);
            i++;
        }
        else if (!strcmp(argv[i], "-seq"))
        {
            byresi_opt = 5;
        }
        else if (!strcmp(argv[i], "-cp"))
        {
            mm_opt = 3;
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
        else if (!strcmp(argv[i], "-mm"))
        {
            if (i >= (argc - 1))
                PrintErrorAndQuit("ERROR! Missing value for -mm");
            mm_opt = atoi(argv[i + 1]);
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
    if (byresi_opt != 0)
    {
        if (i_opt)
            PrintErrorAndQuit("-TMscore >=1 cannot be used with -i or -I");
        if (byresi_opt < 0 || byresi_opt > 7)
            PrintErrorAndQuit("-TMscore can only be 0 to 7");
        if ((byresi_opt == 2 || byresi_opt == 3 || byresi_opt == 6) && ter_opt >= 2)
            PrintErrorAndQuit("-TMscore 2 and 6 must be used with -ter <=1");
    }
    // if (split_opt==1 && ter_opt!=0)
    // PrintErrorAndQuit("-split 1 should be used with -ter 0");
    // else if (split_opt==2 && ter_opt!=0 && ter_opt!=1)
    // PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (split_opt < 0)
        if (byresi_opt == 2 || byresi_opt == 3)
            split_opt = 0;
        else
            split_opt = 2;
    else if (split_opt > 2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");

    if (mm_opt == 3)
    {
        cp_opt = true;
        mm_opt = 0;
    }
    if (cp_opt && i_opt)
        PrintErrorAndQuit("-mm 3 cannot be used with -i or -I");

    if (mirror_opt && het_opt != 1)
        cerr << "WARNING! -mirror was not used with -het 1. "
             << "D amino acids may not be correctly aligned." << endl;

    if (ter_opt < 0)
    {
        if (mm_opt == 1 || mm_opt == 2 || byresi_opt == 2 || byresi_opt == 3 ||
            byresi_opt == 6 || byresi_opt == 7)
            ter_opt = 1;
        else
            ter_opt = 2;
    }

    if (mm_opt)
    {
        if (i_opt) PrintErrorAndQuit("-mm cannot be used with -i or -I");
        if (u_opt) PrintErrorAndQuit("-mm cannot be used with -u or -L");
        //if (cp_opt) PrintErrorAndQuit("-mm cannot be used with -cp");
        if (dir_opt.size() && mm_opt==2) PrintErrorAndQuit("-mm 2 cannot be used with -dir");
        if (byresi_opt) PrintErrorAndQuit("-mm cannot be used with -byresi");
        if (ter_opt>=2 && (mm_opt==1 || mm_opt==2)) PrintErrorAndQuit("-mm 1 or 2 must be used with -ter 0 or -ter 1");
        if (mm_opt==4 && (yname.size() || dir2_opt.size()))
            cerr<<"WARNING! structure_2 is ignored for -mm 4"<<endl;
        if (dirpair_opt.size() && (mm_opt==2 || mm_opt==4))
            PrintErrorAndQuit("-mm 2 or 4 cannot be used with -dirpair");
    }
    else if (full_opt)
        PrintErrorAndQuit("-full can only be used with -mm");

    if (o_opt && ter_opt <= 1 && split_opt == 2)
    {
        if (mm_opt && o_opt == 2)
            cerr << "WARNING! -mm may generate incorrect"
                 << " RasMol output due to limitations in PDB file format. "
                 << "When -mm is used, -o is recommended over -rasmol" << endl;
        else if (mm_opt == 0)
            cerr << "WARNING! Only the superposition of the"
                 << " last aligned structure pair will be generated" << endl;
    }

    if (closeK_opt < 0)
    {
        if (mm_opt == 5)
            closeK_opt = 5;
        else
            closeK_opt = 0;
    }

    if (mm_opt >= 7 && hinge_opt >= 10)
        PrintErrorAndQuit("ERROR! -hinge must be <10");

    if (chainmapfile.size() && mm_opt != 1)
        PrintErrorAndQuit("ERROR! -chainmap must be used with -mm 1");

    /* read initial alignment file from 'align.txt' */
    if (i_opt)
        read_user_alignment(sequence, fname_lign, i_opt);

    if (byresi_opt == 6 || byresi_opt == 7)
        mm_opt = 1;
    else if (byresi_opt)
        i_opt = 3;

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
        if (mm_opt == 2)
            cout << "#Query\tTemplate\tTM" << endl;
        else
            cout << "#PDBchain1\tPDBchain2\tTM1\tTM2\t"
                 << "RMSD\tID1\tID2\tIDali\tL1\tL2\tLali" << endl;
    }

    /* real alignment. entry functions are MMalign_main and
     * TMalign_main */
    if (mm_opt==0) TMalign(xname, yname, fname_super, fname_lign, fname_matrix,
        sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt,
        u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
        split_opt, outfmt_opt, fast_opt, cp_opt, mirror_opt, het_opt,
        atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt,
        dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2,
        byresi_opt, chain1_list, chain2_list, se_opt, do_opt);
    else if (mm_opt==1)
    { 
        if (dir_opt.size()>0 || dir1_opt.size()>0 || dir2_opt.size()>0)
        {
            for (int ii=0; ii<chain1_list.size(); ii++)
            {
                xname = chain1_list[ii];
                vector<string> tmp_vec1(1, xname);
                for (int jj=0; jj<chain2_list.size(); jj++)
                {
                    if (dir_opt.size()>0 && jj<=ii) continue;
                    yname = chain2_list[jj];
                    vector<string> tmp_vec2(1, yname);
                    MMalign(xname, yname, fname_super,
                        fname_lign, fname_matrix, sequence, d0_scale, m_opt, o_opt,
                        a_opt, d_opt, full_opt, TMcut, infmt1_opt, infmt2_opt,
                        ter_opt, split_opt, outfmt_opt, fast_opt, mirror_opt,
                        het_opt, atom_opt, autojustify, mol_opt, 
                        dir_opt+dir1_opt, dir_opt+dir2_opt,
                        chain2parse1, chain2parse2, model2parse1, model2parse2,
                        tmp_vec1, tmp_vec2, byresi_opt, chainmapfile, se_opt);
                    vector<string>().swap(tmp_vec2);
                }
                vector<string>().swap(tmp_vec1);
            }
        }
        else if (dirpair_opt.size()==0) MMalign(xname, yname, fname_super,
            fname_lign, fname_matrix, sequence, d0_scale, m_opt, o_opt,
            a_opt, d_opt, full_opt, TMcut, infmt1_opt, infmt2_opt,
            ter_opt, split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
            atom_opt, autojustify, mol_opt, dir1_opt, dir2_opt,
            chain2parse1, chain2parse2, model2parse1, model2parse2,
            chain1_list, chain2_list, byresi_opt,chainmapfile, se_opt);
        else
        {
            vector<string> tmp_vec1;
            vector<string> tmp_vec2;
            for (i = 0; i < chain1_list.size(); i++)
            {
                xname = chain1_list[i];
                yname = chain2_list[i];
                tmp_vec1.push_back(xname);
                tmp_vec2.push_back(yname);
                MMalign(xname, yname, fname_super, fname_lign, fname_matrix,
                        sequence, d0_scale, m_opt, o_opt, a_opt, d_opt, full_opt,
                        TMcut, infmt1_opt, infmt2_opt, ter_opt, split_opt,
                        outfmt_opt, fast_opt, mirror_opt, het_opt, atom_opt,
                        autojustify, mol_opt, dirpair_opt, dirpair_opt,
                        chain2parse1, chain2parse2, model2parse1, model2parse2,
                        tmp_vec1, tmp_vec2, byresi_opt, chainmapfile, se_opt);
                tmp_vec1[0].clear();
                tmp_vec1.clear();
                tmp_vec2[0].clear();
                tmp_vec2.clear();
            }
        }
        chainmapfile.clear();
    }
    else if (mm_opt == 2)
        MMdock(xname, yname, fname_super,
               fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, o_opt, a_opt,
               u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
               split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
               atom_opt, autojustify, mol_opt, dir1_opt, dir2_opt,
               chain2parse1, chain2parse2, model2parse1, model2parse2,
               chain1_list, chain2_list, do_opt);
    else if (mm_opt == 3)
        ; // should be changed to mm_opt=0, cp_opt=true
    else if (mm_opt == 4)
        mTMalign(xname, yname, fname_super, fname_matrix,
                 sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt, a_opt,
                 u_opt, d_opt, full_opt, TMcut, infmt1_opt, ter_opt,
                 split_opt, outfmt_opt, fast_opt, het_opt,
                 atom_opt, autojustify, mol_opt, dir_opt, byresi_opt, chain1_list,
                 chain2parse1, model2parse1, se_opt);
    else if (mm_opt == 5 || mm_opt == 6)
        SOIalign(xname, yname, fname_super, fname_lign,
                 fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt,
                 a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
                 split_opt, outfmt_opt, fast_opt, cp_opt, mirror_opt, het_opt,
                 atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt,
                 dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2,
                 chain1_list, chain2_list, se_opt, closeK_opt, mm_opt);
    else if (mm_opt == 7)
        flexalign(xname, yname, fname_super, fname_lign,
                  fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt,
                  a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
                  split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
                  atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt,
                  dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2,
                  byresi_opt, chain1_list, chain2_list, hinge_opt, 0);
    else if (mm_opt == 8)
        flexalign(xname, yname, fname_super, fname_lign,
                  fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt,
                  a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
                  split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
                  atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt,
                  dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2,
                  byresi_opt, chain1_list, chain2_list, hinge_opt, 1);
    else if (mm_opt == 9)
        flexalign_best(xname, yname, fname_super, fname_lign,
                       fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt,
                       a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
                       split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
                       atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt,
                       dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2,
                       byresi_opt, chain1_list, chain2_list, hinge_opt);
    else if (mm_opt == 10)
        flexalign_fatcat(xname, yname, fname_super, fname_lign,
                         fname_matrix, sequence, Lnorm_ass, d0_scale, m_opt, i_opt, o_opt,
                         a_opt, u_opt, d_opt, TMcut, infmt1_opt, infmt2_opt, ter_opt,
                         split_opt, outfmt_opt, fast_opt, mirror_opt, het_opt,
                         atom_opt, autojustify, mol_opt, dir_opt, dirpair_opt, dir1_opt,
                         dir2_opt, chain2parse1, chain2parse2, model2parse1, model2parse2,
                         byresi_opt, chain1_list, chain2_list, hinge_opt);
    else
        cerr << "WARNING! -mm " << mm_opt << " not implemented" << endl;

    /* clean up */
    vector<string>().swap(chain1_list);
    vector<string>().swap(chain2_list);
    vector<string>().swap(chain2parse1);
    vector<string>().swap(chain2parse2);
    vector<string>().swap(model2parse1);
    vector<string>().swap(model2parse2);
    vector<string>().swap(sequence);
    vector<pair<string, string>>().swap(chain_pair_list);

    t2 = clock();
    float diff = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;
    if (outfmt_opt < 2)
        printf("#Total CPU time is %5.2f seconds\n", diff);
    return 0;
}
