/*
 * DisCovER.cpp
 *
 * Created  on:
 * 	Author: Sutanu Bhattacharya
 *	@Version: 1.0 (Nov 7, 2020)
 *
 * Purpose:	Distance- and orientation-based protein threading
 *
 * References:	(1) S. Wu, Y. Zhang. MUSTER: Improving protein sequence profile-profile alignments by using multiple sources of structure information. Proteins: Structure, Function, and Bioinformatics 2008
 *		(2) S Ovchinnikov et al. Protein Structure Determination using Metagenome sequence data. (2017) Science. 355(6322):294–8.
 *		(3) D Bhattacharya et al. UniCon3D: de novo protein structure prediction using united-residue conformational search via stepwise, probabilistic sampling. (2016) Bioinformatics. 
 *		(4) D Buchan, D Jones. EigenTHREADER: analogous protein fold recognition by efficient contact map threading. (2017) Bioinformatics. 
 *		(5) Y. Zhang, http://zhanglab.ccmb.med.umich.edu/NW-align
 *
 *
 * Inputs:
 * 	1. seq.fasta
 * 	2. seq.spd33 (output of SPIDER3) for query SS, SA, psi, phi 
 * 	3. library for templates SS, SA, psi, phi, PDB, MTX, DEP and template list
 * 	4. seq.prf (sequence profile) for query protein
 * 	5. seq.mtx (query mtx profile)
 * 	6. query distance map (query.rr), orientation angles (omega.txt, theta.txt, phi.txt)
 *
 * Outputs:
 *	1. top1.fasta : query-template (TOP1) alignment in MODELLER format
 *	
 */


#include <unistd.h>
#include <string>
#include <ctime>
#include <omp.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <limits.h>
#include <cstring>
#include <cmath>
#include <map>
#include <iterator>
#include <set>
#include <algorithm>
#include <functional>
#include<bits/stdc++.h> //for upper and system call
#include <ctype.h>


#define VERSION "V1.0"
using namespace std;


struct OPTS {
	std::string libdir; /* Template library directory */
	std::string list; /* list of template IDs (file) */
	std::string target; /*query*/
	std::string tardir; /* query directory */
	std::string qdist; /*directory containing predicted distance of the query */
	std::string modeller;/*path to MODELLER*/
	int nTop;/*number of TOP templates selceted by pure threadng */
	float idcut; /* sequence identity cutoff to remove homologs */
};

struct ALIGN1_RESULT {
	std::string alignment;//query-template alignment
	int seq_number;
	std::string template_name;
	double score_align;
	int idt;
	double ialg;
	double lena;
	double iii;
};


// for CMO
typedef vector<vector<int> > mtx_int;
typedef vector<vector<double> > mtx_double;
typedef vector<int> vec_int;
typedef vector<double> vec_double;
typedef vector<char> vec_char;
typedef vector<bool> vec_bool;
typedef vector<string> vec_string;

struct ALIGN2_CMO {
	std::string templ; //template
	std::string alignmentQ; //query-template alignment getting from CMO
	std::string alignmentT;
	int aln_len; //alignment length
	float score; //contact score+gap score+(profile score)
	float zscore_thread; //zscore from pure threading
};


bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);
void PrintCap(const OPTS &opts);
string read_fasta(string file);
int get8to3ss(string ssParm);
void read_query_ss_sa_angles(string file, string x, int mysec[], float mysa[], float mypsi[], float myphi[]);
float getSolAcc(char aa, string sa);
void read_mtx(string file, double mtx[][21]);
void read_prf(string file, double mtx[][21]);
void read_template_ss_sa_angles(string dssp, int ss[], float sa[], float psi[], float phi[]);
ALIGN1_RESULT align1(int seq_num, string template_name, string x, string y, double mymtx[][21], double umtxb[][21], int mysec[],
                float mysa[], float mypsi[], float myphi[], int iysec[], float iysa[], float iypsi[], float iyphi[],
                double myprf[][21], double prflib[][21], float p[]);
int getResNum(char a);
void read_DEP_profile(string dep_prf, double dep[][21], string seq);
void Blosum62Matrix(int imut[][24]);
double NeedlemanWunsch(string f1, string f2, int gap_open,int gap_extn);
void smooth_ss(int seq_length,int ss[]);
string trim(const string& str);

//for CMO
double exp_fast(double x){
	x = 1 + x/1024;
    	x *= x; x *= x; x *= x; x *= x;
    	x *= x; x *= x; x *= x; x *= x;
    	x *= x; x *= x;
    	return x;
}

double gaussian(double mean, double stdev, double x){return exp_fast(-pow((x - mean),2)/(2*(pow(stdev,2))));}

vec_int align(vec_double &gap_a, vec_double &gap_b, double &gap_e, mtx_double &sco_mtx, mtx_double &p_sco_mtx);

double Falign(double *sco_mtx, int rows, int cols);

void mod_gap(vec_double &gap, vec_char &ss, double &gap_ss_w)
{
    for(int i = 0; i < ss.size()-1; i++){
        if((ss[i] == 'H' && ss[i+1] == 'H') || (ss[i] == 'E' && ss[i+1] == 'E')){gap[i] *= gap_ss_w;}
    }
}

double sepw(double sep){if(sep <= 4){return 0.50;}else if(sep == 5){return 0.75;}else{return log10(1+sep);}}


struct point3d {
    double x;
    double y;
    double z;
};

struct pdbInfo {
    int res_num;
    point3d ca;
    point3d cb;
    point3d n;
};

struct dist_bin_5_to_20{
	mtx_double m5,m5_5,m6,m6_5,m7,m7_5,m8,m8_5,m9,m9_5,m10,m10_5,m11,m11_5,m12,m12_5,m13,m13_5,m14,m14_5;
	mtx_double m15,m15_5,m16,m16_5,m17,m17_5,m18,m18_5,m19,m19_5,m20,m20_5;
};

struct dist_bin_6_to_14_2{
        mtx_double m6,m8,m10,m12,m14;
};

struct dist_bin_5_to_15_1{
        mtx_double m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15;
};

struct dihedral_angles{
	mtx_double a_15,a_30,a_45,a_60,a_75,a_90,a_105,a_120,a_135,a_150,a_165,a_180,a_195,a_210,a_225,a_240,a_255,a_270,a_285,a_300,a_315,a_330,a_345,a_360;
};
struct planar_angles{
	mtx_double a_15,a_30,a_45,a_60,a_75,a_90,a_105,a_120,a_135,a_150,a_165,a_180;
};

vec_int load_data (string file, int sep_cutoff, mtx_double &mtx, vec_int &vec_div, vec_int &vec, mtx_int &vec_i, mtx_double &prf, int size, vec_char &ss,dist_bin_5_to_20 &rr, double n_bins);
void ini_prf_SCO(mtx_double &P_SCO, double &prf_w, mtx_double &prf_a, vec_char &aa_a, mtx_double &prf_b, vec_char &aa_b);
void ini_SCO(double sep_x, double sep_y, mtx_double &SCO,vec_int &vec_a_div,vec_int &vec_b_div,vec_int &vec_a,vec_int &vec_b,mtx_int &vec_a_i,mtx_int &vec_b_i,mtx_double &mtx_a,mtx_double &mtx_b,dist_bin_5_to_20 query_mtx,dist_bin_5_to_20 temp_mtx,dihedral_angles omega,dihedral_angles theta,planar_angles phi,map<int,pdbInfo> &all_atomsT);
vec_int mod_SCO(double do_it, vec_double &gap_a, vec_double &gap_b, double gap_e, mtx_double &SCO, mtx_double &P_SCO,vec_int &vec_a_div,vec_int &vec_b_div,vec_int &vec_a, vec_int &vec_b,mtx_int &vec_a_i, mtx_int &vec_b_i,mtx_double &mtx_a, mtx_double &mtx_b,dist_bin_5_to_20 query_mtx,dist_bin_5_to_20 temp_mtx,dihedral_angles omega,dihedral_angles theta,planar_angles phi,map<int,pdbInfo> &all_atomsT);
void chk (vec_double &gap_a, vec_double &gap_b, double &gap_e_w, double& con_sco,double& gap_sco,double& prf_sco,vec_int& vec_a_div,mtx_int& vec_a_i,mtx_double& mtx_a,mtx_double& mtx_b,vec_int& a2b,mtx_double &P_SCO,dist_bin_5_to_20 query_mtx,dist_bin_5_to_20 temp_mtx,dihedral_angles omega,dihedral_angles theta,planar_angles phi,map<int,pdbInfo> &all_atomsT);

ALIGN2_CMO cmo(string templ, string qdist, string libdir, float zscore_thread, string target, string seqA, int lenQ, vec_int m2n_a, mtx_double mtx_a, vec_int vec_a,vec_int vec_a_div,mtx_int vec_a_i,mtx_double prf_a, vec_char ss_a, dist_bin_5_to_20 query_mtx,dihedral_angles omega_all_bins, dihedral_angles theta_all_bins, planar_angles phi_all_bins);

void ini_prf_SCO(mtx_double &P_SCO, double &prf_w, mtx_double &prf_a, vec_char &aa_a, mtx_double &prf_b, vec_char &aa_b);

//MODELLER
void same_ali(string model);
void generateModeller_SS_CM_Script(float t);
void generateBasicModellerScript(string temp);
void run_modeller(string top1_temp,string modeller,string tempDIR);


/*inter-residue orientation */
const double PI = std::atan(1.0)*4;
const double RAD2DEG = 180/PI;
const double DEG2RAD = PI/180;

map<int,pdbInfo> pdb2coord(string file);
point3d getDifference(point3d &p1, point3d &p2);
double getDistance(point3d &p1, point3d &p2);
double getDotProduct(point3d &p1, point3d &p2);
point3d getCrossProduct(point3d &p1, point3d &p2);
double getAngle(point3d &p1, point3d &p2, point3d &p3);
double getDihedral(point3d &p1, point3d &p2, point3d &p3, point3d &p4);
void read_orientation_query(string file,dihedral_angles &a,int lenQ); 
void read_orientation_phi_query(string file,planar_angles &phi,int lenQ);

/** main() */
int main(int argc, char *argv[]) {
	OPTS opts = { "", "", "", "", "", "", 50, 0.30 }; //default seq id cutoff 0.30
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts(opts);
		return 1;
	}
	PrintCap(opts);

	int seq_num = 1; //track number of templates processed

	float idcut = opts.idcut; // sequence identity cutoff 
	int nTop = opts.nTop; //number of TOP templates selected by pure threadng

	string libdir = opts.libdir;
	string MTX = libdir + "MTX/";
	string SEQ = libdir + "seq/";
	string DEP = libdir + "DEP/";
	string dssp = libdir + "dssp/";
	string PDB = libdir + "PDB/";
	string qdist = opts.qdist;

	// read the query sequence
	string x = "";
	x = read_fasta("seq.fasta");
	int x_len = x.length();

	// Read query SS, SA, psi, phi from .spd33 file
	int mysec [x_len];
	std::fill_n(mysec, x_len, 3); // default 3
	float mysa [x_len];
	std::fill_n(mysa, x_len, 0); 
	float mypsi [x_len];
	std::fill_n(mypsi, x_len, 0); 	
	float myphi [x_len];
	std::fill_n(myphi, x_len, 0); 
	read_query_ss_sa_angles("seq.spd33",x,mysec, mysa, mypsi, myphi);
	smooth_ss(x_len,mysec);

	//read query MTX file	
	double myprf_mtx [x_len][21];
	memset(myprf_mtx, 0, x_len * 21 * sizeof(double));  
	read_mtx("seq.mtx", myprf_mtx);
		
	//read query frequency profile
	double myprf [x_len][21];
	memset(myprf, 0, x_len * 21 * sizeof(double));  
	read_prf("seq.prf", myprf);


	/* read template IDs */
	std::vector<std::string> listB;
	ifstream in(opts.list.c_str()); 
	if(!in) {
		printf("Error: cannot open template list file '%s'\n", opts.list.c_str());
		return 1;
	}
	string name = "";
	while (std::getline(in, name)) {
		if(name.size() > 0) {
			listB.push_back(name);
		}
	}
	in.close();

	//weights 0,g0,ge,ss,sa,shift,hydro,phi,0,1,depthLib,psi
	float weight[12] = {0, 7.01, 0.55, 0.66, 1.6, -0.99, 0.31, 0.19, 0, 1.0, 0.39, 0.19}; //weights, p[0] is place holder

	//read PDBs one by one and calculate alignments 
	ofstream fw; // open a file in write mode.
        fw.open("align.dat");

	// calculate alignment score for each PDB 
	ofstream fw_rst;
	fw_rst.open("rst.dat");


	map<string,ALIGN1_RESULT> template_align1;//store return values of align1 method for each template

	for (unsigned i = 0; i < listB.size(); i++) { // for schedule(dynamic) will access PDBs one by one as opposed to parallel
		string templatePDB = listB[i]; // make sure whether all libraries has this template
		string dep_file = DEP + templatePDB + ".dep";
		string dssp_fil = dssp + templatePDB + ".dssp";
		string seq_fil = SEQ + templatePDB + ".seq";
		string mtx_file = MTX + templatePDB + ".mtx";
		
		string y = read_fasta(seq_fil);
		int seq_length = y.length();

		//exclude templates with seq id > id cutoff to the query protein
		string f1 = x;
		string f2 = y;
		transform(f1.begin(), f1.end(), f1.begin(), ::toupper);
		transform(f2.begin(), f2.end(), f2.begin(), ::toupper);

		if (idcut<1.0){
			double identity = NeedlemanWunsch(f2,f1,-11,-1); //gap_open=-11,gap_extn=-1;
			if (identity > idcut){continue;} 
		}

		//read template MTX file
		double umtxb [seq_length][21];
		memset(umtxb, 0, seq_length * 21 * sizeof(double));  
		read_mtx(mtx_file, umtxb);

		//read template DEP file
		double dep_prof [seq_length][21];
		memset(dep_prof, 0, seq_length * 21 * sizeof(double));
		read_DEP_profile(dep_file, dep_prof, y);

		//Read template SS, SA, psi, phi from .dssp file
		int libsec [seq_length];
		std::fill_n(libsec, seq_length, 3);
		float libsa [seq_length];
		std::fill_n(libsa, seq_length, 0);
		float libpsi [seq_length];
		std::fill_n(libpsi, seq_length, 0);
		float libphi [seq_length];
		std::fill_n(libphi, seq_length, 0);

		read_template_ss_sa_angles(dssp_fil, libsec, libsa, libpsi, libphi);
		smooth_ss(seq_length,libsec); //smooth SS according to EigenTHREADER

		//calculate alignment score and query-template alignment
		ALIGN1_RESULT result;
		result.template_name = "", result.alignment = "";
		result.seq_number = 0, result.score_align = 0, result.idt = 0, result.ialg = 0, result.lena = 0, result.iii = 0;

		result = align1(seq_num, templatePDB, x, y, myprf_mtx, umtxb, mysec, mysa, mypsi, myphi, libsec, libsa, libpsi, libphi, myprf, dep_prof, weight);
		template_align1[result.template_name] = result;

		seq_num += 1;
	}	

	//calculate z-score for each template 
	//instead of seq_num-1 because 0th index will be ignored
		
	double v1 [seq_num];
	std::fill_n(v1, seq_num, 0); 
	double v2 [seq_num];
	std::fill_n(v2, seq_num, 0);   
	double sequence_id [seq_num];
	std::fill_n(sequence_id, seq_num, 0);  
	double z1 [seq_num];
	std::fill_n(z1, seq_num, 0);  
	double z2 [seq_num];
	std::fill_n(z2, seq_num, 0); 


	string ppa_alignments = "", tempName = "";
	string p [seq_num];
	map<string,ALIGN1_RESULT>::iterator itr; //iterate map
	int i = 0, k;
	for (itr = template_align1.begin(); itr != template_align1.end(); ++itr) {
		tempName = itr->first;
		ALIGN1_RESULT align1_out = itr->second;
		i += 1;
		p[i] = tempName; //template
		v1[i] = align1_out.score_align; //raw score
		v2[i] = 0.000; //reverse raw score
		sequence_id[i] = align1_out.idt; //seq id
		z1[i] = v1[i] / align1_out.lena; //full alignment
		z2[i] = v1[i] / align1_out.iii; //partial alignment

		//save alignment and score details in align.dat and rst.dat files respectively
		ppa_alignments = align1_out.alignment;
		size_t found = ppa_alignments.find(" "); //find space between structure and query alignments
		fw << ">P1;" + tempName + "\n";
		fw << "structureX:" + tempName + "  :1 :  :   :  : t1 : t2 : t3 :\n";
		for (k = 0; k < found; k++){
			fw << ppa_alignments[k];
		}
		fw << "*\n";
		fw << ">P1;query\nsequence:query:1 :  :   :  : q1 : q2 : q3 :\n";
		for (k = found+1; k < ppa_alignments.length(); k++){ //ignoring the space
			fw << ppa_alignments[k];
		}
		fw << "*\n\n";

		fw_rst << std::to_string(i) + "\t" + tempName + "\t" + std::to_string(v1[i]) + "\t" + std::to_string(v2[i]) + "\t" + 
				std::to_string(sequence_id[i]) + "\t" + std::to_string(align1_out.lena) + "\t" + 
				std::to_string(align1_out.iii) + "\n"; 
	}
	fw.close();
	fw_rst.close();

	
	int N_hit = i;
	float z1_a = 0.0, z1_a2 = 0.0, z2_a = 0.0, z2_a2 = 0.0;
	for(i = 1; i<N_hit+1; i++){
		z1_a += z1[i];
		z1_a2 += pow(z1[i], 2);
		z2_a += z2[i];
		z2_a2 += pow(z2[i], 2);
	}

	z1_a /= float(N_hit);
	z1_a2 /= float(N_hit);
	float z1_sd = sqrt(z1_a2 - pow(z1_a, 2));

	map<string,float> z1_zscore;
	map<string,int> TT1;
	for(i = 1; i<N_hit+1; i++){
		z1_zscore[p[i]] = float (z1_a - z1[i]) / z1_sd;
		TT1[p[i]] = i;
	}
	
	//sort map in descending order
	// Declaring the type of Predicate that accepts 2 pairs and return a bool
	typedef std::function<bool(std::pair<std::string, float>, std::pair<std::string, float>)> Comparator;
	// Defining a lambda function to compare two pairs. It will compare two pairs using second field
	Comparator compFunctor =
	 	[](std::pair<std::string, float> elem1 ,std::pair<std::string, float> elem2)
		{
			return elem1.second > elem2.second;
		};		

	//sort z1_zscore in descending order
	std::set<std::pair<std::string, float>, Comparator> z1_zscore_keys(z1_zscore.begin(), z1_zscore.end(), compFunctor);

	//sort z2_zscore in descending order
	z2_a /= float(N_hit);
        z2_a2 /= float(N_hit);
        float z2_sd = sqrt(z2_a2 - pow(z2_a, 2));

        map<string,float> z2_zscore;
        map<string,int> TT2;
        for(i = 1; i<N_hit+1; i++){
                z2_zscore[p[i]] = float (z2_a - z2[i]) / z2_sd;
                TT2[p[i]] = i;
        }
	
	std::set<std::pair<std::string, float>, Comparator> z2_zscore_keys(z2_zscore.begin(), z2_zscore.end(), compFunctor);

	string index1, index2;	
	for (std::pair<std::string, float> element : z1_zscore_keys){
		index1 = element.first;
		break;
	}

	for (std::pair<std::string, float> element : z2_zscore_keys){
                index2 = element.first;
                break;
        }

	int score_flag=1; //always first scheme unless the following:
	if( (sequence_id[TT1[index1]]) + 0.01 <= sequence_id[TT2[index2]] ){
		score_flag=2;
	} else if ( (sequence_id[TT2[index2]]) + 0.01 <= sequence_id[TT1[index1]] ) {
		score_flag=1;
	}

	//cout<<"check...\n";
	//read query rr and orientation angles
	int lowest_bin=5, highest_bin=20;
	double bin_length=0.5, n_bins=31;
	mtx_double mtx_a; vec_int vec_a; vec_int vec_a_div; mtx_int vec_a_i; mtx_double prf_a; vec_char ss_a;	
	dist_bin_5_to_20 query_mtx;	
	//dist_bin_6_to_14_2 query_mtx;
	//dist_bin_5_to_15_1 query_mtx;
	int sep_cutoff = 5;
	//cout<<"read query.rr ";
	vec_int m2n_a = load_data("query.rr",sep_cutoff,mtx_a,vec_a_div,vec_a,vec_a_i,prf_a,x_len,ss_a,query_mtx,n_bins);
	//cout << "done..\nread omega ";
	dihedral_angles om;
	read_orientation_query("omega.txt",om,x_len);
	//cout << "done..\n";
	dihedral_angles th;
	read_orientation_query("theta.txt",th,x_len);
	planar_angles ph;
	read_orientation_phi_query("phi.txt",ph,x_len);
	//cout << "theta and phi are done..\n";
	//find best-fit templates
	int i_t = 0;
	float zscore_value=0.0;
	string top_template_name="";
	string top1_template="";
	map<string,ALIGN2_CMO> template_align2;

	int shrink_pool = nTop;
	for(i = 1; i<shrink_pool+1; i++){
		ALIGN2_CMO result2;
		if(score_flag == 0) {
			continue;
		} else if (score_flag == 1) {
			int index_set = 0;
			for (std::pair<std::string, float> element : z1_zscore_keys){
				top_template_name = element.first;

				if(i==1){top1_template = top_template_name;}

				zscore_value = element.second;
				if (index_set == i-1){
					//cout << top_template_name << "\n";
					result2 = cmo(top_template_name,qdist,libdir,zscore_value,opts.target,x ,x_len,m2n_a,mtx_a,vec_a,vec_a_div,vec_a_i,prf_a,ss_a,query_mtx,om,th,ph);
					break;
				}
				index_set += 1;
			}
			
		} else {
			int index_set = 0;
			for (std::pair<std::string, float> element : z2_zscore_keys){
				top_template_name = element.first;

				if(i==1){top1_template = top_template_name;}

				zscore_value = element.second;
				if (index_set == i-1){
					//cout << top_template_name << "\n";
					result2 = cmo(top_template_name,qdist,libdir,zscore_value,opts.target,x ,x_len,m2n_a,mtx_a,vec_a,vec_a_div,vec_a_i,prf_a,ss_a,query_mtx,om,th,ph);
                                        break;
                                }	
				index_set += 1;
			}
		}

		template_align2[result2.templ] = result2;

	}

	template_align1.clear(); //clear map
	z1_zscore.clear();
	z2_zscore.clear();
	TT1.clear();
	TT2.clear();

	// rank templates based on z_contact_score
	tempName = "";
	
	//calculate z-score for each template using contact score and alignment length
	double v_cmo[shrink_pool+1];
	std::fill_n(v_cmo, shrink_pool+1, 0);
	double z_cmo[shrink_pool+1]; 	
	std::fill_n(z_cmo, shrink_pool+1, 0);	
	double z_thread_sc[shrink_pool+1];
	std::fill_n(z_thread_sc, shrink_pool+1, 0);

	map<string,ALIGN2_CMO>::iterator itr2; //iterate map
	ALIGN2_CMO align2_out;
	int i_cmo = 0;
	string p_cmo[shrink_pool+1];
	for (itr2 = template_align2.begin(); itr2 != template_align2.end(); ++itr2) {
		tempName = itr2->first;
		align2_out = itr2->second;
		if(align2_out.aln_len <= 0){continue;}
		i_cmo += 1;
		p_cmo[i_cmo] = tempName;
		v_cmo[i_cmo] = align2_out.score; //contact score
		z_cmo[i_cmo] = v_cmo[i_cmo] / float(align2_out.aln_len); 
		z_thread_sc[i_cmo] = align2_out.zscore_thread; //z_score from pure threading
	}

	float z_a_cmo = 0.0, z_a2_cmo = 0.0;
	for(i_cmo = 1; i_cmo<shrink_pool+1; i_cmo++){
		z_a_cmo += z_cmo[i_cmo];
		z_a2_cmo += pow(z_cmo[i_cmo], 2);
	}
	
	z_a_cmo /= float(shrink_pool);
	z_a2_cmo /= float(shrink_pool);
	float z_sd_cmo = sqrt(z_a2_cmo - pow(z_a_cmo, 2));

	map<string,float> z_zscore_cmo;
	for(i_cmo = 1; i_cmo<shrink_pool+1; i_cmo++){
		z_zscore_cmo[p_cmo[i_cmo]] = float(z_cmo[i_cmo] - z_a_cmo) / z_sd_cmo;
	}

	//sort map in descending order
	string top1_templ_CMO = "";


	//combine z_score = w1*z_score + w2*z_contact_score
	map<string,float> combine_zscore;
	for(i_cmo = 1; i_cmo<shrink_pool+1; i_cmo++){
		combine_zscore[p_cmo[i_cmo]] = float(1 * z_thread_sc[i_cmo] + 1.0 * z_zscore_cmo[p_cmo[i_cmo]]);
	}

	//sort map in descending order
	std::set<std::pair<std::string, float>, Comparator> sorted_z_zscore_cmo(combine_zscore.begin(), combine_zscore.end(), compFunctor);

	//print score for each template
	cout << "#Template\tZ_final\n";
	for (std::pair<std::string, float> element : sorted_z_zscore_cmo){
		top1_templ_CMO = element.first;
		float z_cmo = element.second;
                align2_out = template_align2[top1_templ_CMO];
		//template\t zscore_thread\t contact score\t alignemnt length\t z_score_cmo
		cout << top1_templ_CMO << "\t\t" << z_cmo << "\n";
	}
	
	// extract query-template alignment for the TOP1 template
	for (std::pair<std::string, float> element : sorted_z_zscore_cmo){
		top1_templ_CMO = element.first;
		align2_out = template_align2[top1_templ_CMO];

		// extract alignment between query and TOP1 template 
		ofstream fw_aln;
		fw_aln.open("top1.fasta");

		ofstream fw_aln2;
		fw_aln2.open("top1_aln.fasta"); 

		string alignmentQ = align2_out.alignmentQ;
		string alignmentT = align2_out.alignmentT;
		
		fw_aln << ">P1;query\nsequence:query:1 :  :   :  : q1 : q2 : q3 :\n";
		fw_aln2 << ">"+opts.target+"\n";
		for(k=0; k<alignmentQ.length(); k++){
			fw_aln << alignmentQ[k];
			fw_aln2 << alignmentQ[k];
		}				
		fw_aln << "*\n";
		fw_aln2 << "\n";
		fw_aln << ">P1;" + top1_templ_CMO + "\n";
		fw_aln2 << ">" + top1_templ_CMO + "\n"; 
		fw_aln << "structureX:" + top1_templ_CMO + "  :1 :  :   :  : t1 : t2 : t3 :\n";
		for(k=0; k<alignmentT.length(); k++){
                        fw_aln << alignmentT[k];
			fw_aln2 << alignmentT[k];  
                }
		fw_aln << "*\n";
		fw_aln2.close();
		fw_aln.close();
		
		break;
	}

	template_align2.clear();
	z_zscore_cmo.clear();
	combine_zscore.clear();
	
	cout << "\n#DisCovER is done. First-ranked query-template alignment is available in 'top1.fasta' file.\n";

	/*run MODELLER if MODELLER path is given*/
	string modeller = opts.modeller;
	if(modeller != ""){
		//cout << "#Running Modeller ...\n";
		//run_modeller(top1_templ_CMO, modeller, PDB);
	}
	//system("rm con*.rr");	
	//system("rm *.dat");
	system("mv top1_aln.fasta top1.fasta");
	
	//finish date/time
	time_t timer;
	time(&timer);
	struct tm* tm_info = localtime(&timer);
	char buf[100];
	strftime(buf, 26, "%Y:%m:%d / %H:%M:%S", tm_info);
	printf("# %s\n", std::string(80, '-').c_str());
	printf("# %20s : %s\n", "end date/time", buf);
	printf("# %s\n", std::string(80, '-').c_str());

	return 0;	
}


/** check whether a file exists */
bool fexists(const std::string& filename) {
  	std::ifstream ifile(filename.c_str());
  	return (bool)ifile;
}

/* Smooth DSSP secondary structure definitions */
void smooth_ss(int seq_length,int ss[]){
	for(int i=2; i<seq_length-1; i++){
		if(ss[i-1] == ss[i+1]){ss[i] = ss[i-1];}
	}
}


/** extract contacts from PDB */
struct AtomRecord
{

	char chainId; /* Chain identifier */
	int resNum; /* Residue sequence number */
	double x, y, z; /* Orthogonal coordinates for X, Y, Z */
};

double eucl_dist(AtomRecord A, AtomRecord B){
        double x = A.x - B.x;
        double y = A.y - B.y;
        double z = A.z - B.z;
        return sqrt(x * x + y * y + z * z);
}
void pdb2con(string file){
	string line;
        ifstream in(file);
	//ofstream fw_5,fw_6,fw_7,fw_8,fw_9,fw_10,fw_11,fw_12,fw_13,fw_14,fw_15;
	AtomRecord atom_rec;
	map<int,AtomRecord> all_atoms;
	ofstream fw;
	fw.open("template.rr");
	//fw_5.open("con5T.rr");fw_7.open("con7T.rr");fw_9.open("con9T.rr");fw_11.open("con11T.rr");fw_13.open("con13T.rr");fw_15.open("con15T.rr");
	//fw_6.open("con6T.rr");fw_8.open("con8T.rr");fw_10.open("con10T.rr");fw_12.open("con12T.rr");fw_14.open("con14T.rr");
        while(getline(in,line)){
		if (trim(line.substr(0,5).c_str()) == "ATOM" && trim(line.substr(12,4).c_str()) == "CA"){
                        atom_rec.chainId = line[21];
                        atom_rec.resNum = stoi(trim(line.substr(22,4).c_str()));
                        atom_rec.x = stof(trim(line.substr(30,8).c_str()));
                        atom_rec.y = stof(trim(line.substr(38,8).c_str()));
                        atom_rec.z = stof(trim(line.substr(46,8).c_str()));
                }
                all_atoms[atom_rec.resNum] = atom_rec;		
	}
	in.close();
        int nRes = all_atoms.size();
	for (int i = 1; i <= nRes-5; i++){
                for (int j = i+5; j <= nRes; j++){
                        if(all_atoms[i].chainId == all_atoms[j].chainId){
                                double dist = eucl_dist(all_atoms[i],all_atoms[j]);
				if (dist <= 5.0){fw << i << " " << j << " 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;}	
				if (dist <= 5.5){fw << i << " " << j << " 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;}
				if (dist <= 6.0){fw << i << " " << j << " 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 6.5){fw << i << " " << j << " 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 7.0){fw << i << " " << j << " 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 7.5){fw << i << " " << j << " 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 8.0){fw << i << " " << j << " 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 8.5){fw << i << " " << j << " 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 9.0){fw << i << " " << j << " 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 9.5){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 10.0){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 10.5){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 11.0){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 11.5){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 12.0){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 12.5){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 13.0){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 13.5){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 14.0){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 14.5){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 15.0){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 15.5){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 16.0){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 16.5){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 17.0){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1\n";continue;} 
				if (dist <= 17.5){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1\n";continue;} 
				if (dist <= 18.0){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1\n";continue;} 
				if (dist <= 18.5){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1\n";continue;} 
				if (dist <= 19.0){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1\n";continue;} 
				if (dist <= 19.5){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1\n";continue;} 
				if (dist <= 20.0){fw << i << " " << j << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1\n";continue;} 
				/*
                                if (dist <= 6.0){fw_6 << i << " " << j << " 0 6 " << dist << "\n";}
				if (dist <= 7.0){fw_7 << i << " " << j << " 0 7 " << dist << "\n";}
				if (dist <= 8.0){fw_8 << i << " " << j << " 0 8 " << dist << "\n";}
				if (dist <= 9.0){fw_9 << i << " " << j << " 0 9 " << dist << "\n";}
				if (dist <= 10.0){fw_10 << i << " " << j << " 0 10 " << dist << "\n";}
				if (dist <= 11.0){fw_11 << i << " " << j << " 0 11 " << dist << "\n";}
				if (dist <= 12.0){fw_12 << i << " " << j << " 0 12 " << dist << "\n";}
				if (dist <= 13.0){fw_13 << i << " " << j << " 0 13 " << dist << "\n";}
				if (dist <= 14.0){fw_14 << i << " " << j << " 0 14 " << dist << "\n";}
				if (dist <= 15.0){fw_15 << i << " " << j << " 0 15 " << dist << "\n";}
				*/
                        }
                }
        }
	/*
	fw_6.close();
        fw_8.close();
        fw_10.close();
        fw_12.close();
        fw_14.close();
	*/
	fw.close();
	all_atoms.clear();
}

//Contact Map Overlap
ALIGN2_CMO cmo(string templ, string qdist, string libdir, float zscore_thread, string target, string x, int lenQ, vec_int m2n_a, mtx_double mtx_a, vec_int vec_a,vec_int vec_a_div,mtx_int vec_a_i,mtx_double prf_a, vec_char ss_a, dist_bin_5_to_20 query_mtx,dihedral_angles omega_all_bins, dihedral_angles theta_all_bins, planar_angles phi_all_bins){

	string SEQ = libdir + "seq/";
        string DEP = libdir + "DEP/";
        string dssp = libdir + "dssp/";
        string PDB = libdir + "PDB/";

	string distT = PDB + templ + ".pdb";
	pdb2con(distT);//create con6T.rr,con8T.rr,con10T.rr,con12T.rr,con14T.rr 

	string seqB = SEQ + templ + ".seq"; //template sequence
	string y = read_fasta(seqB);
        int lenT = y.length(); //eliminating leading *

	map<int,pdbInfo> all_atomsT;
        all_atomsT = pdb2coord(distT);

	//query profile
	double prfA [lenQ][21];
	memset(prfA, 0, lenQ * 21 * sizeof(double));
	read_prf("seq.prf", prfA); 	

	//template profile
	string dep_file = DEP + templ + ".dep";
	double prfB [lenT][21];
	memset(prfB, 0, lenT * 21 * sizeof(double));	
	read_DEP_profile(dep_file, prfB, y);

	double gap_open = -1;
	double gap_ext = -0.01;
	int sep_cutoff = 6;
	int iter = 20;
	bool use_gap_ss = true;
	double gap_ss_w = 2;
	bool use_prf = true;
	double prf_w = 1;


	double gap_ext_w = fabs(gap_ext)/fabs(gap_open);

	int size_a = mtx_a.size();
	prf_a.resize(lenQ,vector<double>(20,0));
	for(int i=1;i<lenQ;i++){
		for (int k = 1; k <= 20; k++) {
			prf_a[i-1][k-1] = prfA[i][k];
		}
	}
	ss_a.resize(lenQ,'X');
	// Read query SS, SA, psi, phi from .spd33 file
	int mysec [lenQ]; std::fill_n(mysec, lenQ, 3); // default 3
	float mysa [lenQ]; std::fill_n(mysa, lenQ, 0);
	float mypsi [lenQ]; std::fill_n(mypsi, lenQ, 0);
	float myphi [lenQ]; std::fill_n(myphi, lenQ, 0);
	read_query_ss_sa_angles("seq.spd33",x,mysec, mysa, mypsi, myphi);
	smooth_ss(lenQ,mysec);
	for(int i=1;i<lenQ;i++){
		if (mysec[i] == 1){ss_a[i-1] = 'H';}
		else if(mysec[i] == 2){ss_a[i-1] = 'E';}
		else{ss_a[i-1] = 'C';}
	}

	// if use_gap_ss on, modify the gap penalities
	vec_double gap_a(size_a,gap_open);if(use_gap_ss == true){mod_gap(gap_a,ss_a,gap_ss_w);}

	
	// load data from contact map B
	mtx_double mtx_b; vec_int vec_b; vec_int vec_b_div; mtx_int vec_b_i; mtx_double prf_b; vec_char ss_b;
	dist_bin_5_to_20 template_mtx;
	//cout << "\ntemplate.rr";
	vec_int m2n_b = load_data("template.rr",sep_cutoff,mtx_b,vec_b_div,vec_b,vec_b_i,prf_b,lenT,ss_b,template_mtx,31);

        int size_b = mtx_b.size();
	prf_b.resize(lenT,vector<double>(20,0));
	for(int i=1;i<lenT;i++){
                for (int k = 1; k <= 20; k++) {
                        prf_b[i-1][k-1] = prfB[i][k];
                }
        }
	ss_b.resize(lenT,'X');
	string dssp_fil = dssp + templ + ".dssp";
	int libsec [lenT]; std::fill_n(libsec, lenT, 3);
	float libsa [lenT]; std::fill_n(libsa, lenT, 0);
	float libpsi [lenT]; std::fill_n(libpsi, lenT, 0);
	float libphi [lenT]; std::fill_n(libphi, lenT, 0);
	read_template_ss_sa_angles(dssp_fil, libsec, libsa, libpsi, libphi);
	smooth_ss(lenT,libsec);
	for(int i=1;i<lenT;i++){
		if (libsec[i] == 1){ss_b[i-1] = 'H';}
                else if(libsec[i] == 2){ss_b[i-1] = 'E';}
                else{ss_b[i-1] = 'C';}
	}

	// if use_gap_ss on, modify the gap penalities
	vec_double gap_b(size_b,gap_open);if(use_gap_ss == true){mod_gap(gap_b,ss_b,gap_ss_w);}

	// if use_prf on, initialize profile SCO matrix
	vec_char aa_a;vec_char aa_b;
	aa_a.resize(lenQ,'X');
	for(int i=1;i<lenQ;i++){
		aa_a[i] = x[i];
	}
	aa_b.resize(lenT,'X');
        for(int i=1;i<lenT;i++){
                aa_b[i] = y[i];
        }
	mtx_double P_SCO;if(use_prf == true){ini_prf_SCO(P_SCO,prf_w,prf_a,aa_a,prf_b,aa_b);}
	//cout << "\nstart wlignmnet";

	// STARTING ALIGNMENT!!!
	// keeping track of the BEST alignment
	int max_sep_x = 0;
	int max_sep_y = 0;
	int max_g_e = 0;
	double con_max = -1;
	double gap_max = 0;
    	double prf_max = 0;
    	vec_int a2b_max;

	// try different sep (sequence seperation difference) penalities
	vec_double sep_x_steps {0,1,2}; // (constant, linear, quadratic)
	for(int sx = 0; sx < sep_x_steps.size(); sx++){
		double sep_x = sep_x_steps[sx];
		//try different scaling factors for sep penalities
		vec_double sep_y_steps {1,2,4,8,16,32};
		for(int sy = 0; sy < sep_y_steps.size(); sy++){
			double sep_y = sep_y_steps[sy];
			// Get initial score matrix
			mtx_double C_SCO(size_a,vector<double>(size_b,0));
			//cout << "\ncall INIT "<<sep_x<<" "<<sep_y;
			//cout << "\ncheck....";
			ini_SCO(sep_x,sep_y,C_SCO,vec_a_div,vec_b_div,vec_a,vec_b,vec_a_i,vec_b_i,mtx_a,mtx_b,query_mtx,template_mtx,omega_all_bins,theta_all_bins,phi_all_bins,all_atomsT);
			//cout << "\nINIT done..";


			// try different gap_ext penalities!
			vec_double gap_e_steps {5,10,100,1000};
			for(int g_e = 0; g_e < gap_e_steps.size(); g_e++){
				//cout << "\ncall INIT "<<sep_x<<" "<<sep_y<<" "<<gap_e_steps[g_e];
				double gap_e = 1/gap_e_steps[g_e];
				// restart SCO matrix
				mtx_double SCO = C_SCO;
				// get alignment (a2b mapping) after X iterations
				vec_int a2b = mod_SCO(iter,gap_a,gap_b,gap_e,SCO,P_SCO,vec_a_div,vec_b_div,vec_a,vec_b,vec_a_i,vec_b_i,mtx_a,mtx_b,query_mtx,template_mtx,omega_all_bins,theta_all_bins,phi_all_bins,all_atomsT);
				// compute number of contacts/gaps made
				double con_sco = 0;
				double gap_sco = 0;
				double prf_sco = 0;
				chk(gap_a,gap_b,gap_ext_w,con_sco,gap_sco,prf_sco,vec_a_div,vec_a_i,mtx_a,mtx_b,a2b,P_SCO,query_mtx,template_mtx,omega_all_bins,theta_all_bins,phi_all_bins,all_atomsT);


				// save if BEST!
				//cout << "if best ";
				if(con_sco+gap_sco+prf_sco > con_max+gap_max+prf_max){
					max_sep_x = sep_x;
					max_sep_y = sep_y;
					max_g_e = g_e;
					con_max = con_sco;
					gap_max = gap_sco;
					prf_max = prf_sco;
					a2b_max = a2b;
				}
				//cout << "save ";
			}
		}
	}
	//cout << "\nReport the BEST score ";
	// Report the BEST score
	//cout << "\ncheck ..";
	ALIGN2_CMO result;
	//cout << "\na2b_max";
	int aln_len = 0;for(int ai = 0; ai < size_a; ai++){int bi = a2b_max[ai];if(bi != -1){aln_len++;}}
	float score = (float) (con_max+gap_max+prf_max); //con_max+gap_max+prf_max
	//cout << "\nscore: "<<score;
	result.aln_len = aln_len;
	result.score = score;
	result.zscore_thread = zscore_thread;
	result.alignmentQ = "";
	result.alignmentT = "";
	result.templ = templ;
	//cout << "\nalignment..";
	// generate alignment
	string seqQ = x;
	string seqT = y;
	string target_seq_dash = "";
	string template_seq_dash = "";
	int j, i, k;
	int prevI = 0; int prevJ = 0;

	for(int a = 0; a < size_a; a++){	
		int b = a2b_max[a];
		if(b != -1){
			i = m2n_a[a]; j = m2n_b[b];
			if (prevI<(i-1)){
				for(int k1 = prevI+1; k1<i; k1++){
					target_seq_dash += seqQ[k1];
					template_seq_dash += "-";
				}
				prevI = i;
			}
			if (prevJ<(j-1)){
				for(int k1 = prevJ+1; k1<j; k1++){
					template_seq_dash += seqT[k1];
					target_seq_dash += "-";
				}
				prevJ = j;
			}
			template_seq_dash += seqT[j];
			target_seq_dash += seqQ[i];
			prevI = i;
			prevJ = j;						
		}
	}
	if (j<seqT.length()){
		for (k = j+1; k < seqT.length(); k++){
			template_seq_dash += seqT[k];
			target_seq_dash += "-";
		}
	}
	if (i<seqQ.length()){
		for (k = i+1; k < seqQ.length(); k++){
			template_seq_dash += "-";
			target_seq_dash += seqQ[k];
		}
	}

	result.alignmentQ = target_seq_dash;
	result.alignmentT = template_seq_dash;

	return result;
}

vec_int align(vec_double &gap_a, vec_double &gap_b, double &gap_e, mtx_double &sco_mtx, mtx_double &p_sco_mtx){
	 // LOCAL_ALIGN
	 // Start	0
	 // [A]lign	1
	 // [D]own	2
	 // [R]ight	3

	double max_sco = 0;
	int rows = sco_mtx.size();
    	int cols = sco_mtx[0].size();

	bool add_prf = true;

	vec_int a2b(rows,-1);

	mtx_double sco(rows+1,vector<double>(cols+1,0));
	mtx_int label(rows+1,vector<int>(cols+1,0));

	int max_i = 0;int max_j = 0;
	for (int i = 1; i <= rows; i++){
		for (int j = 1; j <= cols; j++){
			double A = sco[i-1][j-1] + sco_mtx[i-1][j-1]; if(add_prf == true){A += p_sco_mtx[i-1][j-1];}
			double D = sco[i-1][j];
			double R = sco[i][j-1];

			if(label[i-1][j] == 1){D += gap_b[j-1];}else{D += gap_b[j-1] * gap_e;}
			if(label[i][j-1] == 1){R += gap_a[i-1];}else{R += gap_a[i-1] * gap_e;}

			if(A <= 0 and D <= 0 and R <= 0){label[i][j] = 0;sco[i][j] = 0;}
			else{
				if(A >= R){if(A >= D){label[i][j] = 1;sco[i][j] = A;}else{label[i][j] = 2;sco[i][j] = D;}}
				else{if(R >= D){label[i][j] = 3;sco[i][j] = R;}else{label[i][j] = 2;sco[i][j] = D;}}
				if(sco[i][j] > max_sco){max_i = i;max_j = j;max_sco = sco[i][j];}
			}
		}
	}
	int i = max_i;int j = max_j;
	while(1){
		if(label[i][j] == 0){break;}
		else if(label[i][j] == 1){a2b[i-1] = j-1;i--;j--;}
		else if(label[i][j] == 2){i--;}
		else if(label[i][j] == 3){j--;}
	}
	return(a2b);
}

double Falign(double *sco_mtx, int rows, int cols){
	double max_sco = 0;
	double sco[rows+1][cols+1]; memset(sco, 0, sizeof(sco));
	for (int i = 1; i <= rows; i++){
		for (int j = 1; j <= cols; j++){
			double A = sco[i-1][j-1] + sco_mtx[(i-1)*cols+(j-1)];
			double D = sco[i-1][j];
			double R = sco[i][j-1];

			if(A >= R){if(A >= D){sco[i][j] = A;}else{sco[i][j] = D;}}
			else{if(R >= D){sco[i][j] = R;}else{sco[i][j] = D;}}

			if(sco[i][j] > max_sco){max_sco = sco[i][j];}
		}
	}
	return(max_sco);
}


//load contact map
vec_int load_data (string file, int sep_cutoff, mtx_double &mtx, vec_int &vec_div, vec_int &vec, mtx_int &vec_i, mtx_double &prf, int size, vec_char &ss, dist_bin_5_to_20 &rr, double n_bins){
	vec_int n2m;
    	vec_int m2n;
    	string line="";

	int m = 0;n2m.resize(size,-1);
	for(int n = 0; n < size; n++){
		n2m[n] = m;
		m2n.push_back(n);
		m++;
	}
	mtx.resize(m,vector<double>(m,0));
	rr.m5.resize(m,vector<double>(m,0));rr.m5_5.resize(m,vector<double>(m,0));rr.m6.resize(m,vector<double>(m,0));rr.m6_5.resize(m,vector<double>(m,0));rr.m7.resize(m,vector<double>(m,0));rr.m7_5.resize(m,vector<double>(m,0));rr.m8.resize(m,vector<double>(m,0));rr.m8_5.resize(m,vector<double>(m,0));rr.m9.resize(m,vector<double>(m,0));rr.m9_5.resize(m,vector<double>(m,0));rr.m10.resize(m,vector<double>(m,0));rr.m10_5.resize(m,vector<double>(m,0));rr.m11.resize(m,vector<double>(m,0));rr.m11_5.resize(m,vector<double>(m,0));rr.m12.resize(m,vector<double>(m,0));rr.m12_5.resize(m,vector<double>(m,0));rr.m13.resize(m,vector<double>(m,0));rr.m13_5.resize(m,vector<double>(m,0));rr.m14.resize(m,vector<double>(m,0));rr.m14_5.resize(m,vector<double>(m,0));rr.m15.resize(m,vector<double>(m,0));rr.m15_5.resize(m,vector<double>(m,0));rr.m16.resize(m,vector<double>(m,0));rr.m16_5.resize(m,vector<double>(m,0));rr.m17.resize(m,vector<double>(m,0));rr.m17_5.resize(m,vector<double>(m,0));rr.m18.resize(m,vector<double>(m,0));rr.m18_5.resize(m,vector<double>(m,0));rr.m19.resize(m,vector<double>(m,0));rr.m19_5.resize(m,vector<double>(m,0));rr.m20.resize(m,vector<double>(m,0));	
    	ifstream in(file);
	while(getline(in,line)){
		if(isalpha(line[0])){continue;} //ignore seq header
		istringstream is(line);
		int i, j; is >> i >> j;
		double sco;
		double t[31];int index_t=0; //modify later based on distance bin
		if(abs(j-i) >= sep_cutoff){
			while(is >> sco){
				if(sco<0.2){sco=0;}
				if(sco>1){sco=1;} // for native where distance is calculated (<8A)
				t[index_t++]=sco;	
			}	
			rr.m5[n2m[i]][n2m[j]]=t[0];rr.m5_5[n2m[i]][n2m[j]]=t[1];rr.m6[n2m[i]][n2m[j]]=t[2];rr.m6_5[n2m[i]][n2m[j]]=t[3];rr.m7[n2m[i]][n2m[j]]=t[4];rr.m7_5[n2m[i]][n2m[j]]=t[5];rr.m8[n2m[i]][n2m[j]]=t[6];rr.m8_5[n2m[i]][n2m[j]]=t[7];rr.m9[n2m[i]][n2m[j]]=t[8];rr.m9_5[n2m[i]][n2m[j]]=t[9];
			rr.m10[n2m[i]][n2m[j]]=t[10];rr.m10_5[n2m[i]][n2m[j]]=t[11];rr.m11[n2m[i]][n2m[j]]=t[12];rr.m11_5[n2m[i]][n2m[j]]=t[13];rr.m12[n2m[i]][n2m[j]]=t[14];rr.m12_5[n2m[i]][n2m[j]]=t[15];rr.m13[n2m[i]][n2m[j]]=t[16];rr.m13_5[n2m[i]][n2m[j]]=t[17];
			rr.m14[n2m[i]][n2m[j]]=t[18];rr.m14_5[n2m[i]][n2m[j]]=t[19];rr.m15[n2m[i]][n2m[j]]=t[20];rr.m15_5[n2m[i]][n2m[j]]=t[21];rr.m16[n2m[i]][n2m[j]]=t[22];rr.m16_5[n2m[i]][n2m[j]]=t[23];rr.m17[n2m[i]][n2m[j]]=t[24];rr.m17_5[n2m[i]][n2m[j]]=t[25];
			rr.m18[n2m[i]][n2m[j]]=t[26];rr.m18_5[n2m[i]][n2m[j]]=t[27];rr.m19[n2m[i]][n2m[j]]=t[28];rr.m19_5[n2m[i]][n2m[j]]=t[29];rr.m20[n2m[i]][n2m[j]]=t[30];


			rr.m5[n2m[j]][n2m[i]]=t[0];rr.m5_5[n2m[j]][n2m[i]]=t[1];rr.m6[n2m[j]][n2m[i]]=t[2];rr.m6_5[n2m[j]][n2m[i]]=t[3];rr.m7[n2m[j]][n2m[i]]=t[4];rr.m7_5[n2m[j]][n2m[i]]=t[5];rr.m8[n2m[j]][n2m[i]]=t[6];rr.m8_5[n2m[j]][n2m[i]]=t[7];rr.m9[n2m[j]][n2m[i]]=t[8];rr.m9_5[n2m[j]][n2m[i]]=t[9];
                        rr.m10[n2m[j]][n2m[i]]=t[10];rr.m10_5[n2m[j]][n2m[i]]=t[11];rr.m11[n2m[j]][n2m[i]]=t[12];rr.m11_5[n2m[j]][n2m[i]]=t[13];rr.m12[n2m[j]][n2m[i]]=t[14];rr.m12_5[n2m[j]][n2m[i]]=t[15];rr.m13[n2m[j]][n2m[i]]=t[16];rr.m13_5[n2m[j]][n2m[i]]=t[17];        
                        rr.m14[n2m[j]][n2m[i]]=t[18];rr.m14_5[n2m[j]][n2m[i]]=t[19];rr.m15[n2m[j]][n2m[i]]=t[20];rr.m15_5[n2m[j]][n2m[i]]=t[21];rr.m16[n2m[j]][n2m[i]]=t[22];rr.m16_5[n2m[j]][n2m[i]]=t[23];rr.m17[n2m[j]][n2m[i]]=t[24];rr.m17_5[n2m[j]][n2m[i]]=t[25];        
                        rr.m18[n2m[j]][n2m[i]]=t[26];rr.m18_5[n2m[j]][n2m[i]]=t[27];rr.m19[n2m[j]][n2m[i]]=t[28];rr.m19_5[n2m[j]][n2m[i]]=t[29];rr.m20[n2m[j]][n2m[i]]=t[30];
			mtx[n2m[i]][n2m[j]] = t[6];
			mtx[n2m[j]][n2m[i]] = t[6];
		}
	}
	in.close();
	//cout <<"done rr reading.   ";

	for(int i=0; i < mtx.size(); i++){
		vec_i.push_back(vector<int>());
		for(int j=0; j < mtx.size(); j++){	
			if(i == j){
				if(vec_i[i].empty()){vec_div.push_back(0);}
				else{vec_div.push_back(vec_i[i].size());}
			}
			if(mtx[i][j] > 0){ vec_i[i].push_back(j); }
		}
		if(vec_i[i].size() > 0){vec.push_back(i);}
	}
	//cout <<"done load data. ";
	return m2n;
}

double getDistance(point3d & p1, point3d &p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

double getDotProduct(point3d &p1, point3d &p2) {
    return (p1.x * p2.x + p1.y * p2.y + p1.z * p2.z);
}

point3d getCrossProduct(point3d &p1, point3d &p2) {
    point3d p;
    p.x = p1.y * p2.z - p1.z * p2.y;
    p.y = p1.z * p2.x - p1.x * p2.z;
    p.z = p1.x * p2.y - p1.y * p2.x;
    return p;
}

point3d getDifference(point3d &p1, point3d &p2) {
    point3d p;
    p.x = p2.x - p1.x;
    p.y = p2.y - p1.y;
    p.z = p2.z - p1.z;
    return p;
}

double getAngle(point3d &p1, point3d &p2, point3d &p3) {
    double acc = 0.0;
    acc = (p2.x - p1.x) * (p2.x - p3.x) + (p2.y - p1.y) * (p2.y - p3.y) + (p2.z - p1.z) * (p2.z - p3.z);
    double d1 = getDistance(p1, p2);
    double d2 = getDistance(p3, p2);
    if (d1 == 0 || d2 == 0) {
        return 0;
    }
    acc = acc / (d1 * d2);
    if (acc > 1.0)
        acc = 1.0;
    else if (acc < -1.0)
        acc = -1.0;
    acc = acos(acc);
    return acc;
}

double getDihedral(point3d &p1, point3d &p2, point3d &p3, point3d &p4) {
    point3d q, r, s, t, u, v, z;
    double acc, w;
    z.x = z.y = z.z = 0.0;
    q = getDifference(p1, p2);
    r = getDifference(p3, p2);
    s = getDifference(p4, p3);
    t = getCrossProduct(q, r);
    u = getCrossProduct(s, r);
    v = getCrossProduct(u, t);
    w = getDotProduct(v, r);
    acc = getAngle(t, z, u);
    if (w < 0)
        acc = -acc;
    return (acc);
}

map<int,pdbInfo> pdb2coord(string file){
        string line;
        ifstream in(file);
        pdbInfo atom_rec;
        point3d xyz;
        map<int,pdbInfo> all_atoms;
        while(getline(in,line)){
                if (trim(line.substr(0,5).c_str()) == "ATOM" && trim(line.substr(12,4).c_str()) == "N"){
                        xyz.x = stof(trim(line.substr(30,8).c_str()));
                        xyz.y = stof(trim(line.substr(38,8).c_str()));
                        xyz.z = stof(trim(line.substr(46,8).c_str()));
                        atom_rec.n = xyz;
                        continue;
                }
                if (trim(line.substr(0,5).c_str()) == "ATOM" && trim(line.substr(12,4).c_str()) == "CA"){
                        atom_rec.res_num = stoi(trim(line.substr(22,4).c_str()));
                        xyz.x = stof(trim(line.substr(30,8).c_str()));
                        xyz.y = stof(trim(line.substr(38,8).c_str()));
                        xyz.z = stof(trim(line.substr(46,8).c_str()));
                        atom_rec.ca = xyz;
                        if(trim(line.substr(17,3).c_str()) == "GLY"){
                                xyz.x = -99999.0;
                                xyz.y = -99999.0;
                                xyz.z = -99999.0;
                                atom_rec.cb = xyz;
                                all_atoms[atom_rec.res_num] = atom_rec;
                        }
                        continue;
                 
                }
                if (trim(line.substr(0,5).c_str()) == "ATOM" && trim(line.substr(12,4).c_str()) == "CB"){
                        atom_rec.res_num = stoi(trim(line.substr(22,4).c_str()));
                        xyz.x = stof(trim(line.substr(30,8).c_str()));
                        xyz.y = stof(trim(line.substr(38,8).c_str()));
                        xyz.z = stof(trim(line.substr(46,8).c_str()));
                        atom_rec.cb = xyz;
                        all_atoms[atom_rec.res_num] = atom_rec;
                }
        }
        in.close();
        return all_atoms;
}

void read_orientation_phi_query(string file, planar_angles &a,int lenQ){
        string line;
        ifstream in(file);
	a.a_15.resize(lenQ+1,vector<double>(lenQ+1,0));a.a_30.resize(lenQ+1,vector<double>(lenQ+1,0));
        a.a_45.resize(lenQ+1,vector<double>(lenQ+1,0));a.a_60.resize(lenQ+1,vector<double>(lenQ+1,0));
        a.a_75.resize(lenQ+1,vector<double>(lenQ+1,0));a.a_90.resize(lenQ+1,vector<double>(lenQ+1,0));
	a.a_105.resize(lenQ+1,vector<double>(lenQ+1,0));
	a.a_120.resize(lenQ+1,vector<double>(lenQ+1,0));
	a.a_135.resize(lenQ+1,vector<double>(lenQ+1,0));
	a.a_150.resize(lenQ+1,vector<double>(lenQ+1,0));
	a.a_165.resize(lenQ+1,vector<double>(lenQ+1,0));
	a.a_180.resize(lenQ+1,vector<double>(lenQ+1,0));

        while(getline(in,line)){
                if(isalpha(line[0])){continue;}
                istringstream is(line);
		int i, j;
                double prob1, prob2, prob3, prob4, prob5,prob6,prob7,prob8,prob9,prob10,prob11,prob12;
                is >> i >> j >> prob1 >> prob2 >> prob3 >> prob4 >> prob5 >> prob6>>prob7>>prob8>>prob9>>prob10>>prob11>>prob12;
		a.a_15[i][j] = prob1; a.a_30[i][j] = prob2; a.a_45[i][j] = prob3; a.a_60[i][j] = prob4; a.a_75[i][j] = prob5;
                a.a_90[i][j] = prob6; a.a_105[i][j] = prob7; a.a_120[i][j] = prob8; a.a_135[i][j] = prob9; a.a_150[i][j] = prob10;
		a.a_165[i][j] = prob11; a.a_180[i][j] = prob12;
        }
        in.close();
}

void read_orientation_query(string file,dihedral_angles &a,int lenQ){
	string line;
	ifstream in(file);
	a.a_15.resize(lenQ+1,vector<double>(lenQ+1,0));a.a_30.resize(lenQ+1,vector<double>(lenQ+1,0));
        a.a_45.resize(lenQ+1,vector<double>(lenQ+1,0));a.a_60.resize(lenQ+1,vector<double>(lenQ+1,0));
        a.a_75.resize(lenQ+1,vector<double>(lenQ+1,0));a.a_90.resize(lenQ+1,vector<double>(lenQ+1,0));
        a.a_105.resize(lenQ+1,vector<double>(lenQ+1,0));a.a_270.resize(lenQ+1,vector<double>(lenQ+1,0));
        a.a_120.resize(lenQ+1,vector<double>(lenQ+1,0));a.a_285.resize(lenQ+1,vector<double>(lenQ+1,0));
	a.a_135.resize(lenQ+1,vector<double>(lenQ+1,0));a.a_300.resize(lenQ+1,vector<double>(lenQ+1,0));
        a.a_150.resize(lenQ+1,vector<double>(lenQ+1,0));a.a_315.resize(lenQ+1,vector<double>(lenQ+1,0));
        a.a_165.resize(lenQ+1,vector<double>(lenQ+1,0));a.a_330.resize(lenQ+1,vector<double>(lenQ+1,0));
        a.a_180.resize(lenQ+1,vector<double>(lenQ+1,0));a.a_345.resize(lenQ+1,vector<double>(lenQ+1,0));
        a.a_195.resize(lenQ+1,vector<double>(lenQ+1,0));a.a_360.resize(lenQ+1,vector<double>(lenQ+1,0));
	a.a_210.resize(lenQ+1,vector<double>(lenQ+1,0));
	a.a_225.resize(lenQ+1,vector<double>(lenQ+1,0));
	a.a_240.resize(lenQ+1,vector<double>(lenQ+1,0));
	a.a_255.resize(lenQ+1,vector<double>(lenQ+1,0));

	while(getline(in,line)){
		if(isalpha(line[0])){continue;}
		istringstream is(line);
		//cout << line;
		int i, j;
		double prob1, prob2, prob3, prob4, prob5,prob6,prob7,prob8,prob9,prob10,prob11,prob12,prob13,prob14,prob15,prob16,prob17,prob18,prob19,prob20,prob21,prob22,prob23,prob24; 
		is >> i >> j  >> prob1 >> prob2 >> prob3 >> prob4 >> prob5 >> prob6>>prob7>>prob8>>prob9>>prob10>>prob11>>prob12>>prob13>>prob14>>prob15>>prob16>>prob17>>prob18>>prob19>>prob20>>prob21>>prob22>>prob23>>prob24;
		a.a_15[i][j] = prob1; a.a_30[i][j] = prob2; a.a_45[i][j] = prob3; a.a_60[i][j] = prob4; a.a_75[i][j] = prob5;
		a.a_90[i][j] = prob6; a.a_105[i][j] = prob7; a.a_120[i][j] = prob8; a.a_135[i][j] = prob9; a.a_150[i][j] = prob10;
		a.a_165[i][j] = prob11; a.a_180[i][j] = prob12; a.a_195[i][j] = prob13; a.a_210[i][j] = prob14; a.a_225[i][j] = prob15;
		a.a_240[i][j] = prob16; a.a_255[i][j] = prob17;a.a_270[i][j] = prob18; a.a_285[i][j] = prob19; a.a_300[i][j] = prob20;
		a.a_315[i][j] = prob21; a.a_330[i][j] = prob22; a.a_345[i][j] = prob23; a.a_360[i][j] = prob24; 
	}
	in.close();
}

struct temp_bin{
	int bin1=0,bin6=0,bin7=0,bin8=0,bin9=0,bin10=0,bin11=0,bin12=0,bin13=0,bin14=0,bin15=0,bin16=0,bin17=0;
	int bin2=0,bin18=0,bin19=0,bin20=0,bin21=0,bin22=0,bin23=0,bin24=0;
	int bin3=0;
	int bin4=0;
	int bin5=0;
};

temp_bin get_bin(double angle){
        temp_bin b;
	if(angle<15){b.bin1=1;return b;}
	else if(angle < 30){b.bin2=1;return b;}
        else if(angle < 45){b.bin3=1;return b;}
        else if(angle < 60){b.bin4=1;return b;}
        else if(angle < 75){b.bin5=1;return b;}
        else if(angle < 90){b.bin6=1;return b;}
	else if(angle < 105){b.bin7=1;return b;}
	else if(angle < 120){b.bin8=1;return b;}
	else if(angle < 135){b.bin9=1;return b;}
	else if(angle < 150){b.bin10=1;return b;}
	else if(angle < 165){b.bin11=1;return b;}
	else if(angle < 180){b.bin12=1;return b;}
	else if(angle < 195){b.bin13=1;return b;}
	else if(angle < 210){b.bin14=1;return b;}
	else if(angle < 225){b.bin15=1;return b;}
	else if(angle < 240){b.bin16=1;return b;}
	else if(angle < 255){b.bin17=1;return b;}
	else if(angle < 270){b.bin18=1;return b;}
	else if(angle < 285){b.bin19=1;return b;}
	else if(angle < 300){b.bin20=1;return b;}
	else if(angle < 315){b.bin21=1;return b;}
	else if(angle < 330){b.bin22=1;return b;}
	else if(angle < 345){b.bin23=1;return b;}
	else if(angle < 360){b.bin24=1;return b;}
}

bool check_G(point3d &p){
	if (p.x==-99999 and p.y==-99999 and p.z==-99999){return false;}
	else{return true;}
}

/* weights:
 * 6A: 1/(1+exp_fast(10-6))= 0.0181244; 
 * 7A: 1/(1+exp_fast(10-7))=0.0476244; 
 * 8A: 1/(1+exp_fast(10-8))=0.119408; 
 * 9A: 1/(1+exp_fast(10-9))=0.269037
 * 10A: 1/(1+exp_fast(10-10))= 0.5
 * 11A: (1/(1+exp_fast(11-10)))= 0.269037
 * 12A: (1/(1+exp_fast(12-10)))=0.119408
 * 13A: (1/(1+exp_fast(13-10)))=0.0476244
 * 14A: (1/(1+exp_fast(14-10)))=0.0181244
 */

// INITIATE SCORE MATRIX: function for populating the initial similarity matrix
void ini_SCO(double sep_x, double sep_y, mtx_double &SCO,vec_int &vec_a_div,vec_int &vec_b_div,vec_int &vec_a,vec_int &vec_b,mtx_int &vec_a_i,mtx_int &vec_b_i,mtx_double &mtx_a,mtx_double &mtx_b,dist_bin_5_to_20 query_mtx,dist_bin_5_to_20 temp_mtx, dihedral_angles omega,dihedral_angles theta, planar_angles phi,map<int,pdbInfo> &all_atomsT){
	//cout << " INITIATE  ";	

	// Get initial score matrix
	for(int i=0; i < vec_a.size(); i++){ // go through columns (vec_a) in map_a that has contacts
		int ai = vec_a[i];
		for(int j=0; j < vec_b.size(); j++){ // go through columns (vec_b) in map_b that has contacts
			int bi = vec_b[j];
			int A[2] = {(int)vec_a_div[ai],(int)(vec_a_i[ai].size()-vec_a_div[ai])};
			int B[2] = {(int)vec_b_div[bi],(int)(vec_b_i[bi].size()-vec_b_div[bi])};
			for(int k=0; k <= 1; k++){ // left and right of diagonal
				if(A[k] > 0 and B[k] > 0){
					double M[A[k]*B[k]];double omega_sc=0, theta_sc=0, phi_sc=0;
					for(int n=0; n < A[k]; n++){
						int nn = n; if(k == 1){nn += vec_a_div[ai];}
						int aj = vec_a_i[ai][nn];
						//cout << ai << " "<<aj<<"\n";
						int sep_a = abs(ai-aj);
						for(int m=0; m < B[k]; m++){
							int mm = m; if(k == 1){mm += vec_b_div[bi];}
							int bj = vec_b_i[bi][mm];
							int sep_b = abs(bi-bj);
							int sep_D = abs(sep_a-sep_b);
							double sep_M = min(sep_a,sep_b);
							double sep_std = sep_y*(1+pow(sep_M-2,sep_x));
							if(sep_D/sep_std < 6){
//								M[n*B[k]+m] = (0.1*query_mtx.m6[ai][aj]*temp_mtx.m6[bi][bj] + 0.25*mtx_a[ai][aj]*mtx_b[bi][bj] + 0.3*query_mtx.m10[ai][aj]*temp_mtx.m10[bi][bj] + 0.25*query_mtx.m12[ai][aj]*temp_mtx.m12[bi][bj] + 0.1*query_mtx.m14[ai][aj]*temp_mtx.m14[bi][bj])*sepw(sep_M) * gaussian(0,sep_std,sep_D);			
								M[n*B[k]+m] = (0.0181244*query_mtx.m6[ai][aj]*temp_mtx.m6[bi][bj] + 0.0476244*query_mtx.m7[ai][aj]*temp_mtx.m7[bi][bj] + 0.119408*mtx_a[ai][aj]*mtx_b[bi][bj] + 0.269037*query_mtx.m9[ai][aj]*temp_mtx.m9[bi][bj] + 0.5*query_mtx.m10[ai][aj]*temp_mtx.m10[bi][bj] + 0.269037*query_mtx.m11[ai][aj]*temp_mtx.m11[bi][bj] + 0.119408*query_mtx.m12[ai][aj]*temp_mtx.m12[bi][bj] + 0.0476244*query_mtx.m13[ai][aj]*temp_mtx.m13[bi][bj] + 0.0181244*query_mtx.m14[ai][aj]*temp_mtx.m14[bi][bj])* sepw(sep_M) * gaussian(0,sep_std,sep_D);
								
								//try{
								if(check_G(all_atomsT[bi].cb) and check_G(all_atomsT[bj].cb)){
									double omegaT = 180 + (RAD2DEG * getDihedral(all_atomsT[bi].ca, all_atomsT[bi].cb, all_atomsT[bj].cb, all_atomsT[bj].ca));
									double thetaT = 180 + (RAD2DEG * getDihedral(all_atomsT[bi].n, all_atomsT[bi].ca, all_atomsT[bi].cb, all_atomsT[bj].cb));
									double phiT = RAD2DEG * getAngle(all_atomsT[bi].ca, all_atomsT[bi].cb, all_atomsT[bj].cb);
									temp_bin om = get_bin(omegaT);

									omega_sc = (omega.a_15[ai][aj]*om.bin1 + omega.a_30[ai][aj]*om.bin2 + omega.a_45[ai][aj]*om.bin3 + omega.a_60[ai][aj]*om.bin4 + omega.a_75[ai][aj]*om.bin5+omega.a_90[ai][aj]*om.bin6+omega.a_105[ai][aj]*om.bin7+omega.a_120[ai][aj]*om.bin8+omega.a_135[ai][aj]*om.bin9+omega.a_150[ai][aj]*om.bin10+omega.a_165[ai][aj]*om.bin11+omega.a_180[ai][aj]*om.bin12+omega.a_195[ai][aj]*om.bin13+omega.a_210[ai][aj]*om.bin14+omega.a_225[ai][aj]*om.bin15+omega.a_240[ai][aj]*om.bin16+omega.a_255[ai][aj]*om.bin17+omega.a_270[ai][aj]*om.bin18+omega.a_285[ai][aj]*om.bin19+omega.a_300[ai][aj]*om.bin20+omega.a_315[ai][aj]*om.bin21+omega.a_330[ai][aj]*om.bin22+omega.a_345[ai][aj]*om.bin23+omega.a_360[ai][aj]*om.bin24) * sepw(sep_M);
									om = get_bin(thetaT);

									theta_sc=(theta.a_15[ai][aj]*om.bin1+theta.a_30[ai][aj]*om.bin2+theta.a_45[ai][aj]*om.bin3+theta.a_60[ai][aj]*om.bin4+theta.a_75[ai][aj]*om.bin5+theta.a_90[ai][aj]*om.bin6+theta.a_105[ai][aj]*om.bin7+theta.a_120[ai][aj]*om.bin8+theta.a_135[ai][aj]*om.bin9+theta.a_150[ai][aj]*om.bin10+theta.a_165[ai][aj]*om.bin11+theta.a_180[ai][aj]*om.bin12+theta.a_195[ai][aj]*om.bin13+theta.a_210[ai][aj]*om.bin14+theta.a_225[ai][aj]*om.bin15+theta.a_240[ai][aj]*om.bin16+theta.a_255[ai][aj]*om.bin17+theta.a_270[ai][aj]*om.bin18+theta.a_285[ai][aj]*om.bin19+theta.a_300[ai][aj]*om.bin20+theta.a_315[ai][aj]*om.bin21+theta.a_330[ai][aj]*om.bin22+theta.a_345[ai][aj]*om.bin23+theta.a_360[ai][aj]*om.bin24) * sepw(sep_M);
									om = get_bin(phiT);
									phi_sc=(phi.a_15[ai][aj]*om.bin1+phi.a_30[ai][aj]*om.bin2+phi.a_45[ai][aj]*om.bin3+phi.a_60[ai][aj]*om.bin4+phi.a_75[ai][aj]*om.bin5+phi.a_90[ai][aj]*om.bin6+phi.a_105[ai][aj]*om.bin7+phi.a_120[ai][aj]*om.bin8+phi.a_135[ai][aj]*om.bin9+phi.a_150[ai][aj]*om.bin10+phi.a_165[ai][aj]*om.bin11+phi.a_180[ai][aj]*om.bin12) * sepw(sep_M);									//phi = phi_30[ai][aj]*om.bin1 + phi_60[ai][aj]*om.bin2 + phi_90[ai][aj]*om.bin3 + phi_120[ai][aj]*om.bin4 + phi_150[ai][aj]*om.bin5; 
									M[n*B[k]+m] += (omega_sc+theta_sc+phi_sc);
								}

							}else{M[n*B[k]+m] = 0;}
						}
					}
					SCO[ai][bi] += Falign(M,A[k],B[k]);
				}
			}
		}	
	}
}

/*
// MODIFY SCORE MATRIX: function for modifying the initial similarity matrix
*/
vec_int mod_SCO(double do_it, vec_double &gap_a, vec_double &gap_b, double gap_e, mtx_double &SCO, mtx_double &P_SCO,vec_int &vec_a_div,vec_int &vec_b_div,vec_int &vec_a, vec_int &vec_b,mtx_int &vec_a_i, mtx_int &vec_b_i,mtx_double &mtx_a, mtx_double &mtx_b,dist_bin_5_to_20 query_mtx, dist_bin_5_to_20 temp_mtx, dihedral_angles omega, dihedral_angles theta, planar_angles phi,map<int,pdbInfo> &all_atomsT){
	// iterate
	vec_int a2b_tmp;
	for(int it=0; it < do_it; it++){
		// align
		a2b_tmp = align(gap_a,gap_b,gap_e,SCO,P_SCO);
		// update similarity matrix
		double IT = (double)it + 1;
		double s1 = (IT/(IT+1)); double s2 = (1/(IT+1));
		for(int a=0; a < vec_a.size(); a++){ // go through columns (vec_a) in map_a that has contacts
			int ai = vec_a[a];
			for(int b=0; b < vec_b.size(); b++){ // go through columns (vec_b) in map_b that has contacts
				int bi = vec_b[b];
				double sco_contact = 0;double omega_sc=0,theta_sc=0,phi_sc=0;
				for(int n=0; n < vec_a_i[ai].size(); n++){ // go through contacts in vec_a
					int aj = vec_a_i[ai][n];
					int bj = a2b_tmp[aj]; // get mapping
					if(bj != -1){ // if mapping exists
						if((ai > aj and bi > bj) or (ai < aj and bi < bj)){ // if ai-aj in same direction as bi-bj
							double sep_M = min(abs(ai-aj),abs(bi-bj));
//							sco_contact += (0.1*query_mtx.m6[ai][aj]*temp_mtx.m6[bi][bj] + 0.25*mtx_a[ai][aj]*mtx_b[bi][bj] + 0.3*query_mtx.m10[ai][aj]*temp_mtx.m10[bi][bj] + 0.25*query_mtx.m12[ai][aj]*temp_mtx.m12[bi][bj] + 0.1*query_mtx.m14[ai][aj]*temp_mtx.m14[bi][bj])*sepw(sep_M);
							sco_contact += (0.0181244*query_mtx.m6[ai][aj]*temp_mtx.m6[bi][bj] + 0.0476244*query_mtx.m7[ai][aj]*temp_mtx.m7[bi][bj] + 0.119408*mtx_a[ai][aj]*mtx_b[bi][bj] + 0.269037*query_mtx.m9[ai][aj]*temp_mtx.m9[bi][bj] + 0.5*query_mtx.m10[ai][aj]*temp_mtx.m10[bi][bj] + 0.269037*query_mtx.m11[ai][aj]*temp_mtx.m11[bi][bj] + 0.119408*query_mtx.m12[ai][aj]*temp_mtx.m12[bi][bj] + 0.0476244*query_mtx.m13[ai][aj]*temp_mtx.m13[bi][bj] + 0.0181244*query_mtx.m14[ai][aj]*temp_mtx.m14[bi][bj])* sepw(sep_M);	
							
							//try{
							if(check_G(all_atomsT[bi].cb) and check_G(all_atomsT[bj].cb)){
                                                                double omegaT = 180 + (RAD2DEG * getDihedral(all_atomsT[bi].ca, all_atomsT[bi].cb, all_atomsT[bj].cb, all_atomsT[bj].ca));
                                                                temp_bin om = get_bin(omegaT);
								omega_sc = (omega.a_15[ai][aj]*om.bin1 + omega.a_30[ai][aj]*om.bin2 + omega.a_45[ai][aj]*om.bin3 + omega.a_60[ai][aj]*om.bin4 + omega.a_75[ai][aj]*om.bin5+omega.a_90[ai][aj]*om.bin6+omega.a_105[ai][aj]*om.bin7+omega.a_120[ai][aj]*om.bin8+omega.a_135[ai][aj]*om.bin9+omega.a_150[ai][aj]*om.bin10+omega.a_165[ai][aj]*om.bin11+omega.a_180[ai][aj]*om.bin12+omega.a_195[ai][aj]*om.bin13+omega.a_210[ai][aj]*om.bin14+omega.a_225[ai][aj]*om.bin15+omega.a_240[ai][aj]*om.bin16+omega.a_255[ai][aj]*om.bin17+omega.a_270[ai][aj]*om.bin18+omega.a_285[ai][aj]*om.bin19+omega.a_300[ai][aj]*om.bin20+omega.a_315[ai][aj]*om.bin21+omega.a_330[ai][aj]*om.bin22+omega.a_345[ai][aj]*om.bin23+omega.a_360[ai][aj]*om.bin24) * sepw(sep_M);

								double thetaT = 180 + (RAD2DEG * getDihedral(all_atomsT[bi].n, all_atomsT[bi].ca, all_atomsT[bi].cb, all_atomsT[bj].cb));
								double phiT = RAD2DEG * getAngle(all_atomsT[bi].ca, all_atomsT[bi].cb, all_atomsT[bj].cb);
								om = get_bin(thetaT);
								theta_sc=(theta.a_15[ai][aj]*om.bin1+theta.a_30[ai][aj]*om.bin2+theta.a_45[ai][aj]*om.bin3+theta.a_60[ai][aj]*om.bin4+theta.a_75[ai][aj]*om.bin5+theta.a_90[ai][aj]*om.bin6+theta.a_105[ai][aj]*om.bin7+theta.a_120[ai][aj]*om.bin8+theta.a_135[ai][aj]*om.bin9+theta.a_150[ai][aj]*om.bin10+theta.a_165[ai][aj]*om.bin11+theta.a_180[ai][aj]*om.bin12+theta.a_195[ai][aj]*om.bin13+theta.a_210[ai][aj]*om.bin14+theta.a_225[ai][aj]*om.bin15+theta.a_240[ai][aj]*om.bin16+theta.a_255[ai][aj]*om.bin17+theta.a_270[ai][aj]*om.bin18+theta.a_285[ai][aj]*om.bin19+theta.a_300[ai][aj]*om.bin20+theta.a_315[ai][aj]*om.bin21+theta.a_330[ai][aj]*om.bin22+theta.a_345[ai][aj]*om.bin23+theta.a_360[ai][aj]*om.bin24) * sepw(sep_M);
								om = get_bin(phiT);
								phi_sc=(phi.a_15[ai][aj]*om.bin1+phi.a_30[ai][aj]*om.bin2+phi.a_45[ai][aj]*om.bin3+phi.a_60[ai][aj]*om.bin4+phi.a_75[ai][aj]*om.bin5+phi.a_90[ai][aj]*om.bin6+phi.a_105[ai][aj]*om.bin7+phi.a_120[ai][aj]*om.bin8+phi.a_135[ai][aj]*om.bin9+phi.a_150[ai][aj]*om.bin10+phi.a_165[ai][aj]*om.bin11+phi.a_180[ai][aj]*om.bin12) * sepw(sep_M); 
								sco_contact += (omega_sc+theta_sc+phi_sc);			
							}
						
														
						}
					}
				}
				
                                SCO[ai][bi] = s1*SCO[ai][bi] + s2*sco_contact;
			}
		}
	}
	//cout << "  done.. ";
	return(a2b_tmp);
}

/*
// CHK: compute number of contacts/gaps made
*/

void chk (vec_double &gap_a, vec_double &gap_b, double &gap_e_w, double& con_sco,double& gap_sco,double& prf_sco,vec_int& vec_a_div,mtx_int& vec_a_i,mtx_double& mtx_a,mtx_double& mtx_b,vec_int& a2b,mtx_double &P_SCO,dist_bin_5_to_20 query_mtx,dist_bin_5_to_20 temp_mtx, dihedral_angles omega, dihedral_angles theta, planar_angles phi,map<int,pdbInfo> &all_atomsT){
//	cout << " chk ..";
	int size_a = mtx_a.size();
	int a = 0;int b = 0;
	for(int ai = 0; ai < size_a; ai++){
		int bi = a2b[ai];
		if(bi != -1){
			prf_sco += P_SCO[ai][bi];
			if(a > 0){ // compute number of gaps
				double num_gap_a = ((ai-a)-1); if(num_gap_a > 0){gap_sco += gap_a[ai] + gap_a[ai] * gap_e_w * (num_gap_a-1);}
				double num_gap_b = ((bi-b)-1); if(num_gap_b > 0){gap_sco += gap_b[bi] + gap_b[bi] * gap_e_w * (num_gap_b-1);}
			}
			for(int m=0; m < vec_a_div[ai]; m++){ // compute number of contacts
				int aj = vec_a_i[ai][m];
				int bj = a2b[aj];
				if(bj != -1){
					double sep_M = min(abs(ai-aj),abs(bi-bj));double omega_sc=0,theta_sc=0,phi_sc=0;
//					con_sco += (0.1*query_mtx.m6[ai][aj]*temp_mtx.m6[bi][bj] + 0.25*mtx_a[ai][aj]*mtx_b[bi][bj] + 0.3*query_mtx.m10[ai][aj]*temp_mtx.m10[bi][bj] + 0.25*query_mtx.m12[ai][aj]*temp_mtx.m12[bi][bj] + 0.1*query_mtx.m14[ai][aj]*temp_mtx.m14[bi][bj])*sepw(sep_M);
					con_sco += (0.0181244*query_mtx.m6[ai][aj]*temp_mtx.m6[bi][bj] + 0.0476244*query_mtx.m7[ai][aj]*temp_mtx.m7[bi][bj] + 0.119408*mtx_a[ai][aj]*mtx_b[bi][bj] + 0.269037*query_mtx.m9[ai][aj]*temp_mtx.m9[bi][bj] + 0.5*query_mtx.m10[ai][aj]*temp_mtx.m10[bi][bj] + 0.269037*query_mtx.m11[ai][aj]*temp_mtx.m11[bi][bj] + 0.119408*query_mtx.m12[ai][aj]*temp_mtx.m12[bi][bj] + 0.0476244*query_mtx.m13[ai][aj]*temp_mtx.m13[bi][bj] + 0.0181244*query_mtx.m14[ai][aj]*temp_mtx.m14[bi][bj])* sepw(sep_M);
					//try{
					if(check_G(all_atomsT[bi].cb) and check_G(all_atomsT[bj].cb)){
						 double omegaT = 180 + (RAD2DEG * getDihedral(all_atomsT[bi].ca, all_atomsT[bi].cb, all_atomsT[bj].cb, all_atomsT[bj].ca));  
                                                 temp_bin om = get_bin(omegaT);
						 omega_sc = (omega.a_15[ai][aj]*om.bin1 + omega.a_30[ai][aj]*om.bin2 + omega.a_45[ai][aj]*om.bin3 + omega.a_60[ai][aj]*om.bin4 + omega.a_75[ai][aj]*om.bin5+omega.a_90[ai][aj]*om.bin6+omega.a_105[ai][aj]*om.bin7+omega.a_120[ai][aj]*om.bin8+omega.a_135[ai][aj]*om.bin9+omega.a_150[ai][aj]*om.bin10+omega.a_165[ai][aj]*om.bin11+omega.a_180[ai][aj]*om.bin12+omega.a_195[ai][aj]*om.bin13+omega.a_210[ai][aj]*om.bin14+omega.a_225[ai][aj]*om.bin15+omega.a_240[ai][aj]*om.bin16+omega.a_255[ai][aj]*om.bin17+omega.a_270[ai][aj]*om.bin18+omega.a_285[ai][aj]*om.bin19+omega.a_300[ai][aj]*om.bin20+omega.a_315[ai][aj]*om.bin21+omega.a_330[ai][aj]*om.bin22+omega.a_345[ai][aj]*om.bin23+omega.a_360[ai][aj]*om.bin24) * sepw(sep_M);						 
						 double thetaT = 180 + (RAD2DEG * getDihedral(all_atomsT[bi].n, all_atomsT[bi].ca, all_atomsT[bi].cb, all_atomsT[bj].cb));
						 double phiT = RAD2DEG * getAngle(all_atomsT[bi].ca, all_atomsT[bi].cb, all_atomsT[bj].cb);
						 om = get_bin(thetaT);
						 theta_sc=(theta.a_15[ai][aj]*om.bin1+theta.a_30[ai][aj]*om.bin2+theta.a_45[ai][aj]*om.bin3+theta.a_60[ai][aj]*om.bin4+theta.a_75[ai][aj]*om.bin5+theta.a_90[ai][aj]*om.bin6+theta.a_105[ai][aj]*om.bin7+theta.a_120[ai][aj]*om.bin8+theta.a_135[ai][aj]*om.bin9+theta.a_150[ai][aj]*om.bin10+theta.a_165[ai][aj]*om.bin11+theta.a_180[ai][aj]*om.bin12+theta.a_195[ai][aj]*om.bin13+theta.a_210[ai][aj]*om.bin14+theta.a_225[ai][aj]*om.bin15+theta.a_240[ai][aj]*om.bin16+theta.a_255[ai][aj]*om.bin17+theta.a_270[ai][aj]*om.bin18+theta.a_285[ai][aj]*om.bin19+theta.a_300[ai][aj]*om.bin20+theta.a_315[ai][aj]*om.bin21+theta.a_330[ai][aj]*om.bin22+theta.a_345[ai][aj]*om.bin23+theta.a_360[ai][aj]*om.bin24) * sepw(sep_M);						
						 om = get_bin(phiT);
						 phi_sc=(phi.a_15[ai][aj]*om.bin1+phi.a_30[ai][aj]*om.bin2+phi.a_45[ai][aj]*om.bin3+phi.a_60[ai][aj]*om.bin4+phi.a_75[ai][aj]*om.bin5+phi.a_90[ai][aj]*om.bin6+phi.a_105[ai][aj]*om.bin7+phi.a_120[ai][aj]*om.bin8+phi.a_135[ai][aj]*om.bin9+phi.a_150[ai][aj]*om.bin10+phi.a_165[ai][aj]*om.bin11+phi.a_180[ai][aj]*om.bin12) * sepw(sep_M);
                                                 con_sco += (omega_sc+theta_sc+phi_sc);
                                        }
						
				}
			}
			a = ai;b = bi;
		}
	}
	gap_sco /= 2;
}

// compute profile similarity matrix
void ini_prf_SCO(mtx_double &P_SCO, double &prf_w, mtx_double &prf_a, vec_char &aa_a, mtx_double &prf_b, vec_char &aa_b){
    	int size_a = prf_a.size();
    	int size_b = prf_b.size();
    
    	P_SCO.resize(size_a,vector<double>(size_b,0));

	// compute background frequencies
	vec_double pb(20,0); int prf_size = prf_a[0].size();
	double pb_size = 0;
	for(int ai = 0; ai < size_a; ai++){
		if(aa_a[ai] != 'X'){ // ignore positions that have no identity
			for(int p=0; p < prf_size; p++){pb[p] += prf_a[ai][p];}
			pb_size += 1;
		}
	}
	for(int bi = 0; bi < size_b; bi++){
		if(aa_b[bi] != 'X'){ // ignore positions that have no identity
			for(int p=0; p < prf_size; p++){pb[p] += prf_b[bi][p];}
			pb_size += 1;
		}
	}
	for (int i=0; i < size_a; i++){
		for (int j=0; j < size_b; j++){
			if(aa_a[i] == 'X' or aa_b[j] == 'X'){P_SCO[i][j] = 0;} // if no identity, return score of 0
			else{
				// profile comparison calculation, similar to HHsuite from Soeding.
				double tmp_sco = 0;
				for(int p=0; p < prf_size; p++){tmp_sco += (prf_a[i][p]*prf_b[j][p])/(pb[p]/pb_size);}
				P_SCO[i][j] = log2(tmp_sco)/5 * prf_w;
			}
		}
	}	
}


//Needleman Wunsch global alignment
double NeedlemanWunsch(string f1, string f2, int gap_open,int gap_extn){

        int imut[24][24]={0};
        Blosum62Matrix(imut); // Read Blosum scoring matrix 
        string seqW = "*ARNDCQEGHILKMFPSTWYVBZX"; // Amino acid order in the BLAST's scoring matrix (e.g.,Blosum62). 

	int f1_len = f1.length();
	int f2_len = f2.length();
        int seq1[f1_len]; // seq1 and seq2 are arrays that store the amino acid order numbers of sequence1 and sequence2.
	std::fill_n(seq1, f1_len, 0);
        int seq2[f2_len]; // For example, 1 stand for A, 2 represent R and etc.
	std::fill_n(seq2, f2_len, 0);

        int i,j;
	for(i=1;i<f1_len;i++){
                for(j=1;j<seqW.length();j++){
                        if(f1[i] == seqW[j]){
                                seq1[i]=j;
                        }
                }
        }
        for(i=1;i<f2_len;i++){
                for(j=1;j<seqW.length();j++){
                        if(f2[i] == seqW[j]){
                                seq2[i]=j;
                        }
                }
        }

        int score[f1_len][f2_len]; //// score[i][j] stands for the alignment score that align ith position of the first sequence to the jth position of the second sequence.
	memset(score, 0, f1_len * f2_len * sizeof(int));
        for(i=1;i<f1_len;i++){
                for(j=1;j<f2_len;j++){
                        score[i][j] = imut[seq1[i]][seq2[j]];
                }
        }

        int j2i[f2_len+1];
	std::fill_n(j2i, f2_len+1, 0);

        for(j=1;j<f2_len;j++){j2i[j] = -1;} // !all are not aligned

        int val[f1_len+1][f2_len+1]; // val[][] was assigned as a global variable, and the value could be printed in the final.
	memset(val, 0, (f1_len+1) * (f2_len+1) * sizeof(int));
        int idir[f1_len+1][f2_len+1];
	memset(idir, 0, (f1_len+1) * (f2_len+1) * sizeof(int));
        int preV[f1_len+1][f2_len+1];
	memset(preV, 0, (f1_len+1) * (f2_len+1) * sizeof(int));
        int preH[f1_len+1][f2_len+1];
	memset(preH, 0, (f1_len+1) * (f2_len+1) * sizeof(int));

        int D,V,H;
        bool standard = true;
        if(standard) // This is a standard Needleman-Wunsch dynamic program (by Y. Zhang 2005).
	{
		int jpV[f1_len+1][f2_len+1];
		memset(jpV, 0, (f1_len+1) * (f2_len+1) * sizeof(int));
                int jpH[f1_len+1][f2_len+1];
		memset(jpH, 0, (f1_len+1) * (f2_len+1) * sizeof(int));

                val[0][0]=0;
		for(i=1;i<f1_len;i++){
                        val[i][0] = gap_extn*i;
                        preV[i][0] = val[i][0];
                        idir[i][0] = 0;
                        jpV[i][0] = 1;
                        jpH[i][0] = i;
                }
                for(j=1;j<f2_len;j++){
                        val[0][j] = gap_extn*j;
                        preH[0][j] = val[0][j];
                        idir[0][j] = 0;
                        jpV[0][j] = j;
                        jpH[0][j] = 1;
                }
		//DP
		for(j=1;j<f2_len;j++){
                        for(i=1;i<f1_len;i++){
                                D = val[i-1][j-1] + score[i][j];        // from diagonal, val(i,j) is val(i-1,j-1) 
				jpH[i][j] = 1;
                                int val1 = val[i-1][j] + gap_open;  // gap_open from both D and V
                                int val2 = preH[i-1][j] + gap_extn; // gap_extn from horizontal
                                if(val1>val2){H = val1;}   // last step from D or V             
                                else{H = val2;if(i > 1){jpH[i][j] = jpH[i-1][j] + 1;}} // last step from H and record long-gap
				
				jpV[i][j] = 1;
                                val1 = val[i][j-1] + gap_open;
                                val2 = preV[i][j-1] + gap_extn;
                                if(val1>val2){V = val1;}
                                else{V = val2;if(j > 1){jpV[i][j] = jpV[i][j-1] + 1;}}

				preH[i][j] = H; // unaccepted H
                                preV[i][j] = V; // unaccepted V         
                                if((D>H)&&(D>V)){
                                        idir[i][j]=1;
                                        val[i][j]=D;
                                }
                                else if(H > V){
                                        idir[i][j] = 2;
                                        val[i][j] = H;
                                }
                                else{
                                        idir[i][j] = 3;
                                        val[i][j] = V;
                                }
                        }
                }
		/* traceback the pathway */
		i = f1_len-1;
                j = f2_len-1;
                while((i>0)&&(j>0)){
                        if(idir[i][j]==1)       // from diagonal
                        {
                                j2i[j] = i;
                                i=i-1;
                                j=j-1;
                        }
                        else if(idir[i][j]==2)  // from horizonal
                        {
                                int temp1 = jpH[i][j];
                                for(int me=1;me<=temp1;me++) {
                                        if(i>0){i=i-1;}
                                }
                        }
                        else {
                                int temp2 = jpV[i][j];
                                for(int me=1;me<=temp2;me++) {
                                        if(j>0){j=j-1;}
                                }
                        }
                }
        }		
	int L_id=0;
        int L_ali=0;
        for(j=1;j<f2_len;j++){
                if(j2i[j]>0){
                        i=j2i[j];
                        L_ali=L_ali+1;
                        if(seq1[i]==seq2[j]){L_id=L_id+1;}
                }
        }

        double identity = L_id*1.0/(f2_len-1);
	return identity;
}


string trim(const string& str)
{
    	size_t first = str.find_first_not_of(' ');
    	if (string::npos == first)
    	{
        	return str;
    	}
    	size_t last = str.find_last_not_of(' ');
    	return str.substr(first, (last - first + 1));
}

void read_template_ss_sa_angles(string dssp, int ss[], float sa[], float psi[], float phi[]){
	ifstream in(dssp);
	if(!in) {
		printf("Error: cannot open dssp file '%s'\n");
	}
	char buf[200];
	for(; !in.eof(); ) {
		in.get(buf, 4);
		in.ignore(INT_MAX, '\n');
		if(!strcmp(buf, "  #")) break;
	}
	string line = "", ssSeq = "", saSeq = "";
	char aaSeq;
	string resNum = "";
	int resNumInt = 0;
	while(getline(in,line)){
		resNum = trim(line.substr(5,5).c_str());
		try {
                        resNumInt = stoi(resNum);
		} catch (...){	
			continue;
		}
		aaSeq = line[13];
		ssSeq = line[16];
		saSeq = line.substr(35,3).c_str();
		ss[resNumInt] = get8to3ss(ssSeq);
		sa[resNumInt] = getSolAcc(aaSeq,saSeq);
		phi[resNumInt] = stof(line.substr(103,6).c_str());
		psi[resNumInt] = stof(line.substr(109,6).c_str());
	}
	in.close();
}

// align query-template
ALIGN1_RESULT align1(int seq_num, string template_name, string x, string y, double mymtx[][21], double umtxb[][21], int mysec[],
		float mysa[], float mypsi[], float myphi[], int iysec[], float iysa[], float iypsi[], float iyphi[], 
		double myprf[][21], double prflib[][21], float p[]) {

	// mymtx the MTX profile of query sequence
	// umtxb the MTX profile of template sequence
		
	int x_len = x.length();
	int y_len = y.length();

	double val[x_len+1][y_len+1];
	memset(val, 0, (x_len+1) * (y_len+1) * sizeof(double)); 
	double idir[x_len+1][y_len+1];
	memset(idir, 0, (x_len+1) * (y_len+1) * sizeof(double));
	double preV[x_len+1][y_len+1];
	memset(preV, 0, (x_len+1) * (y_len+1) * sizeof(double));
	double preH[x_len+1][y_len+1];
	memset(preH, 0, (x_len+1) * (y_len+1) * sizeof(double));
	double preD[x_len+1][y_len+1];
	memset(preD, 0, (x_len+1) * (y_len+1) * sizeof(double));
	double idirH[x_len+1][y_len+1];
	memset(idirH, 0, (x_len+1) * (y_len+1) * sizeof(double));
	double idirV[x_len+1][y_len+1];
	memset(idirV, 0, (x_len+1) * (y_len+1) * sizeof(double));

	//cout << "ccc parameters ---------------------->\n";
	double wgap0 = p[1]; // !gap_open
	double wgap1 = p[2]; // !gap_extn
	double wgap2 = p[3]; // !secondary structure match
	double wgap3 = p[4]; // ! solvent accessibility
	double shft = p[5]; // !constant shift
	
	double gap0[4][4] = {0};
	double gap1[4][4] = {0};
	double gap2[4][4] = {0};

	double gap00[3][3] = {0};
	double gap11[3][3] = {0};
	double gap3[3][3] = {0};
	
	int i, j;
	for (i = 1; i <= 3; i++) {
		for (j = 1; j <= 3; j++) {
			gap0[i][j] = wgap0;
			gap1[i][j] = wgap1;
			gap2[i][j] = wgap2;
		}
	}
	for (i = 0; i <= 1; i++) {
		for (j = 0; j <= 1; j++) {
			gap00[i][j] = wgap0;
			gap11[i][j] = wgap1;
			gap3[i][j] = wgap3;
		}
	}

	gap0[1][1] = 100; // !gap_open is big for ss regions
	gap0[2][2] = 100;
	gap1[1][1] = 100;
	gap1[2][2] = 100;

	float mygap = 10;

	gap00[0][0] = 100; // !gap_open is big for sa regions
	gap00[1][1] = 100;
	gap00[2][2] = 100;
	gap11[0][0] = 100; 
	gap11[1][1] = 100;
	gap11[2][2] = 100;

	gap2[1][1] = -wgap2; // !secondary structure match
	gap2[2][2] = -wgap2;
	gap2[3][3] = -wgap2;

	gap3[0][0] = -wgap3; // !solvenet accessibility match
	gap3[1][1] = -wgap3;
	gap3[2][2] = -wgap3;
	gap3[0][1] = -wgap2 * 0.5;
	gap3[1][0] = -wgap2 * 0.5;
	gap3[1][2] = -wgap2 * 0.5;
	gap3[2][1] = -wgap2 * 0.5;

	double gap4 = p[6];
	double w4 = p[7]; // ! for phi angle
	double gap6 = p[8]; // ! not used
	double w1 = p[9]; // !always 1
	double w2 = p[10]; // !weight for depth library
	double w3 = p[11]; // !weight for psi angle

	// cccc hydrophobic matrix
	double Y1[21][21] = {0};
	double gap000[21][21] = {0};
	double gap111[21][21] = {0};

	for (i = 1; i <= 20; i++) {
		for (j = 1; j <= 20; j++) {
			Y1[i][j] = 0;
			if ((i >= 2) && (i <= 8) && (j >= 2) && (j <= 8)) {
				gap000[i][j] = 100;
				gap111[i][j] = 100;
				Y1[i][j] = -gap4;
			} else if (((i == 10) || (i == 20)) && ((j == 10) || (j == 20)) && (i == j)) {
				// ccccccc if GLY and PRO residues are equal
				gap000[i][j] = 100;
				gap111[i][j] = 100;
				Y1[i][j] = -gap4;
			} else if (i == j) {	
				gap000[i][j] = wgap0;
				gap111[i][j] = wgap1;
				Y1[i][j] = -gap4 * 0.7;
			}
		}
	}
	// ccc initializations --------------->
	val[0][0] = 0.0;
	for (i = 1; i < x.length(); i++) {
		val[i][0] = 0;
		idir[i][0] = 0;
		preD[i][0] = 0;
		preH[i][0] = 1000.0;
		preV[i][0] = 1000.0;
	}
	for (j = 1; j < y.length(); j++) {
		val[0][j] = 0;
		idir[0][j] = 0;
		preD[0][j] = 0;
		preH[0][j] = 1000.0;
		preV[0][j] = 1000.0;
	}

	for (i = 1; i < x_len; i++) {
		int is1 = mysec[i];
		float isa1 = mysa[i];   // !predicted solvenet accessibility for query seq
		float ipsi1 = mypsi[i]; // !predicted psi angle from query sequence
		float iphi1 = myphi[i]; // !predicted phi angle from query sequence
		int ir1 = getResNum(x[i]); //check!!
		int ss_score = 0;
		for (j = 1; j < y_len; j++) {
			if (is1 == iysec[j]) {
				ss_score = 1;
			} else {
				ss_score = -1;
			}

			int is2 = iysec[j]; //// !SS of tempates
			float isa2 = iysa[j]; // !solvenent accessibility of tempates
			float ipsi2 = iypsi[j]; // ! psi angel from templates
			float iphi2 = iyphi[j]; // ! phi angel from templates
			int ir2 = getResNum(y[j]);

			double term1 = 0.0;
			double term2 = 0.0;

			for (int k = 1; k <= 20; k++) {
				term1 = term1 - umtxb[j][k] * myprf[i][k]; //myprf is sequence profile (protein.prf)
				term2 = term2 - prflib[j][k] * mymtx[i][k]; //mymtx is query mtx profile (protein.mtx)
			}

			double tmp_psi = abs(ipsi1 - ipsi2);		
			double tmp_phi = abs(iphi1 - iphi2);

			double terms = term1 + term2 * w2 + gap2[is1][is2] + wgap3 * (2 * abs(isa1 - isa2) - 1.0) + Y1[ir1][ir2]
					+ shft + w3 * (2*(tmp_psi/360.0) - 1.0) + w4 * (2 * (tmp_phi / 360.0) - 1.0) - (0.50 * ss_score); // !energy 

			preD[i][j] = val[i - 1][j - 1] + terms; // !from diagonal

			// ccc preH: pre-accepted H ------------->
			double D = 0, H = 0, V = 0;
			if (j == (y_len - 1)) {
				D = preD[i - 1][j];
				H = preH[i - 1][j];
				V = preV[i - 1][j];
			} else {
				D = preD[i - 1][j] + gap0[is1][is2];
				H = preH[i - 1][j] + gap1[is1][is2];
				V = preV[i - 1][j] + gap1[is1][is2];
			}

			if ((D < H) && (D < V)) {		
				preH[i][j] = D;
				idirH[i - 1][j] = 1;
			} else if (H < V) {
				preH[i][j] = H;
				idirH[i - 1][j] = 2;
			} else {
				preH[i][j] = V;
				idirH[i - 1][j] = 3;
			}

			// ccc preV: pre-accepted V ------------>			
			if (i == (x_len - 1)) {
				D = preD[i][j - 1];
				H = preH[i][j - 1];
				V = preV[i][j - 1];
			} else {
				D = preD[i][j - 1] + gap0[is1][is2];
				H = preH[i][j - 1] + gap1[is1][is2];
				V = preV[i][j - 1] + gap1[is1][is2];
			}

			if ((D < H) && (D < V)) {
				preV[i][j] = D;
				idirV[i][j - 1] = 1;
			} else if (H < V) {
				preV[i][j] = H;
				idirV[i][j - 1] = 2;
			} else {
				preV[i][j] = V;
				idirV[i][j - 1] = 3;
			}

			// ccc decide idir(i,j)-------------->
			if ((preD[i][j] < preH[i][j]) && (preD[i][j] < preV[i][j])) {
				idir[i][j] = 1;
				val[i][j] = preD[i][j];
			} else if (preH[i][j] < preV[i][j]) {
				idir[i][j] = 2;
				val[i][j] = preH[i][j];
			} else {
				idir[i][j] = 3;
				val[i][j] = preV[i][j];
			}
		}
	}
	/*connected similarity*/
        for (i = 2; i < x_len-1; i++) {
                for (j = 2; j < y_len-1; j++) {
                        val[i][j] += 0.5*(val[i-1][j-1]+val[i+1][j+1]);
                }
        }
        val[x_len-1][y_len-1] += 0.5*val[x_len-2][y_len-2];

	// ccccc tracing back the pathway:
	string alignX = "";
	string alignY = "";

	int lena = 0; // alignment length
	int iii = 0;
	int nn = 0;

	i = x_len - 1;
	j = y_len - 1;

	while ((i > 0) && (j > 0)) {
		lena++;
		if (idir[i][j] == 1) // !from diagonal
		{
			alignX = x[i] + alignX;
			alignY = y[j] + alignY;
			i = i - 1;
			j = j - 1;
			iii = iii + 1; // !no gaps
			if (nn == 0) {nn = 1;} // !first aligned pair
		} else if (idir[i][j] == 2) // !from horizonal
		{
			alignX = x[i] + alignX;
			alignY = "-" + alignY;
			i = i - 1;
			idir[i][j] = idirH[i][j];
			iii = iii + 1; // !template gaps
		} else {
			alignX = "-" + alignX;
			alignY = y[j] + alignY;
			j = j - 1;
			idir[i][j] = idirV[i][j];
			if (nn == 1) { iii = iii + 1;} // !query gaps
		}
	}

	if ((j <= 0) && (i > 0)) {
		for (int k = i; k >= 1; k--) {
			lena++;
			alignX = x[k] + alignX;
			alignY = "-" + alignY;
			iii = iii + 1;
		}
	}	

	if ((i <= 0) && (j > 0)) {
		for (int k = j; k >= 1; k--) {
			lena++;
			alignX = "-" + alignX;
			alignY = y[k] + alignY;
		}
	}

	double score = val[x_len - 1][y_len - 1];

	double idt = 0.0;
	int ialg = 0;
	for (i = 1; i < alignX.length(); i++) {
		if (alignX[i] == alignY[i]) {
			idt++;
		}
		if ((alignX[i] != '-') && (alignY[i] != '-')) {
			ialg++;
		}
	}

	idt = idt * 100 * 1.0 / (x_len - 1);

	ALIGN1_RESULT result;
	result.alignment = alignY + " " + alignX;
	result.seq_number = seq_num;
	result.template_name = template_name;
	result.score_align = score;
	result.idt = int(idt);
	result.ialg = ialg;
	result.lena = lena;
	result.iii = iii;
	return result;
}	

int getResNum(char a) {
	int rescod[26] = {9, -1, 1, 16, 15, 3, 10, 17, 4, -1, 19, 5,2, 14, -1, 20, 13, 18, 12, 11, -1, 6, 7, -1, 8, -1};
	int ord = int(a); //integer representing the Unicode code 
	int ord_A = int('A');
	return rescod[ord - ord_A];
}

void read_DEP_profile(string dep_prf, double dep[][21], string seq){
	ifstream br(dep_prf); 
	string line = "";
	int count = 0;
	int m_count = 0;
	int r_num = 0;
	int i = 0;
	int seq_index = 1;

	while(getline(br,line)){
                count++;
                if (count == 1)
                {
                        continue;
                } else if (line.find(seq[seq_index]) != std::string::npos){ //return npos if nothing is found
                        r_num++;
                        seq_index++;
                        m_count = 1;
                        istringstream wds(line);
                        istream_iterator<std::string> begin(wds);
                        istream_iterator<std::string> end;
                        std::vector<std::string> vstrings(begin, end);
                        if (vstrings.size() >= 1){
                                for (i = 3; i < vstrings.size(); i++) {
                                        dep[r_num][m_count] = stod(vstrings[i]);
                                        m_count++;
                                }
                        } else {
                                cout << "\nError: " << line;
                                exit(0);
                        }

                } else {
                        istringstream wds(line);
                        istream_iterator<std::string> begin(wds);
                        istream_iterator<std::string> end;
                        std::vector<std::string> vstrings(begin, end);
                        if (vstrings.size() > 0){
                                for (i = 0; i < vstrings.size(); i++) {
                                        dep[r_num][m_count] = stod(vstrings[i]);
                                        m_count++;
                                }
                        } else {
                                continue;
                        }
                }
        }
        br.close();
}

string read_fasta(string file) {
	ifstream in(file); // open a file for reading
	string line = "";
	string x = "";
	while(getline(in,line)){
		if(line[0] == '>') {
			continue;
		} else {
			x = x + line;
		}
	}
	in.close();
	return ("*" + x);
}	
		
int get8to3ss(string ssParm) {
	int eTo3 = 3;
	if (ssParm == "H" || ssParm == "G" || ssParm == "I") {
		eTo3 = 1;
	} else if (ssParm == "E" || ssParm == "B") {
		eTo3 = 2;
	} else {
		eTo3 = 3;
	}
	return eTo3;
}
	
float getSolAcc(char aa, string sa){
	char aaSA[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
	int aaSaVal[] = {115, 135, 150, 190, 210, 75, 195, 175, 200, 170, 185, 160, 145, 180, 225, 115, 140, 155, 255, 230};
	int k = 0;
	float sVal = 0.0; //for testing if dssp file has some entry like '!' (3lqzB.dssp)
	while(k < 20){ // 20 is num of available AA
		if(aa == aaSA[k]){
			sVal = atof(sa.c_str()) / (float) aaSaVal[k];
			break;
		} else{
			k = k + 1;
		}
	}
	return sVal;
}

void read_query_ss_sa_angles(string file, string x, int mysec[], float mysa[], float mypsi[], float myphi[]) { 
	ifstream in(file);
	ofstream fw;
        fw.open("seq.ss");
	string secon_ss = "HEC";
	string line = "";
	getline(in, line); //skip header
	while(getline(in,line)){
		istringstream is(line);
		char aa;
       	 	string ss = "",sa = "",phi = "",psi = "",count = "0";
        	is >> count >> aa >> ss >> sa >> phi >> psi;
		mysec[atoi(count.c_str())] = get8to3ss(ss);
		fw << secon_ss[mysec[atoi(count.c_str())] - 1];
		mysa[atoi(count.c_str())] = getSolAcc(aa,sa);
		myphi[atoi(count.c_str())] = atof(phi.c_str()); //atoi: convert count (string) to int
		mypsi[atoi(count.c_str())] = atof(psi.c_str());
	}
	in.close();
	fw.close();
}

void read_mtx(string file, double mtx[][21]){
	int i, j;	
	double bls[26+1] = {0};
	double mtx0[21] = {0};
	int map2[] = { -1, 4, 13, 7, 10, 12, 20, 21, 23, 2, 8, 19, 18, 16, 14, 6, 5, 9, 17, 11, 15 };
	string line = "";
	int count = 0;
	int num = 1;
	ifstream in(file);
	while(getline(in,line)){
		count++;
		if (count > 14) {
			istringstream is(line);
			string wds;
			for (i = 1; i <= 26; i++) {
				is >> wds;
				bls[i] = atof(wds.c_str());
			}
			for (i = 1; i <= 20; i++) {
				j = map2[i];
				mtx0[i] = bls[j] * 0.01;
				mtx[num][i] = mtx0[i];
			}
			num++;
		}
	}
	in.close();
}

void read_prf(string file, double mtx[][21]){
	string line = "";
	int count = 0;
	int num = 1;
	ifstream in(file);
	while(getline(in,line)){
		count++;
		if (count > 1) {
			istringstream is(line);
			string wds;
			is >> wds >> wds >> wds;
			for (int i = 3; i <= 22; i++) {
				is >> wds;
				mtx[num][i - 2] = atof(wds.c_str());
			}
			num++;
		}
	}
	in.close();
}


bool GetOpts(int argc, char *argv[], OPTS &opts) {
	char tmp;
	while ((tmp = getopt(argc, argv, "hT:L:q:o:n:c:")) != -1) {
		switch (tmp) {
		case 'h': /* help */
			printf("!!! HELP !!!\n");
			return false;
			break;
		case 'T': /* library directory */
			opts.libdir = std::string(optarg);
			break;
		case 'L': /* list file of template IDs */
			opts.list = std::string(optarg);
			break;
		case 'q': /* name of the query */
			opts.target = std::string(optarg);
			break;
		case 'o': /* output directory */
                        opts.tardir = std::string(optarg);
                        break;
		case 'n': /* number of TOP templates*/
			opts.nTop = atoi(optarg);
			break;
		case 'c': /* seq identity cutoff */
                        opts.idcut = atof(optarg);
                        break;
		default:
			return false;
			break;
		}
}

if (opts.libdir == "") {
	printf("Error: Template library directory not specified ('-T')\n");
	return false;
}
if (opts.list == "") {
        printf("Error: template list not specified ('-L')\n");
        return false;
}
if (opts.target == ""){
	printf("Error: query not specified ('-q')\n");
        return false;
}
if (opts.tardir == "") {
        printf("Error: Query directory not specified ('-o')\n");
        return false;
}
if (opts.nTop > 50 || opts.nTop <=0){
	printf("Error: Number of TOP templates should be in a range [1,50] ('-n')\n");
        return false;
}
if (opts.idcut > 1 || opts.idcut <=0){
        printf("Error: sequence identity cutoff should be in a range [0,1] ('-c')\n");
        return false;
}
return true;
}

void PrintOpts(const OPTS &opts) {
	printf("\n-----------------------------------------------------------------------------------\n");
	printf(" DisCovER : Distance- and orientation-based Covariational threadER (%4s)\n",VERSION);
	printf(" For comments, please email to bhattacharyad@auburn.edu\n");
	printf("-----------------------------------------------------------------------------------\n");
	printf("\nUsage:   ./DisCovER [-option] [argument]\n\n");
	printf("Options:  -T path to template library             	- input, required\n");
	printf("          -L .txt file with template IDs 	 	- input, required\n");
	printf("          -q query ID   		          	- input, required\n");
	printf("          -o path to output directory              	- input, required\n");
	printf("          -n number of TOP templates to be selected  	(optional, default 50)\n");
	printf("          -c sequence identity cutoff	             	(optional, default 0.30)\n");
	printf("\n");
}

void PrintCap(const OPTS &opts) {
	time_t timer;
	time(&timer);
	struct tm* tm_info = localtime(&timer);
	char buf[100];
	strftime(buf, 26, "%Y:%m:%d / %H:%M:%S", tm_info);

	printf("# %s\n", std::string(80, '-').c_str());
	printf("# DisCovER : Distance- and orientation-based Covariational threadER (%4s)\n",VERSION);
	//printf("#\tDistance-based Covariational ThreadER (%4s)\n",VERSION);
	printf("# For comments, please email to bhattacharyad@auburn.edu\n");
	printf("# %s\n", std::string(80, '-').c_str());
	printf("# %20s : %s\n", "start date/time", buf);
	printf("# %s\n", std::string(80, '-').c_str());
}

//BLOSUM62 
void Blosum62Matrix(int imut[][24]){
                imut[1][1]=4;                                   // b,z,x are additional
                imut[1][2]=-1;
                imut[1][3]=-2;
                imut[1][4]=-2;
                imut[1][5]=0;
                imut[1][6]=-1;
                imut[1][7]=-1;
                imut[1][8]=0;
                imut[1][9]=-2;
                imut[1][10]=-1;
                imut[1][11]=-1;
                imut[1][12]=-1;
                imut[1][13]=-1;
                imut[1][14]=-2;
                imut[1][15]=-1;
                imut[1][16]=1;
                imut[1][17]=0;
                imut[1][18]=-3;
                imut[1][19]=-2;
                imut[1][20]=0;
                imut[1][21]=-2;
                imut[1][22]=-1;
                imut[1][23]=0;
                imut[2][1]=-1;
                imut[2][2]=5;
                imut[2][3]=0;
                imut[2][4]=-2;
                imut[2][5]=-3;
                imut[2][6]=1;
                imut[2][7]=0;
                imut[2][8]=-2;
                imut[2][9]=0;
                imut[2][10]=-3;
                imut[2][11]=-2;
                imut[2][12]=2;
                imut[2][13]=-1;
                imut[2][14]=-3;
                imut[2][15]=-2;
                imut[2][16]=-1;
                imut[2][17]=-1;
                imut[2][18]=-3;
                imut[2][19]=-2;
                imut[2][20]=-3;
                imut[2][21]=-1;
                imut[2][22]=0;
                imut[2][23]=-1;
                imut[3][1]=-2;
                imut[3][2]=0;
                imut[3][3]=6;
                imut[3][4]=1;
                imut[3][5]=-3;
                imut[3][6]=0;
                imut[3][7]=0;
                imut[3][8]=0;
                imut[3][9]=1;
                imut[3][10]=-3;
                imut[3][11]=-3;
                imut[3][12]=0;
                imut[3][13]=-2;
                imut[3][14]=-3;
                imut[3][15]=-2;
                imut[3][16]=1;
                imut[3][17]=0;
                imut[3][18]=-4;
                imut[3][19]=-2;
                imut[3][20]=-3;
		imut[3][21]=3;
                imut[3][22]=0;
                imut[3][23]=-1;
                imut[4][1]=-2;
                imut[4][2]=-2;
                imut[4][3]=1;
                imut[4][4]=6;
                imut[4][5]=-3;
                imut[4][6]=0;
                imut[4][7]=2;
                imut[4][8]=-1;
                imut[4][9]=-1;
                imut[4][10]=-3;
                imut[4][11]=-4;
                imut[4][12]=-1;
                imut[4][13]=-3;
                imut[4][14]=-3;
                imut[4][15]=-1;
                imut[4][16]=0;
                imut[4][17]=-1;
                imut[4][18]=-4;
                imut[4][19]=-3;
                imut[4][20]=-3;
                imut[4][21]=4;
                imut[4][22]=1;
                imut[4][23]=-1;
                imut[5][1]=0;
                imut[5][2]=-3;
                imut[5][3]=-3;
                imut[5][4]=-3;
                imut[5][5]=9;
                imut[5][6]=-3;
                imut[5][7]=-4;
                imut[5][8]=-3;
                imut[5][9]=-3;
                imut[5][10]=-1;
                imut[5][11]=-1;
                imut[5][12]=-3;
                imut[5][13]=-1;
                imut[5][14]=-2;
                imut[5][15]=-3;
                imut[5][16]=-1;
                imut[5][17]=-1;
                imut[5][18]=-2;
                imut[5][19]=-2;
                imut[5][20]=-1;
                imut[5][21]=-3;
                imut[5][22]=-3;
                imut[5][23]=-2;
                imut[6][1]=-1;
                imut[6][2]=1;
                imut[6][3]=0;
                imut[6][4]=0;
                imut[6][5]=-3;
                imut[6][6]=5;
                imut[6][7]=2;
                imut[6][8]=-2;
                imut[6][9]=0;
                imut[6][10]=-3;
                imut[6][11]=-2;
                imut[6][12]=1;
                imut[6][13]=0;
                imut[6][14]=-3;
                imut[6][15]=-1;
                imut[6][16]=0;
                imut[6][17]=-1;
                imut[6][18]=-2;
                imut[6][19]=-1;
		imut[6][20]=-2;
                imut[6][21]=0;
                imut[6][22]=3;
                imut[6][23]=-1;
                imut[7][1]=-1;
                imut[7][2]=0;
                imut[7][3]=0;
                imut[7][4]=2;
                imut[7][5]=-4;
                imut[7][6]=2;
                imut[7][7]=5;
                imut[7][8]=-2;
                imut[7][9]=0;
                imut[7][10]=-3;
                imut[7][11]=-3;
                imut[7][12]=1;
                imut[7][13]=-2;
                imut[7][14]=-3;
                imut[7][15]=-1;
                imut[7][16]=0;
                imut[7][17]=-1;
                imut[7][18]=-3;
                imut[7][19]=-2;
                imut[7][20]=-2;
                imut[7][21]=1;
                imut[7][22]=4;
                imut[7][23]=-1;
                imut[8][1]=0;
                imut[8][2]=-2;
                imut[8][3]=0;
                imut[8][4]=-1;
                imut[8][5]=-3;
                imut[8][6]=-2;
                imut[8][7]=-2;
                imut[8][8]=6;
                imut[8][9]=-2;
                imut[8][10]=-4;
                imut[8][11]=-4;
                imut[8][12]=-2;
                imut[8][13]=-3;
                imut[8][14]=-3;
                imut[8][15]=-2;
                imut[8][16]=0;
                imut[8][17]=-2;
                imut[8][18]=-2;
                imut[8][19]=-3;
                imut[8][20]=-3;
                imut[8][21]=-1;
                imut[8][22]=-2;
                imut[8][23]=-1;
                imut[9][1]=-2;
                imut[9][2]=0;
                imut[9][3]=1;
                imut[9][4]=-1;
                imut[9][5]=-3;
                imut[9][6]=0;
                imut[9][7]=0;
                imut[9][8]=-2;
                imut[9][9]=8;
                imut[9][10]=-3;
                imut[9][11]=-3;
                imut[9][12]=-1;
                imut[9][13]=-2;
                imut[9][14]=-1;
                imut[9][15]=-2;
                imut[9][16]=-1;
                imut[9][17]=-2;
                imut[9][18]=-2;
                imut[9][19]=2;
                imut[9][20]=-3;
		imut[9][21]=0;
                imut[9][22]=0;
                imut[9][23]=-1;
                imut[10][1]=-1;
                imut[10][2]=-3;
                imut[10][3]=-3;
                imut[10][4]=-3;
                imut[10][5]=-1;
                imut[10][6]=-3;
                imut[10][7]=-3;
                imut[10][8]=-4;
                imut[10][9]=-3;
                imut[10][10]=4;
                imut[10][11]=2;
                imut[10][12]=-3;
                imut[10][13]=1;
                imut[10][14]=0;
                imut[10][15]=-3;
                imut[10][16]=-2;
                imut[10][17]=-1;
                imut[10][18]=-3;
                imut[10][19]=-1;
                imut[10][20]=3;
                imut[10][21]=-3;
                imut[10][22]=-3;
                imut[10][23]=-1;
                imut[11][1]=-1;
                imut[11][2]=-2;
                imut[11][3]=-3;
                imut[11][4]=-4;
                imut[11][5]=-1;
                imut[11][6]=-2;
                imut[11][7]=-3;
                imut[11][8]=-4;
                imut[11][9]=-3;
                imut[11][10]=2;
                imut[11][11]=4;
                imut[11][12]=-2;
                imut[11][13]=2;
                imut[11][14]=0;
                imut[11][15]=-3;
                imut[11][16]=-2;
                imut[11][17]=-1;
                imut[11][18]=-2;
                imut[11][19]=-1;
                imut[11][20]=1;
                imut[11][21]=-4;
                imut[11][22]=-3;
                imut[11][23]=-1;
                imut[12][1]=-1;
                imut[12][2]=2;
                imut[12][3]=0;
                imut[12][4]=-1;
                imut[12][5]=-3;
                imut[12][6]=1;
                imut[12][7]=1;
                imut[12][8]=-2;
                imut[12][9]=-1;
                imut[12][10]=-3;
                imut[12][11]=-2;
                imut[12][12]=5;
                imut[12][13]=-1;
                imut[12][14]=-3;
                imut[12][15]=-1;
                imut[12][16]=0;
                imut[12][17]=-1;
                imut[12][18]=-3;
                imut[12][19]=-2;
                imut[12][20]=-2;
		imut[12][21]=0;
                imut[12][22]=1;
                imut[12][23]=-1;
                imut[13][1]=-1;
                imut[13][2]=-1;
                imut[13][3]=-2;
                imut[13][4]=-3;
                imut[13][5]=-1;
                imut[13][6]=0;
                imut[13][7]=-2;
                imut[13][8]=-3;
                imut[13][9]=-2;
                imut[13][10]=1;
                imut[13][11]=2;
                imut[13][12]=-1;
                imut[13][13]=5;
                imut[13][14]=0;
                imut[13][15]=-2;
                imut[13][16]=-1;
                imut[13][17]=-1;
                imut[13][18]=-1;
                imut[13][19]=-1;
                imut[13][20]=1;
                imut[13][21]=-3;
                imut[13][22]=-1;
                imut[13][23]=-1;
                imut[14][1]=-2;
                imut[14][2]=-3;
                imut[14][3]=-3;
                imut[14][4]=-3;
                imut[14][5]=-2;
                imut[14][6]=-3;
                imut[14][7]=-3;
                imut[14][8]=-3;
                imut[14][9]=-1;
                imut[14][10]=0;
                imut[14][11]=0;
                imut[14][12]=-3;
                imut[14][13]=0;
                imut[14][14]=6;
                imut[14][15]=-4;
                imut[14][16]=-2;
                imut[14][17]=-2;
                imut[14][18]=1;
                imut[14][19]=3;
                imut[14][20]=-1;
                imut[14][21]=-3;
                imut[14][22]=-3;
                imut[14][23]=-1;
                imut[15][1]=-1;
                imut[15][2]=-2;
                imut[15][3]=-2;
                imut[15][4]=-1;
                imut[15][5]=-3;
                imut[15][6]=-1;
                imut[15][7]=-1;
                imut[15][8]=-2;
                imut[15][9]=-2;
                imut[15][10]=-3;
                imut[15][11]=-3;
                imut[15][12]=-1;
                imut[15][13]=-2;
                imut[15][14]=-4;
                imut[15][15]=7;
                imut[15][16]=-1;
                imut[15][17]=-1;
                imut[15][18]=-4;
		imut[15][19]=-3;
                imut[15][20]=-2;
                imut[15][21]=-2;
                imut[15][22]=-1;
                imut[15][23]=-2;
                imut[16][1]=1;
                imut[16][2]=-1;
                imut[16][3]=1;
                imut[16][4]=0;
                imut[16][5]=-1;
                imut[16][6]=0;
                imut[16][7]=0;
                imut[16][8]=0;
                imut[16][9]=-1;
                imut[16][10]=-2;
                imut[16][11]=-2;
                imut[16][12]=0;
                imut[16][13]=-1;
                imut[16][14]=-2;
                imut[16][15]=-1;
                imut[16][16]=4;
                imut[16][17]=1;
                imut[16][18]=-3;
                imut[16][19]=-2;
                imut[16][20]=-2;
                imut[16][21]=0;
                imut[16][22]=0;
                imut[16][23]=0;
                imut[17][1]=0;
                imut[17][2]=-1;
                imut[17][3]=0;
                imut[17][4]=-1;
                imut[17][5]=-1;
                imut[17][6]=-1;
                imut[17][7]=-1;
                imut[17][8]=-2;
                imut[17][9]=-2;
                imut[17][10]=-1;
                imut[17][11]=-1;
                imut[17][12]=-1;
                imut[17][13]=-1;
                imut[17][14]=-2;
                imut[17][15]=-1;
                imut[17][16]=1;
                imut[17][17]=5;
                imut[17][18]=-2;
                imut[17][19]=-2;
                imut[17][20]=0;
                imut[17][21]=-1;
                imut[17][22]=-1;
                imut[17][23]=0;
                imut[18][1]=-3;
                imut[18][2]=-3;
                imut[18][3]=-4;
                imut[18][4]=-4;
                imut[18][5]=-2;
                imut[18][6]=-2;
                imut[18][7]=-3;
                imut[18][8]=-2;
                imut[18][9]=-2;
                imut[18][10]=-3;
                imut[18][11]=-2;
                imut[18][12]=-3;
                imut[18][13]=-1;
                imut[18][14]=1;
                imut[18][15]=-4;
                imut[18][16]=-3;
                imut[18][17]=-2;
                imut[18][18]=11;
                imut[18][19]=2;
		imut[18][20]=-3;
                imut[18][21]=-4;
                imut[18][22]=-3;
                imut[18][23]=-2;
                imut[19][1]=-2;
                imut[19][2]=-2;
                imut[19][3]=-2;
                imut[19][4]=-3;
                imut[19][5]=-2;
                imut[19][6]=-1;
                imut[19][7]=-2;
                imut[19][8]=-3;
                imut[19][9]=2;
                imut[19][10]=-1;
                imut[19][11]=-1;
                imut[19][12]=-2;
                imut[19][13]=-1;
                imut[19][14]=3;
                imut[19][15]=-3;
                imut[19][16]=-2;
                imut[19][17]=-2;
                imut[19][18]=2;
                imut[19][19]=7;
                imut[19][20]=-1;
                imut[19][21]=-3;
                imut[19][22]=-2;
                imut[19][23]=-1;
                imut[20][1]=0;
                imut[20][2]=-3;
                imut[20][3]=-3;
                imut[20][4]=-3;
                imut[20][5]=-1;
                imut[20][6]=-2;
                imut[20][7]=-2;
                imut[20][8]=-3;
                imut[20][9]=-3;
                imut[20][10]=3;
                imut[20][11]=1;
                imut[20][12]=-2;
                imut[20][13]=1;
                imut[20][14]=-1;
                imut[20][15]=-2;
                imut[20][16]=-2;
                imut[20][17]=0;
                imut[20][18]=-3;
                imut[20][19]=-1;
                imut[20][20]=4;
                imut[20][21]=-3;
                imut[20][22]=-2;
                imut[20][23]=-1;
                imut[21][1]=-2;
                imut[21][2]=-1;
                imut[21][3]=3;
                imut[21][4]=4;
                imut[21][5]=-3;
                imut[21][6]=0;
                imut[21][7]=1;
                imut[21][8]=-1;
                imut[21][9]=0;
                imut[21][10]=-3;
                imut[21][11]=-4;
                imut[21][12]=0;
                imut[21][13]=-3;
                imut[21][14]=-3;
                imut[21][15]=-2;
                imut[21][16]=0;
                imut[21][17]=-1;
                imut[21][18]=-4;
                imut[21][19]=-3;
                imut[21][20]=-3;
		imut[21][21]=4;
                imut[21][22]=1;
                imut[21][23]=-1;
                imut[22][1]=-1;
                imut[22][2]=0;
                imut[22][3]=0;
                imut[22][4]=1;
                imut[22][5]=-3;
                imut[22][6]=3;
                imut[22][7]=4;
                imut[22][8]=-2;
                imut[22][9]=0;
                imut[22][10]=-3;
                imut[22][11]=-3;
                imut[22][12]=1;
                imut[22][13]=-1;
                imut[22][14]=-3;
                imut[22][15]=-1;
                imut[22][16]=0;
                imut[22][17]=-1;
                imut[22][18]=-3;
                imut[22][19]=-2;
                imut[22][20]=-2;
                imut[22][21]=1;
                imut[22][22]=4;
                imut[22][23]=-1;
                imut[23][1]=0;
                imut[23][2]=-1;
                imut[23][3]=-1;
                imut[23][4]=-1;
                imut[23][5]=-2;
                imut[23][6]=-1;
                imut[23][7]=-1;
                imut[23][8]=-1;
                imut[23][9]=-1;
                imut[23][10]=-1;
                imut[23][11]=-1;
                imut[23][12]=-1;
                imut[23][13]=-1;
                imut[23][14]=-1;
                imut[23][15]=-2;
                imut[23][16]=0;
                imut[23][17]=0;
                imut[23][18]=-2;
                imut[23][19]=-1;
                imut[23][20]=-1;
                imut[23][21]=-1;
                imut[23][22]=-1;
                imut[23][23]=-1;
}
