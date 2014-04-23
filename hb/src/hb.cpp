#include "mcce.h"

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

const string PDBF = "step2_out.pdb"; //!< step2_out.pdb from MCCE.
const string HBOUT = "hb.dat";		//!< a binary file to save a n*n matrix of info about hbond between two conformers.
const string HBTXT = "hah.txt";		//!< a text file to store the hbond connection between conformers.
const string RESHBOND = "reshbond.txt";	//!< each line represents an edge of hb network of residues.
const string RES_IN_HBNET = "resInHbNet.txt";	//!< all the residues in hb network.

const float DFAR = 3.2;       //!< the distance used to determine a hydrogen bond
const float DNEAR = 1.2;       // 1.2 < d < 3.2
const float PI = 3.1415926;

bool is_hb(const CONF *conf1, const CONF *conf2, ofstream &h_fp);
int load_head3lst(RES &conflist);


int main()
{
	// Load all the tpl files from the path specified in run.prm.
	db_open();
	if (init()) {
		db_close();
		printf("Help message: double check file \"run.prm\" in current directory.\n");
		return USERERR;
	}

	// Load protein from a pdb file.
	FILE *fp;
	PROT prot;
	if (!(fp=fopen(STEP2_OUT, "r"))) {
		printf("   \"No step 2 output \"%s\".\n", STEP2_OUT);
		return USERERR;
	}
	prot = load_pdb(fp);
	if (prot.n_res == 0) {
		printf("   There are errors in pdb file, quiting ...\n");
		return USERERR;
	}

	// get the conformer uniq Id and atom connections
	id_conf(prot);
	get_connect12(prot);


	// load head3.lst to check with the conformers loaded from step2_out.pdb
	RES conflist;
	load_head3lst(conflist);
	int n_conf = conflist.n_conf;

	for (int i=0; i<prot.n_res; i++) {
		for (int j=1; j<prot.res[i].n_conf; j++) {
			strncpy(prot.res[i].conf[j].confName, prot.res[i].conf[j].uniqID, 5);
			prot.res[i].conf[j].confName[5] = '\0';
			// adjust the conf id to avoid id diff between step2 and head3.lst
			for (int k=0; k<n_conf; k++) {
				if (!strcmp(conflist.conf[k].uniqID, prot.res[i].conf[j].uniqID)) {
					prot.res[i].conf[j].iConf = conflist.conf[k].iConf;
					break;
				}
			}
		}
	}


	ofstream h_fp(HBTXT, ios::out);

	// hb_fp is a binary file to store the hbpw matrix, the first 4 bytes give the conformer number n_conf
	ofstream hb_fp(HBOUT, ios::out|ios::binary);
	hb_fp.write((char *)&n_conf, sizeof(int));

	// hbpw is a matrix to store the hb connection between each two conf, the elem is 0 or 1
	ofstream fp_res(RESHBOND, ios::out);
	char *resHb = (char *) calloc(prot.n_res * prot.n_res, sizeof(char));

	char *hbpw = (char *) calloc(n_conf * n_conf, sizeof(char));
	for (int i_res=0; i_res<prot.n_res; i_res++) {
		for (int i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
			for (int j_res=0; j_res<prot.n_res; j_res++) {
				if (j_res == i_res) continue;
				for (int j_conf=1; j_conf<prot.res[j_res].n_conf; j_conf++) {

					/* there is hbond from donor i_conf to accepter j_conf */
					if (is_hb(&prot.res[i_res].conf[i_conf], &prot.res[j_res].conf[j_conf], h_fp)) {
						hbpw[prot.res[i_res].conf[i_conf].iConf * n_conf + prot.res[j_res].conf[j_conf].iConf] = 1;
						resHb[i_res * prot.n_res + j_res] = 1;
					}
				}
			}
		}
	}

	hb_fp.write(hbpw, n_conf*n_conf);
	hb_fp.close();


	set <int> resInHbNet;
	for (int i_res=0; i_res<prot.n_res; i_res++) {
		for (int j_res=0; j_res<prot.n_res; j_res++) {
			if (j_res == i_res) continue;
			if (resHb[i_res * prot.n_res + j_res] == 1) {
				fp_res << prot.res[i_res].resName << prot.res[i_res].chainID;
				fp_res.fill('0');
				fp_res.width(4);
				fp_res << prot.res[i_res].resSeq;
				fp_res << '\t' << prot.res[j_res].resName << prot.res[j_res].chainID;
				fp_res.fill('0');
				fp_res.width(4);
				fp_res << prot.res[j_res].resSeq << endl;

				resInHbNet.insert(i_res);
				resInHbNet.insert(j_res);
			}
		}
	}
	fp_res.close();

	ofstream fp_resInHbNet(RES_IN_HBNET, ios::out);
	for (int i_res=0; i_res<prot.n_res; i_res++) {
		if (resInHbNet.find(i_res) != resInHbNet.end()) {
			fp_resInHbNet << prot.res[i_res].resName << prot.res[i_res].chainID;
			fp_resInHbNet.fill('0');
			fp_resInHbNet.width(4);
			fp_resInHbNet << prot.res[i_res].resSeq << endl;
		}
	}
	fp_resInHbNet.close();

	db_close();

	cout << " Number of residues loaded: " << prot.n_res << endl;



	return 0;
}


/** \brief Determine if there is hydrogen bond between two conformers.
 *
 * @param conf1 first conformer.
 * @param conf2 the second conformer.
 * @param h_fp an out file stream to write the actual hydrogen bond connection in.
 *
 * @return 0 if there is no hydrogen bond, 1 if there is.
 */
bool is_hb(const CONF *conf1, const CONF *conf2, ofstream &h_fp)
{
	STRINGS Datoms, Aatoms;
	int iD, iA, Dseq, Aseq;
	float d=0.0;
	float ang = 0.0;
	bool isHb = false;
	/** To establish hydrogen bond between two conformers, the following criteria must be met:
	 * 1. one conformer is h donor, the other is acceptor.
	 * 2. the distance between one donor atom (H) and one acceptor atom (O, N) is between DNEAR and DFAR.
	 * 3. the angle of h bond should be no less than 90.
	 */
	// the hydrogen donor in hb.tpl file is a hydrogen atom.
	if (!param_get("HDONOR", conf1->confName, "", &Datoms)) {
		if (!param_get("HACCEPT", conf2->confName, "", &Aatoms)) {
			for (iD=0; iD<Datoms.n; iD++) {
				param_get("IATOM", conf1->confName, Datoms.strings[iD], &Dseq);

				for (iA=0; iA<Aatoms.n; iA++) {
					param_get("IATOM", conf2->confName, Aatoms.strings[iA], &Aseq);
					// use the distance between the proton and the acceptor atom.
					d = ddvv(conf1->atom[Dseq].xyz, conf2->atom[Aseq].xyz);

					// use the donor heavy atom and the acceptor atom.
//					d = ddvv(conf1->atom[Dseq].connect12[0]->xyz, conf2->atom[Aseq].xyz);
					if (d > DNEAR * DNEAR && d < DFAR * DFAR) {
						ang = avv(vector_vminusv((conf1->atom[Dseq].connect12[0])->xyz, conf1->atom[Dseq].xyz),
								vector_vminusv(conf2->atom[Aseq].xyz, conf1->atom[Dseq].xyz)) * 180.0 / PI;
						if (abs(ang) < 90.0) continue;
						h_fp << conf1->uniqID << '\t' << conf2->uniqID << '\t' << conf1->atom[Dseq].connect12[0]->name
								<< '~' << conf1->atom[Dseq].name << "--" << conf2->atom[Aseq].name << '\t';
						h_fp.precision(2);
						h_fp << fixed << sqrt(d) << '\t';
						h_fp.precision(0);
						h_fp << fixed << ang << endl;

						isHb = true;
					}
				}
			}
		}
	}

	return isHb;
}


/** \brief Load conformers from head3.lst.
 *
 * To compare the conformers in head3.lst and step2_out.pdb.
 * They should be consistent.
 *
 * @param conflist a list of conformers loaded from head3.lst
 *
 * @return 0 if succeed
 */
int load_head3lst(RES &conflist)
{
	FILE *fp;
	char sbuff[MAXCHAR_LINE];
	char stemp[MAXCHAR_LINE];
	CONF conf_temp;
	int iconf;
	int counter;

	conflist.n_conf = 0;
	conflist.conf   = NULL;

	if (!(fp=fopen(FN_CONFLIST3, "r"))) {
		printf("   FATAL: Can't open file %s\n", FN_CONFLIST3);
		return USERERR;
	}
	fgets(sbuff, sizeof(sbuff), fp); /* skip the first line */
	counter = 0;
	while(fgets(sbuff, sizeof(sbuff), fp)) {
		/* load this line to a conf template */
		if (strlen(sbuff) < 20) continue;
		sscanf(sbuff, "%d %s %c %f %f %f %f %d %d %f %f %f %f %f %f %s", &conf_temp.iConf,
				conf_temp.uniqID,
				&conf_temp.on,
				&conf_temp.occ,
				&conf_temp.netcrg,
				&conf_temp.Em,
				&conf_temp.pKa,
				&conf_temp.e,
				&conf_temp.H,
				&conf_temp.E_vdw0,
				&conf_temp.E_vdw1,
				&conf_temp.E_tors,
				&conf_temp.E_epol,
				&conf_temp.E_dsolv,
				&conf_temp.E_extra,
				conf_temp.history);

		conf_temp.E_TS = 0.0; /* initialize entropy effect at the time of loading conflist */

		/* rescale */
		conf_temp.E_vdw0 *= env.scale_vdw0;
		conf_temp.E_vdw1 *= env.scale_vdw1;
		conf_temp.E_epol *= env.scale_ele;
		conf_temp.E_tors *= env.scale_tor;
		conf_temp.E_dsolv*= env.scale_dsolv;

		strncpy(conf_temp.resName, conf_temp.uniqID, 3); conf_temp.resName[3] = '\0';
		strncpy(conf_temp.confName, conf_temp.uniqID, 5); conf_temp.confName[5] = '\0';
		conf_temp.chainID = conf_temp.uniqID[5];
		strncpy(stemp, conf_temp.uniqID+6, 4); stemp[4] = '\0';
		conf_temp.resSeq = atoi(stemp);
		conf_temp.iCode = conf_temp.uniqID[10];
		conf_temp.n_atom = 0;
		if (conf_temp.on == 't' || conf_temp.on == 'T') conf_temp.on = 't';
		else conf_temp.on = 'f';
		conf_temp.iConf = counter;
		/* creating conflist */
		iconf = ins_conf(&conflist, conflist.n_conf, 0);
		cpy_conf(&conflist.conf[iconf], &conf_temp);
		counter++;
	}

	return 0;
}
