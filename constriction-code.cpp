//#In the name of Allah the most merciful, the most compassionate, who can help people to be and to do the best
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
///
///                      Non-gray BTE-FVM-Heat Transfer Cell_Centered Triangular Cells
///     Constriction Adiabatic (Diffusive,Specular) and Constant Temperature Bc's ; PETSc Read data from the cohll
///     modify the intial condition of e0
///     read initial data from ResultTempCell.dat (if not exist then from inputtemp.dat) but with some problems.
/// $$$ 2017-06-20 $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#ifdef __linux__ // When run in linux, this part is needed.
    #include <petscksp.h>
#endif
#include "mdarray.h"
using namespace std;

int numb, numnode, numcell;
const double PI = 4.0 * atan(1.0);
const double WFACTOR = 2;

// find the label boundary which needs to be omitted 
bool find_omit(int x, const set<int>& omit) {
    return (omit.find(x) != omit.end());
}

// a helper function of union
int Find (int x, int *parent) {
    if (*(parent + x) != x)
       *(parent + x) = Find(*(parent + x), parent);
    return *(parent + x);
}

// union the same boundary which may have different label
void Union (int x, int y, int *parent) {
    int a, b;
    a = Find(x, parent);
    b = Find(y, parent);
    if (a < b)
        *(parent + b) = a;
    else *(parent + a) = b;
}

// get direction of vectors
void get_direction(int nftot, int ntheta, int nphi, double WFACTOR,
    CMDArray<double>& sweight, CMDArray<double>& ss, CMDArray<double>& weight) {
    double dtheta= 0.5 * PI / ntheta;
    double dphi= 0.5 * PI / nphi;
    CMDArray<double> theta(nftot), phi(nftot);
    theta(0) = 0.5 * dtheta;
    phi(0) = 0.5 * dphi;
    for (int nt = 1; nt < ntheta; nt++)
        theta(nt) = theta(nt - 1) + dtheta;
    for (int np = 1; np < 4 * nphi; np++)
        phi(np) = phi(np - 1) + dphi;
    for (int nt=0; nt < ntheta; nt++){
        for (int np=0; np < 4 * nphi; np++){
            int nf = np + nt * 4 * nphi;
            weight(nf) = WFACTOR * 2. * sin(theta(nt)) * sin(0.5 * dtheta) * dphi;
            sweight(nf, 0) = WFACTOR * sin(phi(np)) * sin(0.5 * dphi) * (dtheta - cos(2 * theta(nt)) * sin(dtheta));
            sweight(nf, 1) = WFACTOR * cos(phi(np)) * sin(0.5 * dphi) * (dtheta - cos(2 * theta(nt)) * sin(dtheta));
            sweight(nf, 2) = WFACTOR * 0.5 * dphi * sin(2 * theta(nt)) * sin(dtheta);
            ss(nf, 0) = sin(theta(nt)) * sin(phi(np));
            ss(nf, 1) = sin(theta(nt)) * cos(phi(np));
            ss(nf, 2) = cos(theta(nt));
        }
    }
}

// find the interpolation elements and calculate the center of cells
void find_interpolation_element (CMDArray<int>& elmnod, CMDArray<int>& nelemnode,
    CMDArray<double>& node_r_m1T,CMDArray<double>& Cc,
    CMDArray<double>& p, CMDArray<int>& t) {
    /// Finding center cell and area A of cell [..][2]
    for (int elem1=0; elem1 < numcell; elem1++){
        Cc(elem1, 0) = (p(t(elem1)) + p(t(elem1 + numcell)) + p(t(elem1 + 2 * numcell))) / 3.0;
        Cc(elem1, 1) = (p(t(elem1) + numnode) + p(t(elem1 + numcell) + numnode)+ p(t(elem1 + 2 * numcell) + numnode)) / 3.0;
        Cc(elem1, 2) = 0.5 * (p(t(elem1)) * (p(t(elem1 + numcell) + numnode)-p(t(elem1 + 2 * numcell) + numnode))
        + p(t(elem1 + numcell)) * (p(t(elem1 + 2 * numcell) + numnode)-p(t(elem1) + numnode))
        + p(t(elem1 + 2 * numcell)) *(p(t(elem1) + numnode)-p(t(elem1 + numcell) + numnode)));
    }
    // Here, I deleted the area_node[], as it's useless from the code point of view.
    for (int elemnode0 = 0; elemnode0 < numnode; elemnode0++){
        nelemnode(elemnode0) = 0;
        for (int elemnode1=0; elemnode1 < numcell; elemnode1++){
            if (elemnode0 == t(elemnode1) || elemnode0 == t(elemnode1 + numcell) || elemnode0 == t(elemnode1 + 2 * numcell)){
                elmnod(elemnode0, nelemnode(elemnode0)) = elemnode1;
                nelemnode(elemnode0)++;
                node_r_m1T(elemnode0) += 1.0 / sqrt(pow((Cc(elemnode1, 0) - p(elemnode0)), 2)
                + pow((Cc(elemnode1, 1) - p(elemnode0 + numnode)), 2));
            }
        }
    }
}

// get the order number of boundary elements
void get_boundary_element(CMDArray<int>& elemboundary,
    CMDArray<int>& eboundary,CMDArray<int>& t, int numbDr) {
    for (int nodebndry=0; nodebndry < numbDr; nodebndry++){
        for (int cell1=0; cell1 < numcell; cell1++){
            if (eboundary(nodebndry) == t(cell1) ||
                eboundary(nodebndry) == t(cell1 + numcell) ||
                eboundary(nodebndry) == t(cell1 + 2 * numcell)){
                if (eboundary(nodebndry + numbDr) == t(cell1) ||
                    eboundary(nodebndry + numbDr) == t(cell1 + numcell) ||
                    eboundary(nodebndry + numbDr) == t(cell1 + 2 * numcell)){
                    elemboundary(nodebndry) = cell1;
                }
            }
        }
    }
}

// initialize e0
void initialize_e0 (CMDArray<double>& ee_n, CMDArray<double>& e0_n, CMDArray<double>& Temp,
		    CMDArray<double>& Temp_n, CMDArray<double>& Cc, CMDArray<double>& weight,
		    double* C_n, double Tref, double L_n, int nband, int nftot, bool flag) {
    // get in the data.
    ifstream fin_temp, fin_e0;
    //ofstream fout("test_ee.out"), fout1("test_e0.dat");
    double* node_temp;
    int num = numcell;
    //cout << num << endl;
    if (flag)
    	fin_temp.open("eeinputfilec.dat", ifstream::in);
    else fin_temp.open("eeinputfilew.dat", ifstream::in);
    if (fin_temp.good()) {
       //cout << "Read_ee" << endl;
       for (int i = 0; i < nband; i++)
	   for (int j = 0; j < nftot; j++) {
               for (int k = 0; k < numcell; k++) {
                   fin_temp >> setprecision(12) >> ee_n(k, j, i);
               	   //fout << setprecision(12) << ee_n(k, j, i) << ' ';
               }
               //fout << "\n";
           }
       if (flag)
    	  fin_e0.open("e0inputfilec.dat", ifstream::in);
       else fin_e0.open("e0inputfilew.dat", ifstream::in);
       for (int i = 0; i < nband; i++) {
	   for (int j = 0; j < numcell; j++) {
		fin_e0 >> setprecision(12) >> e0_n(j, i);
		//fout1 << setprecision(12) << e0_n(j, i) << ' ';
           }
	   //fout1 << "\n";
       }
       fin_e0.close();
    }
    else {
        //cout << "no ee" << endl;
        if (flag)
    	    fin_temp.open("ResultTempCell-constriction.dat", ifstream::in);
        else fin_temp.open("ResultTempCell-without.dat", ifstream::in);
    
        if (fin_temp.good()) {
	    node_temp = new double[num * 3];
	    int dump;
	    for (int i = 0; i < num; i++)
	        fin_temp >> dump >> node_temp[i] >> node_temp[i + num] >> node_temp[i + num * 2]; 
        }
        else {
            if (flag)
	    	fin_temp.open("inputtempc.dat", ifstream::in);
            else fin_temp.open("inputtempw.dat", ifstream::in);
            string str, num_str;
            stringstream ss;
            for (int i = 1; i <= 4; i++)
	        getline(fin_temp, str);
            getline(fin_temp, str);
            ss.str(str);
            for (int i = 1; i <= 3; i++)
	        ss >> num_str;
            num = atoi(num_str.c_str());
            for (int i = 1; i <= 4; i++)
	        getline(fin_temp, str);

            node_temp = new double[num * 3];
            for (int i = 0; i < num; i++) {
	        ss.clear();
	        getline(fin_temp, str);
	        ss.str(str);
	        for (int j = 0; j < 3; j++) {
	            ss >> num_str;
	            node_temp[j * num + i] = atof(num_str.c_str());
	        }
	        //cout << i << ' ' << node_temp[i] << ' ' << node_temp[num + i] << ' ' << node_temp[2 * num + i] << endl;
    	    }
        }
        //for (int i = 0; i < num; i++)
	     //cout << node_temp[i] << ' ' << node_temp[i + num] << ' ' << node_temp[i + num * 2] << endl;
        double sum = 0.0;
        for (int k = 0; k < nftot; k++)
	    sum += weight(k);
    
        for (int i = 0; i < numcell; i++){
	    double temp_min = 100.00, dist = 0.0;
	    int min_index = 0;
	    for (int j = 0; j < num; j++){
	        dist = sqrt((node_temp[j] - Cc(i, 0) / L_n) * (node_temp[j] - Cc(i, 0) / L_n)
			    + (node_temp[j + num] - Cc(i, 1) / L_n) * (node_temp[j + num] - Cc(i, 1) / L_n));
	        if(dist < temp_min){
		    temp_min = dist;
		    min_index = j;
	        }
	    }
	    Temp(i) = node_temp[min_index + num * 2];
	    //cout << i << ' ' << temp_min << ' ' << min_index << ' ' << Temp(i) << endl;
	    for (int j = 0; j < nband; j++) {
	        Temp_n(i, j) = Temp(i);
	        e0_n(i, j) = C_n[j] * (Temp(i) - Tref);
	        double temp_e0 = e0_n(i, j) / sum;
	        e0_n(i, j) /= (4 * PI);
	        for (int k = 0; k < nftot; k++)
		    ee_n(i, k, j) = temp_e0;
	    }
        }
        delete[] node_temp;
    }
    fin_temp.close();
    //cout << "AFTER DELETE PASS" << endl;
}

/// get the coefficient of matrix A and vector B in linear system Ax=b
void get_coefficient(const int nedge,int iter, int iband, int inf, int numbDr,
    double xe[], double* mfp_n, double* C_n, double Tleft, double Tright, double Tref,double L_n,
    int ntheta, int nphi, int ntop, int ntopi, int nbottom, int nbottomi, int nleft,
    int nlefti, int nright, int nrighti, int nrighti1, CMDArray<double>& p, CMDArray<int>& t,
    CMDArray<double>& sweight, CMDArray<double>& ss, CMDArray<double>& weight, CMDArray<double>& Cc,
    CMDArray<int>& eboundary, CMDArray<double>& e0_n, CMDArray<double>& e1_n, CMDArray<double>& Ke, CMDArray<double>& Re) {
    double a_f[nedge];
    int nel[3], nfmax = 4 * ntheta * nphi;
    for (int ie = 0; ie < numcell; ie++){
        nel[0] = t(ie);
        nel[1] = t(ie + numcell);
        nel[2] = t(ie + 2 * numcell);
        xe[0]  = p(nel[0]);
        xe[4]  = p(nel[0] + numnode);
        xe[1]  = p(nel[1]);
        xe[5]  = p(nel[1] + numnode);
        xe[2]  = p(nel[2]);
        xe[6]  = p(nel[2] + numnode);
        xe[3]  = xe[0];
        xe[7]  = xe[4];
        /// n1 n2 n3 (unit normal vector)
        double norm[2] = {0};
        double leng_edge;
        double swn;
        int neledge[3][2]={0}, nell[3];;
        neledge[0][0] = nel[0];
        neledge[0][1] = neledge[1][0] = t(ie + numcell);
        neledge[1][1] = neledge[2][0] = t(ie + 2 * numcell);
        neledge[2][1] = t(ie);
        for (int iedge = 0; iedge < nedge; iedge++){
            int ineighb = 0, eneighb = -1;
            /// which element is for this edge
            for (int ie1 = 0; ie1 < numcell; ie1++){
                nell[0] = t(ie1);
                nell[1] = t(ie1 + numcell);
                nell[2] = t(ie1 + 2 * numcell);
                if (ie1 != ie){
                    if (neledge[iedge][0] == nell[0] ||
                        neledge[iedge][0] == nell[1] ||
                        neledge[iedge][0] == nell[2]){
                        if (neledge[iedge][1] == nell[0] ||
                            neledge[iedge][1] == nell[1] ||
                            neledge[iedge][1] == nell[2]){
                            ineighb = 1;
                            eneighb = ie1;
                            break;
                        }
                    }
                    int ib;
                    for (ib = 0; ib < numbDr; ib++){
                        if ((neledge[iedge][0] == eboundary(ib) &&
                            neledge[iedge][1] == eboundary(ib + numb)) || (neledge[iedge][0] == eboundary(ib + numb) &&
                            neledge[iedge][1] == eboundary(ib))){
                            ineighb=0;
                            break;
                        }
                    }
                    if (ib < numbDr)
                    break;
                }
            }
            leng_edge = sqrt(pow(xe[iedge + 1] - xe[iedge], 2) + pow(xe[iedge + 1 + nedge + 1] - xe[iedge + 1 + nedge], 2));
            norm[0] = (xe[iedge + 1 + nedge + 1] - xe[iedge + 1 + nedge]) / leng_edge;
            norm[1] = -(xe[iedge + 1] - xe[iedge]) / leng_edge;

/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$/// find e0 and temperature
            if (ineighb){ //means neighbor cell exists and it is not boundary
                swn = (ss(inf, 0) * norm[0] + ss(inf, 1) * norm[1]);
                a_f[iedge] = mfp_n[iband] *
                (sweight(inf, 0) * norm[0] + sweight(inf, 1) * norm[1]) * leng_edge;
                if (swn>=0)
                    Ke(ie, ie) += a_f[iedge];
                else Ke(ie, eneighb) += a_f[iedge];
            }

/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
///                                                 Boundaries
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            else { //means neighbor cell not exists and it is boundary
                for (int ib = 0; ib < numbDr; ib++){
                    if (((neledge[iedge][0] == eboundary(ib) &&
                        neledge[iedge][1] == eboundary(ib + numb))) || ((neledge[iedge][0] == eboundary(ib + numb) &&
                        neledge[iedge][1] == eboundary(ib)))){
                        swn = (ss(inf, 0) * norm[0] + ss(inf, 1) * norm[1]);
                        ///! outgoing
                        a_f[iedge] = mfp_n[iband] * (sweight(inf, 0) * norm[0] + sweight(inf, 1) * norm[1]) * fabs(leng_edge);
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// Diffusive-Specular boundary                           Top
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// DIFFUSIVE BOUNDARY
                        if (eboundary(ib + numbDr * 2) == ntopi){
                            double einsum = 0.0;
                            for (int nft = 0; nft < nfmax; nft++){
                                swn = (ss(nft, 0) * norm[0] + ss(nft, 1) * norm[1]);
                                if (swn > 0) {
                                    einsum = einsum + e1_n(ie,nft,iband) * sweight(nft, 1);
                                }
                            }
                            einsum = einsum / PI;
                            for (int nft2 = 0; nft2 < nfmax; nft2++){
                                swn = (ss(nft2, 0) * norm[0] + ss(nft2, 1) * norm[1]);
                                if (inf == nft2) {
                                    if (swn >= 0){       ///    outgoing
                                        Ke(ie, ie) += a_f[iedge];
                                    }
                                    else if (swn < 0) {   ///   incoming
                                        Re(ie) -= einsum * mfp_n[iband] * (sweight(inf, 0) * norm[0] + sweight(inf, 1) * norm[1]) * leng_edge;
                                    }
                                }
                            }
                        }
/// SPECULAR BOUNDARY
                        else if (eboundary(ib + numbDr * 2) == ntop){
                            for (int nt = 0; nt < ntheta; nt++){
                                for (int np = 0; np < 4 * nphi; np++){
                                    int nf1 = np + nt * 4 * nphi;
                                    int npref, nfref;
                                    if (inf == nf1) {
                                        cout<<"";
                                        if (swn >= 0){       //         outgoing
                                            Ke(ie, ie) += a_f[iedge];
                                            cout<<"";
                                        }
                                        else if (swn < 0) {  //         incoming
                                            int np2 = np + 1;
                                            int nt2 = nt + 1;
                                            if (np2 > nphi && np2 <= 2 * nphi){
                                                npref = 2 * nphi + 1 - np2;
                                            }
                                            else {
                                                npref = 6 * nphi + 1 - np2;
                                            }
                                            nfref = npref + (nt2 - 1) * 4 * nphi - 1;
                                            Re(ie) -= e1_n(ie, nfref, iband) * mfp_n[iband] * (sweight(inf, 0) * norm[0] + sweight(inf, 1) * norm[1]) * leng_edge;
                                        }
                                    }
                                }
                            }
                        }
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// Diffusive-Specular boundary                        Bottom
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// DIFFUSIVE BOUNDARY
                        else if (eboundary(ib + numbDr * 2) == nbottomi){
                            double einsum = 0.0;
                            for (int nft = 0;nft < nfmax; nft++){
                                swn=(ss(nft, 0) * norm[0] + ss(nft, 1) * norm[1]);
                                if (swn > 0) {
                                    einsum = einsum - e1_n(ie,nft,iband) * sweight(nft, 1);
                                }
                            }
                            einsum = einsum / PI;
                            for (int nft2 = 0; nft2 < nfmax; nft2++){
                                swn=(ss(nft2, 0) * norm[0] + ss(nft2, 1) * norm[1]);
                                if (inf == nft2) {
                                    if (swn >= 0){         ///  outgoing
                                        Ke(ie, ie) += a_f[iedge];
                                    }
                                    else if (swn<0) {      ///  incoming
                                        Re(ie) -= einsum * mfp_n[iband] * (sweight(inf, 0) * norm[0] + sweight(inf, 1) * norm[1]) * leng_edge;
                                    }
                                }
                            }
                        }
/// SPECULAR BOUNDARY
                        else if (eboundary(ib + numbDr * 2) == nbottom){
                            for (int nt = 0; nt < ntheta; nt++){
                                for (int np = 0; np < 4 * nphi; np++){
                                    int nf1 = np + nt * 4 * nphi;
                                    int npref, nfref;
                                    if (inf == nf1) {
                                        cout<<"";
                                        if (swn >= 0){     ///  outgoing
                                            Ke(ie, ie) += a_f[iedge];
                                            cout<<"";
                                        }
                                        else if (swn < 0) {///  incoming
                                            int np2 = np + 1;
                                            int nt2 = nt + 1;
                                            if (np2 <= nphi){
                                                npref = 2 * nphi + 1 - np2;
                                            }
                                            else {
                                                npref = 6 * nphi + 1 - np2;
                                            }
                                            nfref = npref + (nt2 - 1) * 4 * nphi - 1;
                                            Re(ie) -= e1_n(ie, nfref, iband) * mfp_n[iband] * (sweight(inf, 0) * norm[0] + sweight(inf, 1) * norm[1]) * leng_edge;
                                        }
                                    }
                                }
                            }
                        }
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// Diffusive-Specular boundary                  Internal Top
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// DIFFUSIVE BOUNDARY
                        else if (eboundary(ib + numbDr * 2) == ntopi){
                            double einsum = 0.0;
                            for (int nft = 0; nft < nfmax; nft++){
                                swn = ss(nft, 0) * norm[0] + ss(nft, 1) * norm[1];
                                if (swn > 0) {
                                    einsum = einsum - e1_n(ie, nft, iband) * sweight(nft, 1);
                                }
                            }
                            einsum = einsum / PI;
                            for (int nft2 = 0; nft2 < nfmax; nft2++){
                                swn = ss(nft2, 0) * norm[0] + ss(nft2, 1) * norm[1];
                                if (inf == nft2) {
                                    if (swn >= 0){         ///  outgoing
                                        Ke(ie, ie) += a_f[iedge];
                                    }
                                    else if (swn < 0) {    ///  incoming
                                        Re(ie) -= einsum * mfp_n[iband] *
                                                (sweight(inf, 0) * norm[0] +
                                                        sweight(inf, 1) * norm[1]) * leng_edge;
                                    }
                                }
                            }
                        }
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// Diffusive-Specular boundary               Internal Bottom
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// DIFFUSIVE BOUNDARY
                        else if (eboundary(ib + numbDr * 2) == nbottomi){
                            double einsum = 0.0;
                            for (int nft=0; nft < nfmax; nft++){
                                swn = ss(nft, 0) * norm[0] + ss(nft, 1) * norm[1];
                                if (swn > 0) {
                                    einsum = einsum + e1_n(ie, nft, iband) * sweight(nft, 1);
                                }
                            }
                            einsum = einsum / PI;
                            for (int nft2 = 0; nft2 < nfmax; nft2++){
                                swn = ss(nft2, 0) * norm[0] + ss(nft2, 1) * norm[1];
                                if (inf == nft2) {
                                    if (swn >= 0){         ///  outgoing
                                        Ke(ie, ie) += a_f[iedge];
                                    }
                                    else if (swn < 0) {    ///  incoming
                                        Re(ie) -= einsum * mfp_n[iband] * (sweight(inf, 0) * norm[0] + sweight(inf, 1) * norm[1]) * leng_edge;
                                    }
                                }
                            }
                        }
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// Diffusive-Specular boundary                Internal Right
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// DIFFUSIVE BOUNDARY
                        else if (eboundary(ib + numbDr * 2) == 11) {
                            double einsum = 0.0;
                            for (int nft=0; nft < nfmax; nft++){
                                swn = ss(nft, 0) * norm[0] + ss(nft, 1) * norm[1];
                                if (swn >= 0) {
                                    einsum = einsum - e1_n(ie, nft, iband) * sweight(nft, 0);
                                }
                            }
                            einsum = einsum / PI;
                            for (int nft2 = 0; nft2 < nfmax; nft2++){
                                swn = ss(nft2, 0) * norm[0] + ss(nft2, 1) * norm[1];
                                if (inf == nft2) {
                                    if (swn >= 0){         ///  outgoing
                                        Ke(ie, ie) += a_f[iedge];
                                    }
                                    else if (swn < 0) {    ///  incoming
                                        Re(ie) -= einsum * mfp_n[iband] * (sweight(inf, 0) * norm[0] + sweight(inf, 1) * norm[1]) * leng_edge;
                                    }
                                }
                            }
                        }
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// Diffusive-Specular boundary                 Internal Left
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// DIFFUSIVE BOUNDARY
                        else if (eboundary(ib + numbDr * 2) == nlefti) {
                            cout<<"";
                            double einsum = 0.0;
                            for (int nft = 0; nft < nfmax; nft++){
                                swn = ss(nft, 0) * norm[0] + ss(nft, 1) * norm[1];
                                if (swn >= 0) {
                                    if (p(t(ie+2*numcell))>0.5*L_n){
                                       einsum = einsum - e1_n(ie, nft, iband) * sweight(nft, 0);
                                     }
                                    else{

                                        einsum = einsum + e1_n(ie, nft, iband) * sweight(nft, 0);
                                    }
                                }
                            }
                            einsum = einsum / PI;
                            for (int nft2 = 0; nft2 < nfmax; nft2++){
                                swn = ss(nft2, 0) * norm[0] + ss(nft2, 1) * norm[1];
                                if (inf == nft2) {
                                    if (swn >= 0){                      /// !         outgoing
                                        Ke(ie, ie) += a_f[iedge];
                                    }
                                    else if (swn < 0){                  /// !         incoming
                                        Re(ie) -= einsum * mfp_n[iband] * (sweight(inf, 0) * norm[0] + sweight(inf, 1) * norm[1]) * leng_edge;
                                        //double ceo,ceo2;
                                        //ceo=floor(iter/100.);
                                        //ceo2=iter/100.;
                                        //if ((ceo2==ceo)&& (iter>0)){
                                        //cout<<"nlefi-SWN00= "<<ib<<" "<<inf<<" "<<ie<<" t "<< t(ie)<<" "<< t(ie+numcell)<<" "<< t(ie+2*numcell)<<" eb: "<<eboundary(ib+numbDr*2)<<"\n";
                                        //cout<<"nlefi-SWN<0= "<<ib<<" "<<inf<<" "<<ie<<"                  swn "<< swn <<" "<<einsum<<" "<<Re(ie)<<"\n";
                                        //
                                        //}
                                   }
                                }
                            }
                        }
/// SPECULAR BOUNDARY
                        else if (eboundary(ib + numbDr * 2) == 111){
                            for (int nt = 0; nt < ntheta; nt++){
                                for (int np = 0; np < 4 * nphi; np++){
                                    int nf1 = np + nt * 4 * nphi;
                                    int npref, nfref;
                                    if (inf == nf1) {
                                        cout<<"";
                                        if (swn >= 0){     ///  outgoing
                                            Ke(ie, ie) += a_f[iedge];
                                            cout<<"";
                                        }
                                        else if (swn < 0) {///  incoming
                                            int np2 = np + 1;
                                            int nt2 = nt + 1;
                                            if (np2 <= nphi){
                                                npref = 2 * nphi + 1 - np2;
                                            }
                                            else {
                                                npref = 6 * nphi + 1 - np2;
                                            }
                                            nfref = npref + (nt2 - 1) * 4 * nphi - 1;
                                            Re(ie) -= e1_n(ie, nfref, iband) * mfp_n[iband] * (sweight(inf, 0) * norm[0] + sweight(inf, 1) * norm[1]) * leng_edge;
                                        }
                                    }
                                }
                            }
                        }
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// Temperature Constant                         Left & Right
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

/// Incoming
                        else if (swn < 0) {
/// Boundary right
                            if (eboundary(ib + numbDr * 2) == nright || eboundary(ib + numbDr * 2) == nlefti * 2 ||
                                eboundary(ib + numbDr * 2) == 40 || eboundary(ib + numbDr * 2) == 50){
                                Re(ie) -= C_n[iband] / (4 * PI) * (Tright - Tref) * mfp_n[iband] * (sweight(inf, 0) * norm[0] + sweight(inf, 1) * norm[1]) * leng_edge;
                            }
/// Boundary left
                            else if (eboundary(ib + numbDr * 2) == nleft || eboundary(ib + numbDr * 2) == nbottomi){
                                Re(ie) -= C_n[iband] / (4 * PI) * (Tleft - Tref) * mfp_n[iband] * (sweight(inf, 0) * norm[0] + sweight(inf, 1) * norm[1]) * leng_edge;
                            }
                        }
                        else if (swn >= 0){
/// Outgoing
                            Ke(ie, ie) += a_f[iedge]  ;
                        }
                    }
                }
            }
        }//iedge
        Ke(ie, ie) += Cc(ie, 2) * weight(inf);
        Re(ie) += e0_n(ie,iband) * Cc(ie, 2) * weight(inf);
    }//ie
}

// PETSc part which can only be run in Linux.
#ifdef __linux__
vector<double> solve_matrix (CMDArray<double>& Ke, CMDArray<double>& Re) {
    //cout << "Linux" << endl;
    //ofstream fout("solution.out");
    int* nnz = new int[numcell], count;
    for (int i = 0; i < numcell; i++) {
        count = 0;
        for (int j = 0; j < numcell; j++)
            if (Ke(i, j) != 0)
                count++;
        *(nnz + i) = count;
    }
    //for (int i = 0; i < numcell; i++){
    	//for (int j = 0; j < numcell; j++)
             //fout << Ke(i, j) << ' ';
        //fout << Re(i) << "\n";
    //}
    static char help[] = "Solving matrix.\n\n";
    const double offset = 1e16;
    PetscInitialize(NULL, NULL, (char*)0, help);
    Mat	 A;        //matrix definition
    Vec xxx, bbb;   //vector definition
    PetscScalar  v, u, yyy[numcell]; //double in c++
    PetscInt iii, jjj, one, ix[numcell]; //integer in c++
    KSP ksp;   //solver of matrix
    MatCreateSeqAIJ(PETSC_COMM_SELF, numcell, numcell, 3, nnz, &A);
    //numm=numcell; 3 maximun non zero ; nnz : number of non zero in each row
    for (iii = 0;iii < numcell; iii++){
        for (jjj = 0; jjj < numcell; jjj++){
            if (Ke(iii, jjj) != 0) {
                v = Ke(iii, jjj) * offset; //make bigger digits
                MatSetValues(A, 1, &iii, 1, &jjj, &v, ADD_VALUES);
                // copy Ke to Petsc; 1 1 for 1 elements in matrix
            }
        }
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    //MatView(A,PETSC_VIEWER_STDOUT_WORLD);//print the matrix
    VecCreate(PETSC_COMM_WORLD, &bbb); // create bbb
    VecSetSizes(bbb, PETSC_DECIDE, numcell); //set the size for bbb (right hand side value)
    VecSetFromOptions(bbb);
    VecDuplicate(bbb, &xxx); // create another solution vector (xxx)
    for (iii = 0; iii < numcell; iii++){
        u = Re(iii) * offset;// make bigger digits
        VecSetValues(bbb, 1, &iii, &u, ADD_VALUES);//set value for bbb just one value in each time
    }
    //VecView(bbb,PETSC_VIEWER_STDOUT_WORLD);//print vector
///create solver of the system
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A); //format of petsc ;
    KSPSetFromOptions(ksp);

    KSPSolve(ksp,bbb,xxx); ///solve the matrix
    //VecView(xxx,PETSC_VIEWER_STDOUT_WORLD); //the results is xxx

    one = numcell;
    vector<double> x(numcell);

    for (iii = 0; iii < numcell; iii++)
        ix[iii]=iii;

    VecGetValues(xxx, one, ix, yyy);
    // send value from petsc to petsc another format ( petsc scaler : yyy) (xxx:vector . can not do operation on that)
    for (iii = 0; iii < numcell; iii++)
        x[iii]=yyy[iii];
    delete[] nnz;
    VecDestroy(&xxx);
    VecDestroy(&bbb);
    MatDestroy(&A);
    KSPDestroy(&ksp);
    return x;
/// finished petsc 
}
#endif

// Simple gauss method which used in Windows system.
#if (defined _WIN32) || (defined _WIN64)
vector<double> gauss(vector< vector<double> > Km_Re)
{
    int n = Km_Re.size();
    for (int i=0; i<n; i++)
    {
        // Search for maximum in this column
        double maxEl = fabs(Km_Re[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++)
        {
            if (fabs(Km_Re[k][i]) > maxEl)
            {
                maxEl = fabs(Km_Re[k][i]);
                maxRow = k;
            }
        }
        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1; k++)
        {
            double tmp = Km_Re[maxRow][k];
            Km_Re[maxRow][k] = Km_Re[i][k];
            Km_Re[i][k] = tmp;
        }
        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++)
        {
            double c = -Km_Re[k][i]/Km_Re[i][i];
            for (int j=i; j<n+1; j++)
            {
                if (i==j)
                {
                    Km_Re[k][j] = 0;
                }
                else
                {
                    Km_Re[k][j] += c * Km_Re[i][j];
                }
            }
        }
    }
    vector<double> x(n); // Solve equation Ax=b for an upper triangular matrix A
    for (int i=n-1; i>=0; i--)
    {
        x[i] = Km_Re[i][n]/Km_Re[i][i];
        for (int k=i-1; k>=0; k--)
        {
            Km_Re[k][n] -= Km_Re[k][i] * x[i];
        }
    }
    return x;
}

vector<double> solve_matrix (CMDArray<double>& Ke, CMDArray<double>& Re) {
    CMDArray<double> Km_active(numcell * numcell);
    CMDArray<double> Re_active(numcell);
    vector<double> line(numcell + 1,0);
    vector< vector<double> > Km_Re(numcell, line);
    for (int i1 = 0; i1 < numcell; i1++) {
        for (int i2 = 0; i2 < numcell; i2++)
            Km_Re[i1][i2] = Ke(i1, i2);
            Km_Re[i1][numcell] = Re(i1);
    }
    vector<double> x(numcell) ;
    x = gauss(Km_Re);
    return x;
}
#endif

// get the temperature of cells
void get_cell_temp (int iband, int nftot, CMDArray<double>& Cc, CMDArray<double>& ee_n,
		    CMDArray<double>& weight, CMDArray<double>& Temp_n, double* C_n, double Tref, ofstream& fout) {
    CMDArray<double> e0(numcell);
    e0.set_zero();
    int width_numcell = 0, num_temp = numcell;
    while (num_temp > 0) {
        num_temp /= 10;
	    width_numcell ++;
    }
    int double_width = 12;
    for (int ic = 0; ic < numcell; ic++){
        for (int iinf = 0; iinf < nftot; iinf++){
            e0(ic) += ee_n(ic, iinf, iband) * weight(iinf);
         }
        Temp_n(ic, iband) = e0(ic) / C_n[iband] + Tref;
        e0(ic) = e0(ic) / (4 * PI);
    }//ic
    for (int i2 = 0; i2 < numcell; i2++) {
	fout << setprecision(6) << setiosflags(ios::left);
        fout << setw(width_numcell) << i2 << ' ' << setw(double_width) << Cc(i2, 0) << ' '
	     << setw(double_width) << Cc(i2, 1) << ' ' << setw(7) << Temp_n(i2,iband) << "\n";
    }
}

// recover temperature from energy density
void recover_temp(CMDArray<double>& e0_n,
    CMDArray<double>& Temp_n, CMDArray<double>& Temp,
    int nband, double* R_n, double* C_n, double Tref) {
    double RRn, Rnn;
    for (int icc = 0; icc < numcell; icc++){
        RRn = Rnn = 0.0;
        for (int ibnd1 = 0; ibnd1 < nband; ibnd1++){
            RRn += (R_n[ibnd1] * Temp_n(icc,ibnd1));
            Rnn += R_n[ibnd1];
        }
        Temp(icc) = RRn / Rnn;
    }//icc
    for (int icc2 = 0; icc2 < numcell; icc2++)
        for (int iibnd = 0; iibnd < nband; iibnd++)
            e0_n(icc2,iibnd) = C_n[iibnd] / (4 * PI) * (Temp(icc2) - Tref);
}

// interpolation
void interpolation(int nband, CMDArray<double>& Temp, CMDArray<double>& Temnode_n,
    CMDArray<double>& node_r_m1T, CMDArray<int>& elmnod, CMDArray<int>& nelemnode,
    CMDArray<double>& Cc, CMDArray<double>& p,CMDArray<int>& t, ofstream& fout) {
    int width_numnode = 0, num_temp = numnode;
    while (num_temp > 0) {
    	num_temp /= 10;
	width_numnode ++;
    }
    int double_width = 12;
    /// Finding center cell and area A of cell [..][2]
    for (int elem1=0; elem1 < numcell; elem1++){
        Cc(elem1, 0) = (p(t(elem1)) + p(t(elem1 + numcell))+ p(t(elem1 + 2 * numcell))) / 3.0;
        Cc(elem1, 1) = (p(t(elem1) + numnode) + p(t(elem1 + numcell) + numnode)+ p(t(elem1 + 2 * numcell) + numnode)) / 3.0;
        Cc(elem1, 2) = 0.5 * (p(t(elem1)) * (p(t(elem1 + numcell) + numnode)-p(t(elem1 + 2 * numcell) + numnode))+ p(t(elem1 + numcell)) * (p(t(elem1 + 2 * numcell) + numnode)
        -p(t(elem1) + numnode))+ p(t(elem1 + 2 * numcell)) *(p(t(elem1) + numnode)-p(t(elem1 + numcell) + numnode)));
    }
    for (int iband = 0; iband < nband; iband++){
        for (int inode = 0; inode < numnode; inode++){
            Temnode_n(inode) = 0;
            for (int ino = 0; ino < nelemnode(inode); ino++){
                Temnode_n(inode) += Temp(elmnod(inode, ino))* (1.0 / (sqrt((pow((Cc(elmnod(inode, ino), 0)
                - p(inode)), 2)+ pow((Cc(elmnod(inode, ino), 1)-p(inode + numnode)), 2)))))/node_r_m1T(inode);
            }
 	    fout << setprecision(6) << setiosflags(ios::left);
            fout << setw(width_numnode) << inode << ' ' << setw(double_width) <<p(inode) << " "
		    << setw(double_width) << p(inode + numnode)
            << " " << setw(7) <<Temnode_n(inode) << "\n";
        }
    }
}

// calculate the heat transfer flux
void get_heat_transfer_flux(double* C_n,double Tleft,double Tright,double Tref,
    CMDArray<double>& sweight,CMDArray<int>& eboundary,CMDArray<int>& elemboundary,
    CMDArray<double>& ee_n, CMDArray<double>& weight,CMDArray<double>& ss, double* R_n, double* vg_n,
    int nband, int nftot, int numbDr, int ntop,int nbottom, int nleft, int nright, ofstream& fout,ofstream& fout2) {
    double *qleft = new double[numb];
    double *qleft2 = new double[numb];
    double *qright = new double[numb];
    //double *qrightsum = new double[numb];
    double *qright2 = new double[numb];
    double *qtop = new double[numb];
    double *qbot = new double[numb];
    //double *qleftsum = new double[numb];
    int width_numbDr = 1, num_temp = numbDr;
    while (num_temp > 0) {
    	num_temp /= 10;
	width_numbDr ++;
    }
    int width_numcell = 0; num_temp = numcell;
    while (num_temp > 0) {
    	num_temp /= 10;
	width_numcell ++;
    }
    int double_width = 15;

/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
///                                   Calculation of Heat Transfer flux
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    double qleftsum=0,  qrightsum=0;
    int counterr=0,counterl=0;
    fout <<setiosflags(ios::scientific)<<setiosflags(ios::left)<<setprecision(5);
    for (int iboun = 0; iboun < numbDr; iboun++){
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// Heat flux                                            Left
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if (eboundary(iboun + 2 * numbDr) == nleft){
            qleft[iboun] =qleft2[iboun] =  qright[iboun] = qtop[iboun] = qbot[iboun] = 0;
            for (int ibandc = 0; ibandc < nband; ibandc++){
                for (int infec = 0; infec < nftot; infec++){
                    qleft[iboun] += R_n[ibandc] * ee_n(elemboundary(iboun), infec, ibandc) * weight(infec) * ss(infec, 0) * vg_n[ibandc];
                    if (sweight(infec, 0)>0){
                        qleft2[iboun]=qleft2[iboun]+R_n[ibandc]*C_n[ibandc]/(4*PI)*(Tleft-Tref)*sweight(infec, 0)*vg_n[ibandc];
                    }
                    else if (sweight(infec, 0)>0){
                        qleft2[iboun]=qleft2[iboun] +R_n[ibandc]*ee_n(elemboundary(iboun),infec,ibandc)*sweight(infec, 0)*vg_n[ibandc];
                    }
                    if ((ibandc == nband - 1) && (infec == nftot - 1)){
                        qleftsum=qleftsum+qleft[iboun];
                        counterl=counterl+1;
                        fout<<"qleft:  "<<setw(width_numbDr)<<iboun<<" "<<setw(double_width)<<qleft[iboun]<<"\n";
			//fout <<setiosflags(ios::scientific)<<setiosflags(ios::left)<<setprecision(5);
                        fout2<<"qleft2: "<<setw(width_numbDr)<<iboun<<" "<<setw(double_width)<<qleft2[iboun]<<"\n";
                        //fout <<setiosflags(ios::scientific)<<setprecision(5);
                        //fout << setw(width_numbDr) << iboun <<" " << setw(width_numcell) << elemboundary(iboun)<<" "<<infec<<" "<<" eb: "<<setw(width_numcell)<<eboundary(iboun)<<" "<<setw(width_numcell) <<eboundary(iboun+numb)<<" "<<setw(width_numcell)<<eboundary(iboun+2*numb)<<" e: "<<setw(double_width)<<ee_n(elemboundary(iboun), infec, ibandc)<< " "<<setw(double_width)<<
                        //ee_n(elemboundary(iboun), infec, ibandc)* weight(infec)* ss(infec, 0) * vg_n[ibandc]<<"  qlef:    "<<setw(double_width)<<qleft[iboun]<<"\n" ;
		            }
                }
            }
        }

/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// Heat flux                                           Right
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if (eboundary(iboun + 2 * numbDr) == nright){
            qleft[iboun] = qright[iboun] = qright2[iboun] = qtop[iboun] = qbot[iboun] = 0;
            for (int ibandc = 0; ibandc < nband; ibandc++){
                for (int infec = 0; infec < nftot; infec++){
                    qright[iboun] += R_n[ibandc]*
                    ee_n(elemboundary(iboun), infec, ibandc)* weight(infec) * ss(infec, 0) * vg_n[ibandc];
                    if (sweight(infec, 0)<=0){
                        qright2[iboun]=qright2[iboun]+R_n[ibandc]*C_n[ibandc]/(4*PI)*(Tright-Tref)*sweight(infec, 0)*vg_n[ibandc];
                    }
                    else if (sweight(infec, 0)>0){
                        qright2[iboun]=qright2[iboun]+R_n[ibandc]*ee_n(elemboundary(iboun),infec,ibandc)*sweight(infec, 0)*vg_n[ibandc];
                    }
                    if ((ibandc == nband - 1) && (infec == nftot - 1)){
                        qrightsum=qrightsum+qright[iboun];
                        counterr=counterr+1;
                        fout<<"qright: "<<setw(width_numbDr)<<iboun<<" "<<setw(double_width)<<qright[iboun]<<"\n";
			//fout <<setiosflags(ios::scientific)<<setiosflags(ios::left)<<setprecision(5);
                        fout2<<"qright2:"<<setw(width_numbDr)<<iboun<<" "<<setw(double_width)<<qright2[iboun]<<"\n";
                        //fout <<setiosflags(ios::scientific)<<setprecision(5);
                        //fout << setw(width_numbDr) << iboun <<" " << setw(width_numcell) << elemboundary(iboun)<<" "<<infec<<" "<<" eb: "<<setw(width_numcell)<<eboundary(iboun)<<" "<<setw(width_numcell) <<eboundary(iboun+numb)<<" "<<setw(width_numcell)<<eboundary(iboun+2*numb)<<" e: "<<setw(double_width)<<ee_n(elemboundary(iboun), infec, ibandc)<< " "<<setw(double_width)<<
                        //ee_n(elemboundary(iboun), infec, ibandc)* weight(infec)* ss(infec, 0) * vg_n[ibandc]<<"  qlef:    "<<setw(double_width)<<qleft[iboun]<<"\n" ;
		            }
		        }
            }
        }

/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// Heat flux                                             Top
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if (eboundary(iboun + 2 * numbDr) == ntop){
            qleft[iboun] = qright[iboun] = qtop[iboun] = qbot[iboun] = 0;
            for (int ibandc = 0; ibandc < nband; ibandc++){
                for (int infec = 0; infec < nftot; infec++){
                    if (sweight(infec,1)>=0){
                        qtop[iboun] += R_n[ibandc] * ee_n(elemboundary(iboun), infec, ibandc) * weight(infec) * ss(infec, 1) * vg_n[ibandc];
                    }
                    else if (sweight(infec,1)<0){
                        int ntheta=1, nphi=1;
                        for (int nt = 0; nt < ntheta; nt++){
                            for (int np = 0; np < 4 * nphi; np++){
                                int nf1 = np + nt * 4 * nphi;
                                int npref, nfref;
                                if (infec == nf1) {
                                    cout<<"";
                                    int np2 = np + 1;
                                    int nt2 = nt + 1;
                                    if (np2 > nphi && np2 <= 2 * nphi){
                                        npref = 2 * nphi + 1 - np2;
                                    }
                                    else {
                                        npref = 6 * nphi + 1 - np2;
                                    }
                                    nfref = npref + (nt2 - 1) * 4 * nphi - 1;
                                    qtop[iboun] += R_n[ibandc] * ee_n(elemboundary(iboun), nfref, ibandc) * weight(infec) * ss(infec, 1) * vg_n[ibandc];
                                }
                            }
                        }
                    }
                    if ((ibandc == nband - 1) && (infec == nftot - 1)) {
			   	        //fout<<setiosflags(ios::scientific)<<setprecision(4);
                        //fout << setw(width_numbDr) <<iboun <<" " << setw(width_numcell) <<elemboundary(iboun)<<" "<<infec<<" "<<" eb:"<<setw(width_numcell)<<eboundary(iboun)<<" " <<setw(width_numcell)<<eboundary(iboun+numb)<<" "<<setw(width_numcell)<<eboundary(iboun+2*numb)<<"e"<<setw(double_width)<<ee_n(elemboundary(iboun), infec, ibandc)<< " "<<setw(double_width)<<
                        //ee_n(elemboundary(iboun), infec, ibandc)* weight(infec)* ss(infec, 1) * vg_n[ibandc]<<"  qtop: "<<setw(double_width)<<qtop[iboun]<<"\n" ;
			        }
		        }
            }
        }
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// Heat flux                                          Bottom
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if (eboundary(iboun + 2 * numbDr) == nbottom){
            qleft[iboun] = qright[iboun] = qtop[iboun] = qbot[iboun] = 0;
            for (int ibandc = 0; ibandc < nband; ibandc++){
                for (int infec = 0; infec < nftot; infec++){
                    if (sweight(infec,1)>=0){
                        qbot[iboun] += R_n[ibandc] * ee_n(elemboundary(iboun), infec, ibandc) * weight(infec) * ss(infec, 1) * vg_n[ibandc];
                    }
                    else if (sweight(infec,1)<0){
                        int ntheta=1, nphi=1;
                        for (int nt = 0; nt < ntheta; nt++){
                            for (int np = 0; np < 4 * nphi; np++){
                                int nf1 = np + nt * 4 * nphi;
                                int npref, nfref;
                                if (infec == nf1) {
                                    cout<<"";
                                    int np2 = np + 1;
                                    int nt2 = nt + 1;
                                    if (np2 > nphi && np2 <= 2 * nphi){
                                        npref = 2 * nphi + 1 - np2;
                                    }
                                    else {
                                        npref = 6 * nphi + 1 - np2;
                                    }
                                    nfref = npref + (nt2 - 1) * 4 * nphi - 1;
                                    qbot[iboun] += R_n[ibandc] * ee_n(elemboundary(iboun), nfref, ibandc) * weight(infec) * ss(infec, 1) * vg_n[ibandc];
                                }
                            }
                        }
                    }
                    if ((ibandc == nband - 1) && (infec == nftot - 1)) {
			   	        //fout<<setiosflags(ios::scientific)<<setprecision(4);
                        //fout << setw(width_numbDr) <<iboun <<" " << setw(width_numcell) <<elemboundary(iboun)<<" "<<infec<<" "<<" eb:"<<setw(width_numcell)<<eboundary(iboun)<<" " <<setw(width_numcell)<<eboundary(iboun+numb)<<" "<<setw(width_numcell)<<eboundary(iboun+2*numb)<<"e"<<setw(double_width)<<ee_n(elemboundary(iboun), infec, ibandc)<< " "<<setw(double_width)<<
                        //ee_n(elemboundary(iboun), infec, ibandc)* weight(infec)* ss(infec, 1) * vg_n[ibandc]<<"  qbot: "<<setw(double_width)<<qbot[iboun]<<"\n" ;
 			        }
		        }
            }
        }
    }
    delete[] qtop;
    delete[] qbot;
    delete[] qleft;
    delete[] qright;
    fout<<"qleftsum:   "<< qleftsum/counterl<<"\n";
    fout<<"qrightsum:  "<< qrightsum/counterr<<"\n";
    fout<<"Resistance: "<< (Tleft-Tright)/(((qleftsum/counterl)+(qrightsum/counterr))/2)<<"\n";

}

// calculate the error between results from
// current iteration and the previous iteration
double get_error(int nband, int nftot, CMDArray<double>& e1_n, CMDArray<double>& e2_n) {
    double error=0;
    for (int iband = 0; iband < nband; iband++){
        for (int inf10 = 0; inf10 < nftot; inf10++){
            for (int inn10 = 0; inn10 < numcell; inn10++)
                error += pow((e2_n(inn10, inf10, iband) - e1_n(inn10, inf10, iband)), 2);
        }
    }
    return (error / nband);
}

// check the convergence
bool check_convergence(double error, int nftot) {
    return (error <= numnode * nftot * 0.01);
}

// output ee
void output_ee(CMDArray<double>& ee_n, CMDArray<double>& e0_n, int flag, int nband, int nftot) {
    ofstream fout, fout1;
    //cout << nband << endl;
    if (flag) {
    	fout.open("eeinputfilec.dat");
        fout1.open("e0inputfilec.dat");
    }
    else {
	fout.open("eeinputfilew.dat");
	fout1.open("e0inputfilew.dat");
    }
    ostringstream output, output1;
    int i, j, k;
    for (i = 0; i < nband; i++)
	for (j = 0; j < nftot; j++) {
            for (k = 0; k < numcell; k++)
                output << setprecision(12) << ee_n(k, j, i) << ' ';
    	output << "\n";
    }
    fout << output.str();
    fout.close();
    for (i = 0; i < nband; i++) {
        for (k = 0; k < numcell; k++)
             output1 << setprecision(12) << e0_n(k, i) << ' ';
	output1 << "\n";
    }
    fout1 << output1.str();
    fout1.close();
}

/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
///                                                 MAIN CODE
/// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
int main (int argc, char* argv[]) {
    if (argc != 2) {
	cout << "Wrong number of arguments!." << endl;
	return 0;
    }
    // Get constants.
    double L_n, Tleft, Tright, Tref, WFACTOR;
    int ntheta, nphi, max_iter, nband;
    string str;
    char new_line;
    ifstream fin_const("inputdata.dat");
    ofstream fout_ret, fout_ret2, fout_node, fout_cell;
    bool flag = (((string)argv[1]) == "1");
    getline(fin_const, str);
    fin_const >> L_n >> new_line;
    getline(fin_const, str);
    fin_const >> ntheta >> new_line;
    getline(fin_const, str);
    fin_const >> nphi >> new_line;
    getline(fin_const, str);
    fin_const >> WFACTOR >> new_line;
    getline(fin_const,str);
    fin_const >> max_iter >> new_line;
    getline(fin_const,str);
    fin_const >> Tleft >> new_line;
    getline(fin_const,str);
    fin_const >> Tright >> new_line;
    getline(fin_const,str);
    fin_const >> Tref >> new_line;
    // Non-gray BTE - Different Bands.
    getline(fin_const, str);
    fin_const >> nband >> new_line;
    //cout << nband << "\n";
    double *tau_n = new double[nband];
    double *vg_n = new double[nband];
    double *Kn_n = new double[nband];
    double *R_n = new double[nband];
    double *C_n = new double[nband];
    double *mfp_n = new double[nband];
    string buf;
    stringstream buffer;
    getline(fin_const, str);
    //cout << str << "\n";
    getline(fin_const, buf);
    //cout << buf << "\n";
    buffer.str(buf);
    for (int i = 0; i < nband; i++)
        buffer >> *(vg_n + i);
    getline(fin_const, str);
    //cout << str << "\n";
    getline(fin_const, buf);
    //cout << buf << "\n";
    buffer.clear();
    buffer.str(buf);
    for (int i = 0; i < nband; i++)
        buffer >> *(R_n + i);
    getline(fin_const, str);
    getline(fin_const, buf);
    buffer.clear();
    buffer.str(buf);
    for (int i = 0; i < nband; i++)
        buffer >> *(tau_n + i);
    getline(fin_const, str);
    getline(fin_const, buf);
    buffer.clear();
    buffer.str(buf);
    for (int i = 0; i < nband; i++)
        buffer >> *(C_n + i);
    
    for (int iibb = 0;iibb < nband; iibb++) {
	mfp_n[iibb] = vg_n[iibb] * tau_n[iibb];         
	Kn_n[iibb] = mfp_n[iibb] / L_n;
    }
    // Boundary numbers.
    int ntop, ndifftop;
    int nbottom, ndiffbottom;
    int nright;
    int nleft;
    int ntopi, ndifftopi;
    int nbottomi, ndiffbottomi;
    int nrighti, ndiffrighti;
    int nrighti1, ndiffrighti1;
    int nlefti, ndifflefti;
    getline(fin_const,str);
    fin_const >> ntop >> ndifftop >> new_line;
    getline(fin_const,str);
    fin_const >> nbottom >> ndiffbottom >> new_line;
    getline(fin_const,str);
    fin_const >> nright >> new_line;
    getline(fin_const,str);
    fin_const >> nleft >> new_line;
    getline(fin_const,str);
    fin_const >> ntopi >> ndifftopi >> new_line;
    getline(fin_const,str);
    fin_const >> nbottomi >> ndiffbottomi >> new_line;
    getline(fin_const,str);
    fin_const >> nrighti >> ndiffrighti >> new_line;
    getline(fin_const,str);
    fin_const >> nrighti1 >> ndiffrighti1 >> new_line;
    getline(fin_const,str);
    fin_const >> nlefti >> ndifflefti >> new_line;
    double x_l = L_n;//y_l = 1.0 * x_l;
    ifstream fin_mesh("inputmesh.dat");
    // get node position
    for (int i = 0; i < 18; i++)
        getline(fin_mesh, str);
    fin_mesh >> numnode;
    CMDArray<double> p(2 * numnode);
    for (int i = 0; i < 4; i++)
        getline(fin_mesh, str);
    for (int i = 0; i < numnode; i++) {
        getline(fin_mesh, str);
        stringstream sss;
        sss.str(str);
        string coord_x, coord_y;
        sss >> coord_x >> coord_y;
        p(i) = strtod(coord_x.c_str(), NULL) * x_l;
        p(i + numnode) = strtod(coord_y.c_str(), NULL) * x_l;
    }
    // get node position
    int numb_ori;
    for (int i = 0; i < 37; i++)
        getline(fin_mesh, str);
    fin_mesh >> numb_ori;
    for (int i = 0; i < 2; i++)
        getline(fin_mesh, str);
    set<int> omit;
    if (!flag) {
	while (getline(fin_const, str)) {
	    if (str == "without")
		break;
	}
    }
    getline(fin_const, str);
    getline(fin_const, str);
    stringstream sss;
    sss.str(str);
    int num;
    while (sss >> num) {
        omit.insert(num);
	//cout << num << endl;
    }
    int* equi = new int[numb_ori];
    for (int i = 0; i < numb_ori; i++)
        *(equi + i) = i;
    getline(fin_const, str);
    while (getline(fin_const, str)) {
        stringstream sse;
        sse.str(str);
        int x, y;
        sse >> x;
        while (sse >> y) {
            Union(x, y, equi);
        }
    }
    CMDArray<int> eboundary_ori(3 * numb_ori);
    for (int i = 0; i < numb_ori; i++) {
        int x, y;
        fin_mesh >> x >> y;
        eboundary_ori(i) = x; eboundary_ori(i + numb_ori) = y;
    }
    for (int i = 0; i < 4; i++)
        getline(fin_mesh, str);
    numb = 0;
    for (int i = 0; i < numb_ori; i++) {
        getline(fin_mesh, str);
        int x = atoi(str.c_str());
        if (find_omit(x, omit))
            eboundary_ori(i + 2 * numb_ori) = -1;
        else {
            numb ++;
            eboundary_ori(i + 2 * numb_ori) = equi[x];
        }
    }
    CMDArray<int> eboundary(3 * numb);
    for (int i = 0, ii = 0; i < numb_ori; i++)
        if (eboundary_ori(i + 2 * numb_ori) != -1) {
            eboundary(ii) = eboundary_ori(i);
            eboundary(ii + numb) = eboundary_ori(i + numb_ori);
            eboundary(ii + 2 * numb) = eboundary_ori(i + 2 * numb_ori);
            ii++;
        }
    // get cell position
    for (int i = 0; i < 7; i++)
        getline(fin_mesh, str);
    fin_mesh >> numcell;
    CMDArray<int> t(3 * numcell);
    for (int i = 0; i < 2; i++)
        getline(fin_mesh, str);
    for (int i = 0; i < numcell; i++) {
        int x, y, z;
        fin_mesh >> x >> y >> z;
        t(i) = x; t(i + numcell) = y; t(i + 2 * numcell) = z;
    }
    fin_const.close();
    fin_mesh.close();
    delete[] equi;
    int numbDr = numb;
    int nftot = 4 * ntheta * nphi;
    CMDArray<double> sweight(nftot, 3), ss(nftot, 3), weight(nftot);
    get_direction(nftot, ntheta, nphi, WFACTOR, sweight, ss, weight);
    CMDArray<int> elmnod(numnode, numcell), nelemnode(numnode);
    CMDArray<double> node_r_m1T(numnode), Cc(numcell, 3);
    find_interpolation_element(elmnod, nelemnode, node_r_m1T, Cc, p, t);
    CMDArray<int> elemboundary(numb);
    get_boundary_element(elemboundary, eboundary, t, numbDr);
    CMDArray<double> ee_n(numcell, nftot, nband);
    CMDArray<double> e1_n(numcell, nftot, nband);
    CMDArray<double> e2_n(numcell, nftot, nband);
    CMDArray<double> Temp_n(numcell,nband);
    CMDArray<double> e0_n(numcell,nband), Temp(numcell);
    ee_n.set_zero();
    initialize_e0(ee_n, e0_n, Temp, Temp_n, Cc, weight, C_n, Tref, L_n, nband, nftot, flag);
    CMDArray<double> Temnode(numnode), Temnode_n(numnode);
    const int nedge = 3;
    cout << "Ncell: " << numcell << " Kn: " << Kn_n[0] << " nt: " << ntheta << "\n";
    clock_t t_start = clock();
    // Iteration
    int iter_num;
    for (iter_num = 0; iter_num < max_iter; iter_num ++){
        // Initialization
        for (int ibnde = 0; ibnde < nband; ibnde++)
            for (int infe = 0; infe < nftot; infe++)
                for (int inn = 0; inn < numcell; inn++)
                    e1_n(inn,infe,ibnde) = ee_n(inn,infe,ibnde);
        ee_n.set_zero();
        e2_n.set_zero();
        // solve for each band
        for (int iband = 0; iband < nband; iband++){
            for (int inf = 0; inf < nftot; inf++) {
                CMDArray<double> Ke(numcell, numcell), Re(numcell);
                double xe[2 * 4];
                get_coefficient(nedge, iter_num, iband, inf, numbDr, xe, mfp_n, C_n,
                                Tleft, Tright, Tref,L_n, ntheta, nphi, ntop, ntopi,
                                nbottom, nbottomi, nleft, nlefti, nright, nrighti,
                                nrighti1, p, t, sweight, ss, weight, Cc, eboundary,
                                e0_n, e1_n, Ke, Re);
                vector<double> sol = solve_matrix(Ke, Re);
                for (int id = 0; id < numcell; id++) {
                    ee_n(id, inf, iband) = sol[id];
                    e2_n(id, inf, iband) = sol[id];
                }
            }
	    if (flag)
		fout_cell.open("ResultTempCell-constriction.dat", ofstream::out);
    	    else fout_cell.open("ResultTempCell-without.dat", ofstream::out);
            get_cell_temp(iband, nftot, Cc, ee_n, weight, Temp_n, C_n, Tref, fout_cell);
	    fout_cell.close();
        }//iband
        // post-processing
        recover_temp(e0_n, Temp_n, Temp, nband, R_n, C_n, Tref);
	if (flag) {
	    fout_ret.open("Results-constriction.dat", ofstream::out);
	    fout_ret2.open("Results2-constriction.dat", ofstream::out);
	    fout_node.open("ResultTempCell-constriction.dat", ofstream::out);
	}
    	else {
	    fout_ret.open("Results-without.dat", ofstream::out);
	    fout_ret2.open("Results2-without.dat", ofstream::out);
	    fout_node.open("ResultTempCell-without.dat", ofstream::out);
	}
        interpolation(nband, Temp, Temnode_n, node_r_m1T, elmnod, nelemnode, Cc, p, t, fout_node);
        get_heat_transfer_flux(C_n,Tleft,Tright,Tref,sweight, eboundary, elemboundary, ee_n, weight, ss, R_n, vg_n,
        nband, nftot, numbDr, ntop, nbottom, nleft, nright, fout_ret, fout_ret2);
        fout_ret.close();
        fout_ret2.close();
        fout_node.close();
        
        double error = get_error(nband, nftot, e1_n, e2_n);
        cout << "Iter,Err_all = " << iter_num << " " << error << "\n";
        if (check_convergence(error, nftot)) {
	    output_ee(ee_n, e0_n, flag, nband, nftot);
            break;
        }
        if ((iter_num + 1) % 500 == 0)
	    output_ee(ee_n, e0_n, flag, nband, nftot);
    }//iter_num
    if (iter_num == max_iter)
	output_ee(ee_n, e0_n, flag, nband, nftot);
    clock_t t_end = clock();
    cout << "Time used: " << (double)(t_end - t_start) / CLOCKS_PER_SEC << "\n";
    delete[] tau_n;
    delete[] vg_n;
    delete[] Kn_n;
    delete[] R_n;
    delete[] C_n;
    delete[] mfp_n;
    return 0;
}

