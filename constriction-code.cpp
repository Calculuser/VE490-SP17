//#In the name of Allah the most merciful, the most compassionate, who can help people to be and to do the best
///#
///##
///#### Non-gray BTE-FVM-Heat Transfer Cell_Centered Triang-Cells
///#####         Constriction ADIABATIC Diffusive Bc's
///######                   2017-3-6 -Comsol-mesh3-452 cells
///#######                  Modified by Jiaqi Zuo
///##################################################################
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <petscksp.h>
#include "mdarray.h"
using namespace std;

vector<double> gauss(vector< vector<double> > A);
int numb, numnode, numcell;
const double PI = 4.0 * atan(1.0);

void get_inputfile_data (CMDArray<double>& data_matrix) {
    ifstream fin("inputfile.dat");
    int in_row=0;
    string str;
    while (getline(fin, str)){
        stringstream ss;
        ss.str(str);
        if(str.at(0) >= 'a' && str.at(0) <= 'z') {
            char c;
            ss >> c;
            switch (c){
                case 'e':
                    ss >> numb;
                    break;
                case 'p':
                    ss >> numnode;
                    break;
                case 't':
                    ss >> numcell;
                    break;
                default:
                    break;

            }
        }
        else {
            int in_col = 0;
            double data_in;
            while (ss >> data_in) {
                data_matrix(in_row, in_col) = data_in;
                ++in_col;
            }
        }
        ++in_row;
    }
    //cout << numb << ' ' << numnode << ' ' << numcell << endl;
    fin.close();
}

void get_node_position (CMDArray<double>& p,
                        CMDArray<double>& data_matrix,
                        const double x_l) {
    int iprow=0;
    for (int ipn = 9; ipn < 11; ipn++){
        for (int ip = 0; ip < numnode; ip++){
            int ipnode=ip + numnode * iprow;
            p(ipnode) = data_matrix(ipn, ip) * x_l;
        }
        ++iprow;
    }
}

int get_boundary_position (CMDArray<double>& data_matrix,
                            CMDArray<int>& BcDirich,
                            CMDArray<int>& eboundary,
                            CMDArray<int>& t) {
    int numbDr = 0, ip1 = 0, ipDr = 0, xxa;
    ///.Boundary Condition numbers For Drichlet Boundary condition can be 1,2..,..
    for (int ip = 0; ip < numb; ip++){
        if (((int)data_matrix(5, ip)) >= 1 &&
                ((int)data_matrix(5, ip)) <= 10){
            numbDr++;
            xxa = (int)data_matrix(1, ip);
            BcDirich(xxa) = 1 ;
            xxa = (int)data_matrix(2, ip);
            BcDirich(xxa)=1;
            //cout<<" "<<xxa<< " "<<BcDirich[xxa]<<endl;
        }
    }
    for (int ip = 0; ip < numbDr; ip++){
        ///. input file for Drichlet Boundary Condition numbers can be 1,2,3,4
        if (((int)data_matrix(5, ip)) >= 1 &&
            ((int)data_matrix(5, ip)) <= 10){
            ip1 = ipDr + numbDr * 0;
            eboundary(ip1) = (int)data_matrix(1, ip);
            ip1 = ipDr + numbDr * 1;
            eboundary(ip1) = (int)data_matrix(2, ip);
            ip1 = ipDr + numbDr * 2;
            eboundary(ip1) = (int)data_matrix(5, ip);
            ip1 = ipDr + numbDr * 3;
            eboundary(ip1) = (int)data_matrix(7, ip);
            ++ipDr;
        }
        ///. input file for Boundary Condition numbers For Neumann  Boundary condition are 0,0
        ///. input file for Neumann Boundary Condition numbers are 0,0
    }
    for (int ip = 0; ip < numcell; ip++){
        ip1 = ip + numcell * 0;
        t(ip1) = (int)data_matrix(12, ip);
        ip1 = ip + numcell * 1;
        t(ip1) = (int)data_matrix(13, ip);
        ip1 = ip + numcell * 2;
        t(ip1) = (int)data_matrix(14, ip);
        ip1 = ip + numcell * 3;
        t(ip1) = (int)data_matrix(15, ip);
    }
    return numbDr;
}

void get_direction(int nftot, int ntheta, int nphi, double WFACTOR,
                   CMDArray<double>& sweight, CMDArray<double>& ss,
                   CMDArray<double>& weight) {
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
            sweight(nf, 0) = WFACTOR * sin(phi(np)) * sin(0.5 * dphi)
                             * (dtheta - cos(2 * theta(nt)) * sin(dtheta));
            sweight(nf, 1) = WFACTOR * cos(phi(np)) * sin(0.5 * dphi)
                             * (dtheta - cos(2 * theta(nt)) * sin(dtheta));
            sweight(nf, 2) = WFACTOR * 0.5 * dphi * sin(2 * theta(nt)) * sin(dtheta);
            ss(nf, 0) = sin(theta(nt)) * sin(phi(np));
            ss(nf, 1) = sin(theta(nt)) * cos(phi(np));
            ss(nf, 2) = cos(theta(nt));
        }
    }
}

void find_interpolation_element (CMDArray<int>& elmnod,
                                 CMDArray<int>& nelemnode,
                                 CMDArray<double>& node_r_m1T,
                                 CMDArray<double>& Cc,
                                 CMDArray<double>& p,
                                 CMDArray<int>& t) {
    /// Finding center cell and area A of cell [..][2]
    for (int elem1=0; elem1 < numcell; elem1++){
        Cc(elem1, 0) = (p(t(elem1)) + p(t(elem1 + numcell))
                        + p(t(elem1 + 2 * numcell))) / 3.0;
        Cc(elem1, 1) = (p(t(elem1) + numnode) + p(t(elem1 + numcell) + numnode)
                        + p(t(elem1 + 2 * numcell) + numnode)) / 3.0;
        Cc(elem1, 2) = 0.5 * (p(t(elem1)) * (p(t(elem1 + numcell) + numnode)
                                             -p(t(elem1 + 2 * numcell) + numnode))
                              + p(t(elem1 + numcell)) * (p(t(elem1 + 2 * numcell) + numnode)
                                                         -p(t(elem1) + numnode))
                              + p(t(elem1 + 2 * numcell)) *(p(t(elem1) + numnode)
                                                            -p(t(elem1 + numcell) + numnode)));
    }
    // Here, I delete the area_node[], as it's useless from the code point of view.
    for (int elemnode0 = 0; elemnode0 < numnode; elemnode0++){
        nelemnode(elemnode0) = 0;
        for (int elemnode1=0; elemnode1 < numcell; elemnode1++){
            if (elemnode0 == t(elemnode1)
                || elemnode0 == t(elemnode1 + numcell)
                || elemnode0 == t(elemnode1 + 2 * numcell)){
                elmnod(elemnode0, nelemnode(elemnode0)) = elemnode1;
                nelemnode(elemnode0)++;
                node_r_m1T(elemnode0) += 1.0 / sqrt(pow((Cc(elemnode1, 0) - p(elemnode0)), 2)
                                                     + pow((Cc(elemnode1, 1) - p(elemnode0 + numnode)), 2));
                // cout<< elemnode0 <<" "<<elemnode1 <<" "<< t[elemnode1] <<" "<<pow((Cc[elemnode1][0]-p[elemnode0]),2) <<" "<<pow((Cc[elemnode1][1]-p[elemnode0+numnode]),2)<< " "<<node_r_m1T[elemnode0]<<endl;
            }
        }
    }
}

void get_boundary_element(CMDArray<int>& elemboundary,
                          CMDArray<int>& eboundary,
                          CMDArray<int>& t, int numbDr) {
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

void initialize_e0 (CMDArray<double>& e0, CMDArray<double>& Temp,
                    double Tleft, double Tright) {
    for (int i = 0; i < numcell; i++){
        Temp(i) = Tleft;
        e0(i) = 0;                                 //ctot*(Temp[i]-Tref)/(4*PI);
    }
    Temp(1) = Tright; Temp(4) = Tright; Temp(5) = Tright;
    Temp(2) = Tright; Temp(9) = Tright; Temp(8) = Tright;
    Temp(6) = Tright; Temp(7) = Tright;
}

void get_coefficient(const int nedge, int iband, int inf, int numbDr,
                     double xe[], double* vg_n, double* C_n, double* tau_n,
                     double Tleft, double Tright, double Tref,
                     int ntheta, int nphi, int ntop, int ntopi,
                     int nbottom, int nbottomi, int nleft,
                     int nlefti, int nright, int nrighti, int nrighti1,
                     CMDArray<double>& p, CMDArray<int>& t,
                     CMDArray<double>& sweight, CMDArray<double>& ss,
                     CMDArray<double>& weight, CMDArray<double>& Cc,
                     CMDArray<int>& eboundary, CMDArray<double>& e0_n,
                     CMDArray<double>& e1_n, CMDArray<double>& Ke,
                     CMDArray<double>& Re) {

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
        /// cout<<endl;
        ///  cout<<"ie   :"<<ie<<endl;
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
                                neledge[iedge][1] == eboundary(ib + numb))||
                            (neledge[iedge][0] == eboundary(ib + numb) &&
                                    neledge[iedge][1] == eboundary(ib))){
                            ineighb=0;
                            break;
                            //cout<<" rrrr "<<eboundary[ib+2*numb]<<" "<<eboundary[ib+2*numb]<<" "<< eboundary[ib+2*numb]<<endl;
                        }
                    }
                    if (ib < numbDr)
                        break;
                }
            }
            leng_edge = sqrt(pow(xe[iedge + 1] - xe[iedge], 2) +
                                     pow(xe[iedge + 1 + nedge + 1] - xe[iedge + 1 + nedge], 2));
            norm[0] = (xe[iedge + 1 + nedge + 1] - xe[iedge + 1 + nedge]) / leng_edge;
            norm[1] = -(xe[iedge + 1] - xe[iedge]) / leng_edge;
            //cout<<ie<<" "<<iedge<< " "<<xe[iedge+1]<<" "<<xe[iedge]<<" "<<xe[iedge+1+nedge+1]<<" "<<xe[iedge+1+nedge]<<" " <<leng_edge<<endl;
/// find e0 and temperature
            if (ineighb){ //means neighbor cell exists and it is not boundary
                swn = (ss(inf, 0) * norm[0] + ss(inf, 1) * norm[1]);
                a_f[iedge] = vg_n[iband] *
                        (sweight(inf, 0) * norm[0] + sweight(inf, 1) * norm[1]) * leng_edge;
                if (swn>=0)
                    Ke(ie, ie) += a_f[iedge];
                else Ke(ie, eneighb) += a_f[iedge];
            }
//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/// Boundaries
//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            else { //means neighbor cell not exists and it is boundary
                for (int ib = 0; ib < numbDr; ib++){
                    if (((neledge[iedge][0] == eboundary(ib) &&
                            neledge[iedge][1] == eboundary(ib + numb))) ||
                            ((neledge[iedge][0] == eboundary(ib + numb) &&
                                    neledge[iedge][1] == eboundary(ib)))){
                        swn = (ss(inf, 0) * norm[0] + ss(inf, 1) * norm[1]);
///! outgoing
                        a_f[iedge] = vg_n[iband] *
                                (sweight(inf, 0) * norm[0] +
                                        sweight(inf, 1) * norm[1]) * fabs(leng_edge);

/// ******************************************
/// Specular boundary y = L  (ebound 2) top
/// ******************************************
/// DIFFUSIVE BOUNDARY
                        //cout<<ib<< " " <<eboundary[ib+numbDr*2]<< " "<<endl;
                        if (eboundary(ib + numbDr * 2) == ntopi){
                            //cout<<"ntop"<<ib<<" "<<ie<<" "<<eboundary[ib+numbDr*2]<<endl;
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
                                    if (swn >= 0){                    ///!         outgoing
                                        Ke(ie, ie) += a_f[iedge];
                                    }
                                    else if (swn < 0) {                ///!         incoming
                                        Re(ie) -= einsum * vg_n[iband] *
                                                (sweight(inf, 0) * norm[0] +
                                                        sweight(inf, 1) * norm[1]) * leng_edge;
                                    }
                                }
                            }

                        }
                        //cout<<ie<< " "<< ib<< " " <<eboundary[ib+numbDr*2]<< " "<<endl;
                        else if (eboundary(ib + numbDr * 2) == ntop){
                            //cout<<"ntop"<<ib<<" "<<eboundary[ib+numbDr*2]<<endl;
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
                                            Re(ie) -= e1_n(ie, nfref, iband) *
                                                    vg_n[iband] * (sweight(inf, 0) * norm[0] +
                                                    sweight(inf, 1) * norm[1]) * leng_edge;
                                        }
                                    }
                                }
                            }
                        }
/// ******************************************
/// Specular boundary y = 0 (ebound 7) bottom
/// ******************************************
/// DIFFUSIVE BOUNDARY
                        else if (eboundary(ib + numbDr * 2) == nbottomi){
                            // cout<<"nboti"<<ib<<" "<<eboundary[ib+numbDr*2]<<endl;
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
                                    if (swn >= 0){                      ///!         outgoing
                                        Ke(ie, ie) += a_f[iedge];
                                    }
                                    else if (swn<0) {                ///!         incoming
                                        Re(ie) -= einsum*vg_n[iband] *
                                                (sweight(inf, 0) * norm[0] +
                                                        sweight(inf, 1) * norm[1]) * leng_edge;
                                    }
                                }
                            }
                        }
/// SPECULAR BOUNDARY
                        else if (eboundary(ib + numbDr * 2) == nbottom){
                            //cout<<"nbot"<<ib<<" "<<eboundary[ib+numbDr*2]<<endl;
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
                                            if (np2 <= nphi){
                                                npref = 2 * nphi + 1 - np2;
                                            }
                                            else {
                                                npref = 6 * nphi + 1 - np2;
                                            }
                                            nfref = npref + (nt2 - 1) * 4 * nphi - 1;
                                            Re(ie) -= e1_n(ie, nfref, iband) *
                                                    vg_n[iband] * (sweight(inf, 0) * norm[0] +
                                                    sweight(inf, 1) * norm[1]) * leng_edge;
                                        }
                                    }
                                }
                            }
                        }
/// *************************************************************************************************************

/// *******************************************************
///Diffusive OR Specular boundary  (ebound 4 internal-top)
/// *******************************************************
                        else if (eboundary(ib + numbDr * 2) == ntopi){
// DIFFUSIVE BOUNDARY
                            // cout<<"ntopi"<<ib<<" "<<eboundary[ib+numbDr*2]<<endl;
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
                                    if (swn >= 0){                      ///!         outgoing
                                        Ke(ie, ie) += a_f[iedge];
                                    }
                                    else if (swn < 0) {                ///!         incoming
                                        Re(ie) -= einsum * vg_n[iband] *
                                                (sweight(inf, 0) * norm[0] +
                                                        sweight(inf, 1) * norm[1]) * leng_edge;
                                    }
                                }
                            }
// SPECULAR BOUNDARY
//                                       for (int nt=0;nt<ntheta;nt++){
//                                          for (int np=0;np<4*nphi;np++){
//                                             int nf1 = np +(nt)*4*nphi;
//                                             int npref,nfref;
//                                             if (inf==nf1) {
//                                                cout<<"";
//                                                if (swn>=0){       //         outgoing
//                                                   Ke[ie][ie]+=a_f[iedge];
//                                                   cout<<"";
//                                                }
//                                                else if (swn<0) {  //         incoming
//                                                   int np2=np+1;
//                                                   int nt2=nt+1;
//                                                   if ((np2)<=nphi){
//                                                      npref = 2*nphi+1 -np2;
//                                                   }
//                                                   else {
//                                                      npref = 6*nphi+1-np2;
//                                                   }
//                                                   nfref=npref +(nt2-1)*4*nphi-1;
//                                                   Re[ie]+=-e1_n(ie,nfref,iband)*vg_n[iband]*(sweight[inf][0]*norm[0]+sweight[inf][1]*norm[1])*leng_edge;
//                                                }
//                                             }
//                                          }
//                                       }
                        }
/// *******************************************************
///Diffusive OR Specular boundary  (ebound 5 internal-bottom)
/// *******************************************************

                        else if (eboundary(ib + numbDr * 2) == nbottomi){
// DIFFUSIVE BOUNDARY
                            // cout<<"nboti"<<ib<<" "<<eboundary[ib+numbDr*2]<<endl;
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
                                    if (swn >= 0){                      ///!         outgoing
                                        Ke(ie, ie) += a_f[iedge];
                                    }
                                    else if (swn < 0) {                ///!         incoming
                                        Re(ie) -= einsum * vg_n[iband] *
                                                (sweight(inf, 0) * norm[0] +
                                                        sweight(inf, 1) * norm[1]) * leng_edge;
                                    }
                                }
                            }
// SPECULAR BOUNDARY

//                                       for (int nt=0;nt<ntheta;nt++){
//                                          for (int np=0;np<4*nphi;np++){
//                                             int nf1 = np +(nt)*4*nphi;
//                                             int npref,nfref;
//                                             if (inf==nf1) {
//                                                cout<<"";
//                                                if (swn>=0){       //         outgoing
//                                                   Ke[ie][ie]+=a_f[iedge];
//                                                   cout<<"";
//                                                }
//                                                else if (swn<0) {  //         incoming
//                                                   int np2=np+1;
//                                                   int nt2=nt+1;
//                                                   if (np2>nphi&&np2<=2*nphi){
//                                                      npref = 2*nphi+1 -np2;
//                                                   }
//                                                   else {
//                                                       npref = 6*nphi+1-np2;
//                                                   }
//                                                   nfref=npref +(nt2-1)*4*nphi-1;
//                                                   Re[ie]+=-e1_n(ie,nfref,iband)*vg_n[iband]*(sweight[inf][0]*norm[0]+sweight[inf][1]*norm[1])*leng_edge;
//                                                }
//                                             }
//                                          }
//                                       }
                        }
/// *******************************************************
///Diffusive OR Specular boundary  (ebound 3 internal-right)
/// *******************************************************
                        else if (eboundary(ib + numbDr * 2) == 11) {
// DIFFUSIVE BOUNDARY
                            /// cout<<"nrigi"<<ib<<" "<<eboundary[ib+numbDr*2]<<endl;
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
                                    if (swn >= 0){                      ///!         outgoing
                                        Ke(ie, ie) += a_f[iedge];
                                    }
                                    else if (swn < 0) {                ///!         incoming
                                        Re(ie) -= einsum * vg_n[iband] *
                                                (sweight(inf, 0) * norm[0] +
                                                        sweight(inf, 1) * norm[1]) * leng_edge;
                                    }
                                }
                            }
// SPECULAR BOUNDARY

//                                       for (int nt=0;nt<ntheta;nt++){
//                                          for (int np=0;np<4*nphi;np++){
//                                             int nf1 = np +(nt)*4*nphi;
//                                             int npref,nfref;
//                                             if (inf==nf1) {
//                                                cout<<"";
//                                                if (swn>=0){       //         outgoing
//                                                   Ke[ie][ie]+=a_f[iedge];
//                                                   cout<<"";
//                                                }
//                                                else if (swn<0) {  //         incoming
//                                                   int np2=np+1;
//                                                   int nt2=nt+1;
//                                                   if (np2>nphi&&np2<=2*nphi){
//                                                      npref = 2*nphi+1 -np2;
//                                                   }
//                                                   else {
//                                                      npref = 6*nphi+1-np2;
//                                                   }
//                                                   nfref=npref +(nt2-1)*4*nphi-1;
//                                                   Re[ie]+=-e1_n(ie,nfref,iband)*vg_n[iband]*(sweight[inf][0]*norm[0]+sweight[inf][1]*norm[1])*leng_edge;
//                                                   cout<<"";
//                                                }
//                                             }
//                                          }
//                                       }

                        }

/// *******************************************************
///Diffusive OR Specular boundary  (ebound 8 internal-left)
/// *******************************************************
                        else if (eboundary(ib + numbDr * 2) == nlefti) {
// DIFFUSIVE BOUNDARY
                            ///cout<<"nlefi"<<ib<<" "<<eboundary[ib+numbDr*2]<<endl;
                            double einsum = 0.0;
                            for (int nft = 0; nft < nfmax; nft++){
                                swn = ss(nft, 0) * norm[0] + ss(nft, 1) * norm[1];
                                if (swn >= 0) {
                                    einsum = einsum + e1_n(ie, nft, iband) * sweight(nft, 0);
                                }
                            }
                            einsum = einsum / PI;
                            for (int nft2 = 0; nft2 < nfmax; nft2++){
                                swn = ss(nft2, 0) * norm[0] + ss(nft2, 1) * norm[1];
                                if (inf == nft2) {
                                    if (swn >= 0){                      ///!         outgoing
                                        Ke(ie, ie) += a_f[iedge];
                                    }
                                    else if (swn < 0){                ///!         incoming
                                        Re(ie) -= einsum * vg_n[iband] *
                                                (sweight(inf, 0) * norm[0] +
                                                        sweight(inf, 1) * norm[1]) * leng_edge;
                                    }
                                }
                            }
// SPECULAR BOUNDARY

//                                       for (int nt=0;nt<ntheta;nt++){
//                                          for (int np=0;np<4*nphi;np++){
//                                             int nf1 = np +(nt)*4*nphi;
//                                             int npref,nfref;
//                                             if (inf==nf1) {
//                                                cout<<"";
//                                                if (swn>=0){       //         outgoing
//                                                   Ke[ie][ie]+=a_f[iedge];
//                                                   cout<<"";
//                                                }
//                                                else if (swn<0) {  //         incoming
//                                                   int np2=np+1;
//                                                   int nt2=nt+1;
//                                                   if (np2>nphi&&np2<=2*nphi){
//                                                      npref = 2*nphi+1 -np2;
//                                                    }
//                                                   else {
//                                                      npref = 6*nphi+1-np2;
//                                                   }
//                                                   nfref=npref +(nt2-1)*4*nphi-1;
//                                                   Re[ie]+=-e1_n(ie,nfref,iband)*vg_n[iband]*(sweight[inf][0]*norm[0]+sweight[inf][1]*norm[1])*leng_edge;
//                                                }
//                                             }
//                                          }
///////////////////                                       }
                            //       cout<<ie<<" "<<einsum<<" "<<iband <<"   e-e0= "<<ie<< " "<<inf<<" " <<e1_n(ie,inf,iband)<<" "<<einsum<<" "<<e1_n(ie,inf,iband)*weight[inf]<<" "<<e0[ie]<<endl;
//////////////////
                        }
//Temperature constant boundary

/// incoming
                        else if (swn < 0) {
/// boundary 1  right
                            //cout<<"nr"<<nright<<" "<<nleft<<" "<<eboundary[ib+numbDr*2]<<endl;
                            if (eboundary(ib + numbDr * 2) == nright ||
                                    eboundary(ib + numbDr * 2) == nlefti * 2 ||
                                    eboundary(ib + numbDr * 2) == 40 ||
                                    eboundary(ib + numbDr * 2) == 50){
                                //  cout<<"nrig"<<ib<<" "<<eboundary[ib+numbDr*2]<<endl;
                                //cout<<"right"<<ib<< " " <<eboundary[ib+numbDr*2]<< " "<<endl;
                                Re(ie) -= C_n[iband] / (4 * PI) * (Tright - Tref) *
                                        vg_n[iband] * (sweight(inf, 0) * norm[0] +
                                        sweight(inf, 1) * norm[1]) * leng_edge;
                            }
///!         boundary 6 left
                            else if (eboundary(ib + numbDr * 2) == nleft ||
                                    eboundary(ib + numbDr * 2) == nbottomi){
                                //     cout<<"nlef"<<ib<<" "<<eboundary[ib+numbDr*2]<<endl;
                                //cout<<"left"<<ib<< " " <<eboundary[ib+numbDr*2]<< " "<<endl;
                                Re(ie) -= C_n[iband] / (4 * PI) * (Tleft - Tref) *
                                        vg_n[iband] * (sweight(inf, 0) * norm[0] +
                                        sweight(inf, 1) * norm[1]) * leng_edge;
                            }
                        }
                        else if (swn >= 0){
/// outgoing
                            Ke(ie, ie) += a_f[iedge]  ;
                        }
                    }
                }
            }
        }//iedge
        Ke(ie, ie) += 1 / tau_n[iband] * Cc(ie, 2) * weight(inf);
        Re(ie) += e0_n(ie,iband) / tau_n[iband] * Cc(ie, 2) * weight(inf);
    }//ie
}

vector<double> solve_matrix_gauss (CMDArray<double>& Ke, CMDArray<double>& Re) {

    CMDArray<double> Km_active(numcell * numcell);
    CMDArray<double> Re_active(numcell);

    vector<double> line(numcell + 1,0);
    vector< vector<double> > Km_Re(numcell, line);
    for (int i1 = 0; i1 < numcell; i1++) {
        for (int i2 = 0; i2 < numcell; i2++)
            Km_Re[i1][i2] = Ke(i1, i2);
            /// if (Ke[i1][i2]>0||Ke[i1][i2]<0)fout<<"KE= "<<i1<<" "<<i2<<" "<<Ke[i1][i2]<<endl;
            /// if (Ke[i1][i2]>0||Ke[i1][i2]<0)cout<<"KE= "<<i1<<" "<<i2<<" "<<Ke[i1][i2]<<endl;
        Km_Re[i1][numcell] = Re(i1);
    }
    vector<double> x(numcell) ;
    x = gauss(Km_Re);
    return x;
}

vector<double> solve_matrix_PETSc (CMDArray<double>& Ke, CMDArray<double>& Re) {
    int* nnz = new int[numcell], count;
    for (int i = 0; i < numcell; i++) {
        count = 0;
        for (int j = 0; j < numcell; j++)
            if (Ke(i, j) != 0)
                count++;
        *(nnz + i) = count;
    }
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
    VecDuplicate(bbb, &xxx); // create anothe solution vector (xxx)
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
    return x;
/// finished petsc
}

void get_cell_temp (int iband, int nftot, CMDArray<double>& Cc, CMDArray<double>& e0,
                    CMDArray<double>& ee_n, CMDArray<double>& weight,
                    CMDArray<double>& Temp_n, double* C_n, double Tref, ofstream& fout) {
    Temp_n.set_zero();
    e0.set_zero();
    /*fout << "iband: " << iband << endl;
    fout << "ee_n:" << endl;
    for (int i = 0; i < numcell; i++){
        for (int j = 0; j < nftot; j++)
            fout << ee_n(i, j, iband) << ' ';
        fout << endl;
    }*/

    for (int ic = 0; ic < numcell; ic++){
        for (int iinf = 0; iinf < nftot; iinf++){
            e0(ic) += ee_n(ic, iinf, iband) * weight(iinf);
            //fout << "ic: " << ic << ' ' << "iinf: " << iinf << endl;
            //fout << e0(ic) << ' ' << ee_n(ic, iinf, iband) << endl;
            /*
            if (ic==10||ic==327||ic==55||ic==80){
//////////////////
                // cout<<ic<<" "<<iband <<"   e-e0= "<<ic<< " "<<iinf<<" " <<ee_n(ic,iinf,iband)<<" "<<weight[iinf]<<" "<<ee_n(ic,iinf,iband)*weight[iinf]<<" "<<e0[ic]<<endl;
/////////////////
            }
             */
            /// fout<<"e-e0= "<<ic<< " "<<iinf<<" " <<ee[ic][iinf]<<endl;
        }
        Temp_n(ic, iband) = e0(ic) / C_n[iband] + Tref;
        /*
        if (ic==10||ic==327||ic==55||ic==80){
////////////////
                //    cout<<ic<<" "<<iband <<" "<<e0[ic]  <<" "<<Temp_n(ic,iband)<<endl;
////////////////
        }*/
        e0(ic) = e0(ic) / (4 * PI);
        ///fout_n<<"bnd-T"<<"   "<<ic<<"  "<<iband<<"  "<<e0[ic]<<"  "<<Temp_n(ic,iband)<<endl;
    }//ic
    /*
    fout << "e0:" << endl;
    for (int i = 0; i < numcell; i++)
        fout << e0(i) << ' ';
    fout << endl;*/
    /*fout << "Temp_n:" << endl;*/
    for (int i2 = 0; i2 < numcell; i2++)
        fout << i2 << ' ' << Cc(i2, 0) << ' '  << Cc(i2, 1)
             << ' ' << Temp_n(i2,iband) << endl;
///non_Gray-non_Gray-non_Gray-non_Gray-non_Gray-non_Gray-non_Gray-non_Gray-non_Gray-non_Gray-non_Gray-\\\.
}

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

void interpolation(int nband, CMDArray<double>& Temp,
                   CMDArray<double>& Temnode_n,
                   CMDArray<double>& node_r_m1T,
                   CMDArray<int>& elmnod, CMDArray<int>& nelemnode,
                   CMDArray<double>& Cc, CMDArray<double>& p,
                   CMDArray<int>& t, ofstream& fout) {
    /// Finding center cell and area A of cell [..][2]
    for (int elem1=0; elem1 < numcell; elem1++){
        Cc(elem1, 0) = (p(t(elem1)) + p(t(elem1 + numcell))
                        + p(t(elem1 + 2 * numcell))) / 3.0;
        Cc(elem1, 1) = (p(t(elem1) + numnode) + p(t(elem1 + numcell) + numnode)
                        + p(t(elem1 + 2 * numcell) + numnode)) / 3.0;
        Cc(elem1, 2) = 0.5 * (p(t(elem1)) * (p(t(elem1 + numcell) + numnode)
                                             -p(t(elem1 + 2 * numcell) + numnode))
                              + p(t(elem1 + numcell)) * (p(t(elem1 + 2 * numcell) + numnode)
                                                         -p(t(elem1) + numnode))
                              + p(t(elem1 + 2 * numcell)) *(p(t(elem1) + numnode)
                                                            -p(t(elem1 + numcell) + numnode)));
    }

    for (int iband = 0; iband < nband; iband++){
        for (int inode = 0; inode < numnode; inode++){
            Temnode_n(inode) = 0;
            //if (BcDirich[inode]!=1){
            for (int ino = 0; ino < nelemnode(inode); ino++){
                Temnode_n(inode) += Temp(elmnod(inode, ino))
                                    * (1.0 / (sqrt((pow((Cc(elmnod(inode, ino), 0)
                                                         - p(inode)), 2)
                                                    + pow((Cc(elmnod(inode, ino), 1)
                                                           -p(inode + numnode)), 2)))))
                                    /node_r_m1T(inode);
                //cout<<inode<<" "<<ino<<" "<<elmnod[inode][ino]<< " "<<Temp[elmnod[inode][ino]]<<" "<<Temnode[inode]<<endl;
                //cout<<(1/(   sqrt(  ( pow((Cc[elmnod[inode][ino]][0]-p[inode]),2)+pow((Cc[elmnod[inode][ino]][1]-p[inode+numnode]),2)  )   )))<<" "<<node_r_m1T[inode]<<" "<<endl;
                //cout<<"";
            }
            fout << inode << ' ' << p(inode) << " " << p(inode + numnode)
                 << " " << Temnode_n(inode) << endl;
            //}
            // else if (BcDirich[inode]==1){
            //     Temnode[inode]=TempDrich[inode];
            // cout<<p[inode]<<" "<<p[inode+numnode]<<" "<<Temnode[inode]<<endl;
            //}
        }
    }
}

void get_heat_transfer_flux(CMDArray<int>& eboundary,
                            CMDArray<int>& elemboundary,
                            CMDArray<double>& ee_n, CMDArray<double>& weight,
                            CMDArray<double>& ss, double* R_n, double* vg_n,
                            int nband, int nftot, int numbDr, int ntop,
                            int nbottom, int nleft, int nright, ofstream& fout) {
    double *qleft = new double[numb];
    double *qright = new double[numb];
    double *qtop = new double[numb];
    double *qbot = new double[numb];
///.//////////////////////////////////////////////////////////////////////////
///.//////////////////////////////////////////////////////////////////////////
///.//////// Calculation of Heat Transfer flux //////////.///
///1 top
    for (int iboun = 0; iboun < numbDr; iboun++)
        if (eboundary(iboun + 2 * numbDr) == ntop){
            qleft[iboun] = qright[iboun] = qtop[iboun] = qbot[iboun] = 0;
            for (int ibandc = 0; ibandc < nband; ibandc++){
                for (int infec = 0; infec < nftot; infec++){
                    // if (sweight[infec][1]>0){
                    //cout<<"toppp:"<<qtop[iboun]<<endl;
                    qtop[iboun] += R_n[ibandc] *
                            ee_n(elemboundary(iboun), infec, ibandc) *
                            weight(infec)* ss(infec, 0) * vg_n[ibandc];
                    //cout<<"top : "<<elemboundary[iboun]<<" "<<infec<<" "<<ibandc<<" "<<ee_n(elemboundary[iboun],infec,ibandc)<<" "<<qtop[iboun]<<" "<<ee_n(elemboundary[iboun],infec,ibandc)*weight[infec]* ss[infec][0]*vg_n[ibandc]<<endl;
                    if ((ibandc == nband - 1) && (infec == nftot - 1))
                        fout << "top : " << elemboundary(iboun) << " "
                             << infec << " " << ibandc << " "
                             << ee_n(elemboundary(iboun), infec, ibandc)
                             << " " << qtop[iboun] << " "
                             << R_n[ibandc] * ee_n(elemboundary(iboun), infec, ibandc) *
                                weight(infec) * ss(infec, 0) * vg_n[ibandc] << endl;
// }
                    // else{
                    //     qtop[iboun]=qtop[iboun]+C_n[ibandc]/(4*PI)*(Tright-Tref)*sweight[infec][1]*vg_n[ibandc]*L_n/10/(4*PI);
                    //     cout<<"1 en: "<<infec<<" "<<elemboundary[iboun]<<" "<<qtop[iboun]<<" "<<C_n[ibandc]/(4*PI)*(Tright-Tref)*sweight[infec][1]*vg_n[ibandc]*L_n/10/(4*PI)<<endl;
                    // }
                }
                //         cout<<"" ;
                //         cout<<"1 iboun= "<<iboun<<" ; ibandc= "<<ibandc<<" ;qtop= "<<qtop[iboun]<<endl;
            }
        }
///2 right
    for (int iboun = 0; iboun < numbDr; iboun++){
        if (eboundary(iboun + 2 * numbDr) == nright){
            qleft[iboun] = qright[iboun] = qtop[iboun] = qbot[iboun] = 0;
            for (int ibandc = 0; ibandc < nband; ibandc++){
                for (int infec = 0; infec < nftot; infec++){
//                     if (sweight[infec][0]>0){
                    qright[iboun] += R_n[ibandc] *
                            ee_n(elemboundary(iboun), infec, ibandc) *
                            weight(infec) * ss(infec, 0) * vg_n[ibandc];
                    //cout<<"righ: "<<elemboundary[iboun]<<" "<<infec<<" "<<ibandc<<" "<<ee_n(elemboundary[iboun],infec,ibandc)<<" "<<qright[iboun]<<" "<<ee_n(elemboundary[iboun],infec,ibandc)*weight[infec]* ss[infec][0]*vg_n[ibandc]<<endl;
                    if ((ibandc == nband - 1) && (infec == nftot - 1))
                        fout<< "righ: " << elemboundary(iboun)
                            << " " << infec << " " << ibandc << " "
                            << ee_n(elemboundary(iboun), infec, ibandc)
                            << " " << qright[iboun] << " "
                            << R_n[ibandc] * ee_n(elemboundary(iboun), infec, ibandc) *
                                    weight(infec) * ss(infec, 0) * vg_n[ibandc]<<endl;
//                     }
//                     else{
//                        qright[iboun]=qright[iboun]+C_n[ibandc]/(4*PI)*(Tright-Tref)*sweight[infec][0]*vg_n[ibandc];
                    //              cout<<"2  en: "<<infec<<" "<<elemboundary[iboun]<<" "<<qright[iboun]<<endl;
//                        fout<<"2  en: "<<infec<<" "<<elemboundary[iboun]<<" "<<qright[iboun]<<endl;
//                     }
                }
                //         cout<<"" ;
                //         cout<<"2 iboun= "<<iboun<<" ; ibandc= "<<ibandc<<" ;qright= "<<qright[iboun]<<endl;
                //fout<<"2 iboun= "<<iboun<<" ; ibandc= "<<ibandc<<" ;qright= "<<qright[iboun]<<endl;

            }
        }
    }
///3 bottom
    for (int iboun = 0; iboun < numbDr; iboun++){
        if (eboundary(iboun + 2 * numbDr) == nbottom){
            qleft[iboun] = qright[iboun] = qtop[iboun] = qbot[iboun] = 0;
            for (int ibandc = 0; ibandc < nband; ibandc++){
                for (int infec = 0; infec < nftot; infec++){
                    //if (sweight[infec][1]<0){
                    qbot[iboun] += R_n[ibandc] *
                            ee_n(elemboundary(iboun), infec, ibandc) *
                            weight(infec) * ss(infec, 0) * vg_n[ibandc];
                    //           cout<<"3    : "<<infec<<" "<<ee_n(elemboundary[iboun],infec,ibandc)<<" "<<elemboundary[iboun]<<" "<<qtop[iboun]<<" "<<ee_n(elemboundary[iboun],infec,ibandc)*sweight[infec][1]*vg_n[ibandc]<<endl;
                    if ((ibandc==nband-1) && (infec== nftot - 1))
                        fout << "bot : " << elemboundary(iboun) << " "
                             << infec << " " << ibandc << " "
                             << ee_n(elemboundary(iboun), infec, ibandc)
                             << " " << qbot[iboun] << " "
                             << R_n[ibandc] * ee_n(elemboundary(iboun), infec, ibandc) *
                                     ss(infec, 0) * vg_n[ibandc]<<endl;
                    //}
                    //else{
                    //   qbot[iboun]=qbot[iboun]+C_n[ibandc]/(4*PI)*(Tright-Tref)*sweight[infec][1]*vg_n[ibandc]*L_n/10/(4*PI);
                    //      cout<<"3 en: "<<iboun <<" "<<infec <<" "<<elemboundary[iboun]<<" "<<qbot[iboun]<<endl;
                    //}
                }
                //      cout<<"" ;
                //       // cout<<"3 iboun= "<<iboun<<" ; ibandc= "<<ibandc<<" ;qbot= "<<qbot[iboun]<<endl;

            }
        }
    }
///4 left
    for (int iboun = 0; iboun < numbDr; iboun++){
        if (eboundary(iboun + 2 * numbDr) == nleft){
            qleft[iboun] = qright[iboun] = qtop[iboun] = qbot[iboun] = 0;
            for (int ibandc = 0; ibandc < nband; ibandc++){
                for (int infec = 0; infec < nftot; infec++){
//                     if (sweight[infec][0]<0){
                    qleft[iboun] += R_n[ibandc] *
                            ee_n(elemboundary(iboun), infec, ibandc) *
                            weight(infec) * ss(infec, 0) * vg_n[ibandc];
                    //           //cout<<"4  ex: "<<infec<<" "<<ee_n(elemboundary[iboun],infec,ibandc)<<" "<<elemboundary[iboun]<<" "<<qtop[iboun]<<" "<<ee_n(elemboundary[iboun],infec,ibandc)*sweight[infec][1]*vg_n[ibandc]<<endl;
                    if ((ibandc == nband - 1) && (infec == nftot - 1))
                        fout << "left: " << elemboundary(iboun) << " "
                             << infec << " " << ibandc << " "
                             << ee_n(elemboundary(iboun), infec, ibandc)
                             << " " << qleft[iboun] << " "
                             << R_n[ibandc] * ee_n(elemboundary(iboun), infec, ibandc) *
                                     ss(infec, 0) * vg_n[ibandc]<<endl;
//                     }
//                     else{
//                        qleft[iboun]=qleft[iboun]+C_n[ibandc]/(4*PI)*(Tleft-Tref)*sweight[infec][0]*vg_n[ibandc];
                    //            //cout<<"4  en: "<<infec<<" "<<elemboundary[iboun]<<" "<<qleft[iboun]<<endl;
//                        fout<<"4  en: "<<infec<<" "<<elemboundary[iboun]<<" "<<qleft[iboun]<<endl;

//                     }
                }
                //       //cout<<"" ;
                //       //cout<<"4 iboun= "<<iboun<<" ; ibandc= "<<ibandc<<" ;qleft= "<<qleft[iboun]<<endl;
                // fout<<"4 iboun= "<<iboun<<" ; ibandc= "<<ibandc<<" ;qleft= "<<qleft[iboun]<<endl;

            }
        }
    }
    delete[] qtop;
    delete[] qbot;
    delete[] qleft;
    delete[] qright;
}

double get_error(int nband, int nftot,
               CMDArray<double>& e1_n, CMDArray<double>& e2_n) {
    double error=0;
    for (int iband = 0; iband < nband; iband++){
        for (int inf10 = 0; inf10 < nftot; inf10++){
            for (int inn10 = 0; inn10 < numcell; inn10++)
                error += pow((e2_n(inn10, inf10, iband) -
                                e1_n(inn10, inf10, iband)), 2);
                /// fout <<e1[inn10][inf10]<<"     "<<e2[inn10][inf10]<<"     "<<errorall<<endl;
                ///fout_n<<"Err"<<" "<<iband<<" "<<inf10<<" "<<inn10<<e2_n(inn10,inf10,iband)<<" "<< e1_n(inn10,inf10,iband)<<" "<<errorall<<endl;
        }
    }
    return (error / nband);
}

bool check_convergence(double error, int nftot) {
    return (error <= numnode * nftot * 0.01);
}

int main (void) {
    // Get constants.
    double L_n, Tleft, Tright, Tref, WFACTOR;
    int ntheta, nphi, max_iter, nband;
    string str;
    char new_line;
    ifstream fin_const("inputconst.dat");
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
    double *tau_n = new double[nband];
    double *vg_n = new double[nband];
    double *Kn_n = new double[nband];
    double *R_n = new double[nband];
    double *C_n = new double[nband];
    string dump;
    getline(fin_const, str);
    for (int i = 0; i < nband; i++)
        fin_const >> *(vg_n + i);
    fin_const >> new_line;
    if (new_line != '\n') {
        getline(fin_const, dump);
        fin_const >> new_line;
    }
    getline(fin_const, str);
    for (int i = 0; i < nband; i++)
        fin_const >> *(R_n + i);
    fin_const >> new_line;
    if (new_line != '\n') {
        getline(fin_const, dump);
        fin_const >> new_line;
    }
    getline(fin_const, str);
    for (int i = 0; i < nband; i++)
        fin_const >> *(tau_n + i);
    fin_const >> new_line;
    if (new_line != '\n') {
        getline(fin_const, dump);
        fin_const >> new_line;
    }
    getline(fin_const, str);
    for (int i = 0; i < nband; i++)
        fin_const >> *(C_n + i);
    fin_const >> new_line;
    if (new_line != '\n') {
        getline(fin_const, dump);
        fin_const >> new_line;
    }
    for (int iibb = 0;iibb < nband; iibb++)
        Kn_n[iibb] = vg_n[iibb] * tau_n[iibb] / L_n;

    // Boundary numbers.
    int ntop, ndifftop; //ndifftop=1 : Boundary diffusive
    int nbottom, ndiffbottom;
    int nright;
    int nleft;
    int ntopi, ndifftopi;
    int nbottomi, ndiffbottomi;
    int nrighti, ndiffrighti;
    int nrighti1, ndiffrighti1;
    int nlefti, ndifflefti;

    getline(fin_const, str);
    fin_const >> ntop >> ndifftop >> new_line;
    getline(fin_const, str);
    fin_const >> nbottom >> ndiffbottom >> new_line;
    getline(fin_const, str);
    fin_const >> nright >> new_line;
    getline(fin_const, str);
    fin_const >> nleft >> new_line;
    getline(fin_const, str);
    fin_const >> ntopi >> ndifftopi >> new_line;
    getline(fin_const, str);
    fin_const >> nbottomi >> ndiffbottomi >> new_line;
    getline(fin_const, str);
    fin_const >> nrighti >> ndiffrighti >> new_line;
    getline(fin_const, str);
    fin_const >> nrighti1 >> ndiffrighti1 >> new_line;
    getline(fin_const, str);
    fin_const >> nlefti >> ndifflefti >> new_line;

    int eboundary_num;
    getline(fin_const, str);
    fin_const >> eboundary_num >> new_line;

    bool is_Linux;
    getline(fin_const, str);
    fin_const >> is_Linux >> new_line;
    vector<double> (*solve_matrix)(CMDArray<double>&, CMDArray<double>&);
    if (is_Linux)
        solve_matrix = solve_matrix_PETSc;
    else solve_matrix = solve_matrix_gauss;

    CMDArray<double> data_matrix(20, eboundary_num);
    get_inputfile_data(data_matrix);
    CMDArray<double> p(4 * numnode);
    CMDArray<int> BcDirich(numnode);
    CMDArray<int> eboundary(eboundary_num);
    CMDArray<int> t(4 * numcell);
    double x_l = L_n, y_l = 1.0 * x_l;
    int numbDr;
    get_node_position(p, data_matrix, x_l);

    numbDr = get_boundary_position(data_matrix, BcDirich, eboundary, t);
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
    CMDArray<double> e0_n(numcell,nband);
    CMDArray<double> e0(numcell), Temp(numcell);
    initialize_e0(e0, Temp, Tleft, Tright);
    ee_n.set_zero();

    CMDArray<double> Temnode(numnode), Temnode_n(numnode);
    const int nedge = 3;

    cout << "Kn: " << Kn_n[0] << " nt: " << ntheta << endl;
    clock_t t_start = clock();
    // Iteration
    for (int iter=0; iter < max_iter; iter++){
        ofstream fout_ret("Results.dat");
        ofstream fout_node("ResultTempNode.dat");
        ofstream fout_cell("ResultTempCell.dat");
        // Initialization
        for (int ibnde = 0; ibnde < nband; ibnde++)
            for (int infe = 0; infe < nftot; infe++)
                for (int inn = 0; inn < numcell; inn++)
                    e1_n(inn,infe,ibnde) = ee_n(inn,infe,ibnde);
        ee_n.set_zero();
        e2_n.set_zero();
        // solve for each band.
        for (int iband = 0; iband < nband; iband++){
            for (int inf = 0; inf < nftot; inf++) {
                //cout << "iband: " << iband << ' ' << "inf: " << inf << endl;
                CMDArray<double> Ke(numcell, numcell), Re(numcell);
                double xe[2 * 4];
                get_coefficient(nedge, iband, inf, numbDr, xe, vg_n, C_n, tau_n,
                                Tleft, Tright, Tref, ntheta, nphi, ntop, ntopi,
                                nbottom, nbottomi, nleft, nlefti, nright, nrighti,
                                nrighti1, p, t, sweight, ss, weight, Cc, eboundary,
                                e0_n, e1_n, Ke, Re);
                /*cout << "Re:" << endl;
                for (int i = 0; i < numcell; i++)
                    cout << Re(i) << endl;*/
                //cout << "Solution:" << endl;
                vector<double> sol = solve_matrix(Ke, Re);
                for (int id = 0; id < numcell; id++) {
                    ee_n(id, inf, iband) = sol[id];
                    //cout << ee_n(id, inf, iband) << ' ';
                    e2_n(id, inf, iband) = sol[id];
                }
                //cout << endl;
            }
            /*ofstream fout("Other.dat");
            fout << "iband: " << iband << endl;
            cout << "ee_n:(before)" << endl;
            for (int i = 0; i < numcell; i++){
                for (int j = 0; j < nftot; j++)
                    cout << ee_n(i, j, iband) << ' ';
                cout << endl;
            }*/
            get_cell_temp(iband, nftot, Cc, e0, ee_n, weight, Temp_n, C_n, Tref, fout_cell);
        }//iband
        recover_temp(e0_n, Temp_n, Temp, nband, R_n, C_n, Tref);
        interpolation(nband, Temp, Temnode_n, node_r_m1T, elmnod, nelemnode, Cc, p, t, fout_node);
        /*fout_ret << "iter: " << iter << " ee_n:" << endl;
        for (int i = 0; i < numcell; i++) {
            for (int j = 0; j < nftot; j++)
                fout_ret << ee_n(numcell, nftot, 0) << ' ';
            fout_ret << endl;
        }*/
        get_heat_transfer_flux(eboundary, elemboundary, ee_n, weight, ss, R_n, vg_n,
                               nband, nftot, numbDr, ntop, nbottom, nleft, nright, fout_ret);
        fout_ret.close();
        fout_node.close();
        fout_cell.close();
        double error = get_error(nband, nftot, e1_n, e2_n);
        cout << "Iter,Err_all = " << iter << " " << error << endl;
        if (check_convergence(error, nftot))
            break;
    }//iter
    clock_t t_end = clock();
    cout << "Time used: " << (double)(t_end - t_start) / CLOCKS_PER_SEC
         << " sec" <<endl;
    delete[] tau_n;
    delete[] vg_n;
    delete[] Kn_n;
    delete[] R_n;
    delete[] C_n;
    return 0;
}

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
            //      cout<<maxEl<<"  "<<maxRow<<" "<<k<<" "<<Km_Re[k][i]<<endl;
            if (fabs(Km_Re[k][i]) > maxEl)
            {
                maxEl = fabs(Km_Re[k][i]);
                maxRow = k;
                //          cout<<maxEl<<" "<<maxRow<<" "<<k<<" "<<Km_Re[k][i]<<endl;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1; k++)
        {
            double tmp = Km_Re[maxRow][k];
            Km_Re[maxRow][k] = Km_Re[i][k];
            Km_Re[i][k] = tmp;
            //         cout<<k<<"  "<<Km_Re[i][k]<<" "<<Km_Re[i][i]<<" "<<endl;
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
//               cout<<"  "<<Km_Re[k][i]<<" "<<Km_Re[i][i]<<" "<<c<<endl;
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    vector<double> x(n);
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


///////////////////////////////////////////////////////////////////////////////////////////////////
/*
void print(vector< vector<double> > A)
{
    int n = A.size();
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n+1; j++)
        {
            //          cout << A[i][j] << "\t";
            if (j == n-1)
            {
                //              cout << "| ";
            }
        }
        //      cout << "\n";
    }
    cout << endl;
}
*/
