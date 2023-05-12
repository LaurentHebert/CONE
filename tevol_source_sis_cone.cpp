/**
 * @file   tevol_source_sis_cone.cpp
 * @brief  ODE system for SIS dynamics on the CONE
 *
 * Source code. Single parameter passed as argument, precision and output format hardcoded.
 * g++ -std=c++11 -O3 -o tevol_source_sis_cone ./tevol_source_sis_cone.cpp $(gsl-config --cflags) $(gsl-config --libs)
 *
 * @author  LHD
 * @since   2021-02-06
 */

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>

#include <boost/multi_array.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "dyn_sis_cone.hpp"
#include "undirected_graph_t.hpp" //to extract matrices

double binom(int n, int k, double p);

using namespace std;

typedef vector< int > v;
typedef vector< double > vd;
typedef vector< v > vv;
typedef vector< vd > vvd;

// special comparison
struct classcomp {
  bool operator() (const pair< pair< pair<int,int>, pair<int,int> >, int >& lhs, const pair< pair< pair<int,int>, pair<int,int> >, int >& rhs) const
  {return lhs.first<rhs.first;}
};

int main(int argc, const char *argv[]) {
		 
	//Model parameters	
	double lambda = atof(argv[1]); //transmission rate
    double epsilon = 1e-4; //initial condition

    if(argc<1) {cerr << "Requires one parameter:\n(1) normalized transmission rate (lambda)"<< endl; return 0;}
  
    // Read the nodetype file and create the system
    ifstream input0("./matrices/Pnkl.txt");
    int ntemp, ktemp, ltemp, ctemp; double Ptemp; long long int N = 0;
    if (!input0.is_open()) {
	    cerr<<"expected node type matrix not found."<<endl;;
 	    return EXIT_FAILURE;
    } //end if
    else while ( input0 >> ntemp >> ws >> ktemp >> ws >> ltemp >> ws >> ctemp >> ws >> Ptemp ) N+= Ptemp;
    ifstream Pnkl("./matrices/Pnkl.txt");
    vector< pair < tuple<int,int,int>,double > > NodeTypes; 
    vector<int> c;  
    vector<double> Nl;
    map< pair<int,int>, int > Nlk;
    // Temporary read to prepare dimensions
    int nb_types=0; int maxl=0; int maxn=0;
    if (!Pnkl.is_open()) {
	    cerr<<"expected node type matrix not found."<<endl;;
 	    return EXIT_FAILURE;
    } //end if
    else {
	    while ( Pnkl >> ntemp >> ws >> ktemp >> ws >> ltemp >> ws >> ctemp >> ws >> Ptemp ) {
            ++nb_types; 
            if(ltemp>maxl) {maxl = ltemp; c.resize(maxl+1); Nl.resize(maxl+1);}
            if(ntemp>maxn) {maxn = ntemp;}
            c[ltemp]=ctemp;
            Nl[ltemp]+=Ptemp;
            Nlk[make_pair(ltemp,ktemp)] = Nlk[make_pair(ltemp,ktemp)] + Ptemp;
            NodeTypes.push_back(make_pair(make_tuple(ntemp,ktemp, ltemp),1.0*Ptemp/N));
        }//end while
    }//end else
    maxl = maxl + 1;
    maxn = maxn + 1;
    //cout << "read Pnkl.dat" << endl;

    // Read and create the matrix
    ifstream LMat("./matrices/Lmat.txt");    
    map< tuple<int,int,int,int>, int> L;
    int l1, k1, l2, k2; int Lval;
    vector<int> sumg(nb_types,0); vector<int> sumr(nb_types,0);
    set< pair<int,int> > typeset;
    if (!LMat.is_open()) {
	    cerr<<"expected link matrix not found."<<endl;;
 	    return EXIT_FAILURE;
    } //end if
    else {
	    while ( LMat >> k1 >> ws >> l1 >> ws >> k2 >> ws >> l2 >> ws >> Lval ) {
            //Assume L is a matrix in # of edges (not stubs)
            if(l1==l2 && k1==k2) L[make_tuple(l1,k1,l2,k2)] = 2*Lval;
            else {
                L[make_tuple(l1,k1,l2,k2)] = Lval;
                L[make_tuple(l2,k2,l1,k1)] = Lval;
            }
            typeset.insert(make_pair(l1,k1)); typeset.insert(make_pair(l2,k2));
        }
    }//end else

    // Calculate link normalisations
   vector< vector<double> > NormsRGB(nb_types, vector<double>(3,0.0));
   for(int it = 0; it < NodeTypes.size(); ++it) {
        tuple<int,int,int> nodetype = NodeTypes[it].first; Ptemp = NodeTypes[it].second;
        k1 = get<1>(nodetype); l1 = get<2>(nodetype);
        //for(int it2 = 0; it2 < NodeTypes.size(); ++it2) {
        for(set< pair<int,int> >::iterator it2 = typeset.begin(); it2 != typeset.end(); ++it2) {
            //tuple<int,int,int> nodetype2 = NodeTypes[it2].first;
            //k2 = get<1>(nodetype2); l2 = get<2>(nodetype2);
            k2 = (*it2).second; l2 = (*it2).first;
            try{
                L.at(make_tuple(l1,k1,l2,k2)) == 0;
            } catch(const out_of_range &e) {
              continue;
            }
            if(l2<l1-1) NormsRGB[it][0] += L[make_tuple(l1,k1,l2,k2)];
            if(l2==l1-1) NormsRGB[it][1] += L[make_tuple(l1,k1,l2,k2)];
            if(l2>=l1) NormsRGB[it][2] += L[make_tuple(l1,k1,l2,k2)];
            if(l2<l1-1) sumg[it] += L[make_tuple(l1,k1,l2,k2)];
            if(l2>=l1) sumr[it] += L[make_tuple(l1,k1,l2,k2)];
        }
    }

    // Computer anchoring probabilities
    vector< tuple<double,double,double> > degrees; degrees.resize(nb_types);
    for(int it = 0; it < nb_types; ++it) {
        tuple<int,int,int> nodetype = NodeTypes[it].first; Ptemp = NodeTypes[it].second;
        int k = get<1>(nodetype); int l = get<2>(nodetype);
        double pr = 1.0*sumr[it]/(Nlk[make_pair(l,k)]*c[l]);
        double pg = 0.0;
        if(k-c[l]>0) pg = 1.0*sumg[it]/(Nlk[make_pair(l,k)]*(k-c[l]));
        double kr = 1.0*c[l]*pr;
        double kg, kb;
        if(k-c[l]<=1) {
            kg = 1.0*(k-c[l])*pg - 1.0*int(c[l]==c[l-1])*pow(pr,c[l])*(k-c[l])*pg;
            kb = 1.0*c[l]*(1-pr) + 1.0*(k-c[l])*(1-pg) + 1.0*int(c[l]==c[l-1])*pow(pr,c[l])*(k-c[l])*pg;
        }
        else {
            kg = 1.0*(k-c[l])*pg;
            kb = 1.0*c[l]*(1-pr) + 1.0*(k-c[l])*(1-pg);
        }
        degrees[it] = make_tuple(kr, kg, kb);    }
    /*cout << "computed average degrees" << endl;
    for(int it = 0; it < nb_types; ++it) {
        cout << get<0>(degrees[it]) << " " << get<1>(degrees[it]) << " "<< get<2>(degrees[it]) << endl;
    }*/
    

    const int dim = nb_types;
    const int dim2 = maxn;
    Sparam param = {lambda, c, L, degrees, NormsRGB, NodeTypes, maxn, dim};
    //cout << "preparing " << dim2 << " times " << dim << " equations" << endl;

    // Integrator parameters
    double t = 0;
    double dt = 1e-5;
    double t_step = 1.0;
    const double eps_abs = 1e-9;
    const double eps_rel = 1e-9;

    // Setting initial conditions
    typedef boost::multi_array<double,2> mat_type;
    typedef mat_type::index index;
    mat_type y(boost::extents[dim2][dim]);
    fill(y.data(),y.data()+y.num_elements(),0.0);

    // Initial conditions
    int current_type = 0; double lastI = 0.0;
    for(int it = 0; it < NodeTypes.size(); ++it) {
        tuple<int,int,int> nodetype = NodeTypes[it].first; Ptemp = NodeTypes[it].second;
        int n = get<0>(nodetype);
        for(int i=0; i<=n; ++i) { y[i][current_type] = binom(n,i,epsilon)*Ptemp; lastI += i*y[i][current_type]; }
        ++current_type;
    }//end for node types

    double newI = 0.0;
    double newS = 1.0;

    //initial state output
    /*cout << t;
    for(int it = 0; it < nb_types; ++it) {
        tuple<int,int,int> nodetype = NodeTypes[it].first; Ptemp = NodeTypes[it].second;
        int n = get<0>(nodetype);
        for(int i=0; i<=n; ++i) cout << " " << y[i][it];
        cout << "\n";
    }*/

    // Define GSL odeiv parameters
    const gsl_odeiv_step_type * step_type = gsl_odeiv_step_rkf45;
    gsl_odeiv_step * step = gsl_odeiv_step_alloc (step_type, maxn*dim);
    gsl_odeiv_control * control = gsl_odeiv_control_y_new (eps_abs,eps_rel);
    gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc (maxn*dim);
    gsl_odeiv_system sys = {dydt, NULL, maxn*dim, &param};

	//Integration
    int status(GSL_SUCCESS);
    double diff1 = 1.0; double Ptemp1 = NodeTypes[0].second; double lastI1 = y[1][0]/Ptemp1;
    double diff2 = 1.0; double Ptemp2 = NodeTypes[NodeTypes.size()-1].second; double lastI2 = y[1][NodeTypes.size()-1]/Ptemp2;
    //for (double t_target = t+t_step; t_target < 1000; t_target += t_step ) { //stop by time
    for (double t_target = t+t_step; abs(diff1) > 1e-8 || abs(diff2) > 1e-8; t_target += t_step ) { //stop by difference
        while (t < t_target) {
            status = gsl_odeiv_evolve_apply (evolve,control,step,&sys,&t,t_target,&dt,y.data());
            if (status != GSL_SUCCESS) {
				cout << "SNAFU" << endl;
                break;
			}
        } // end while
        /*        
        //current state output
        for(int it = 0; it < NodeTypes.size(); ++it) {
            tuple<int,int,int> nodetype = NodeTypes[it].first; Ptemp = NodeTypes[it].second;
            int n = get<0>(nodetype); int k = get<1>(nodetype); int l = get<2>(nodetype);
            for(int i=0; i<=n; ++i) cout << t << " " << n << " " << k << " " << l << " " << i << " " << y[i][it]/Ptemp << "\n";
        } //end output 
        */
        //cout << t << " " << y[0][0]/(NodeTypes[0].second) << "\n";
        tuple<int,int,int> nodetype1 = NodeTypes[0].first; tuple<int,int,int> nodetype2 = NodeTypes[NodeTypes.size()-1].first;
        int n1 = get<0>(nodetype1); int n2 = get<0>(nodetype2);
        Ptemp1 = NodeTypes[0].second; diff1 = -lastI1; lastI1 = 0.0; for(int i=0; i<=n1; ++i) lastI1 += (y[i][0]*i)/(Ptemp1*n1); diff1 += lastI1;
        Ptemp2 = NodeTypes[NodeTypes.size()-1].second; diff2 = -lastI2; lastI2 = 0.0; for(int i=0; i<=n2; ++i) lastI2 += (y[i][NodeTypes.size()-1]*i)/(Ptemp2*n2); diff2 += lastI2;
	} //end while

    //final state output
    cout << lambda << " ";
    for(int it = 0; it < NodeTypes.size(); ++it) {
        tuple<int,int,int> nodetype = NodeTypes[it].first; Ptemp = NodeTypes[it].second;
        int n = get<0>(nodetype); int k = get<1>(nodetype); int l = get<2>(nodetype);
        //for(int i=0; i<=n; ++i) cout << lambda << " " << n << " " << k << " " << l << " " << i << " " << y[i][it]/Ptemp << "\n";
        double localS = 0.0;
        for(int i=0; i<=n; ++i) localS += (y[i][it]*i)/(Ptemp*n);
        cout << localS << " ";
        //cout << lambda << " " << n << " " << k << " " << l << " " << localS << "\n";
    }
    cout << "\n";

    cout.flush();

    // Free memory
    gsl_odeiv_evolve_free(evolve);
    gsl_odeiv_control_free(control);
    gsl_odeiv_step_free(step);
    
    return 0;
}

double binom(int n, int k, double p)
	{
	//check
	if(k>n || k<0 || n<0) return 0.0;
		
	double result=0.L;
	double postfix = pow(p,k)*pow(1.L-p,n-k);
	if(n<2*k) k=n-k;
	
	if(k<=100&&n<=200)
		{
		double numerator=1;
		for(int i=0; i<k; i++) numerator = numerator*((double)(n-i));
		double denominator=1;
		for(int i=2; i<=k; i++) denominator = denominator*((double)i);
		result = (double)(numerator/denominator)*postfix;
		}
	else
		{
		if(k<100||n<300)
			{
			result=1.L;
			for(int i=1; i<=k; i++) result = result*((double)(n-k+i))/((double)i);
			result = result*postfix;
			}
		else if(n<1000)
			{
			double nd=(double)n, kd=(double)k, nkd=(double)(n-k);
			result = nd*log(nd)-nd+0.5L*log(2.L*M_PI*nd)+1.L/(12.L*nd)-1.L/(360.L*pow(nd,3))+1.L/(1260.L*pow(nd,5));
			result -= kd*log(kd)-kd + 0.5L*log(2.L*M_PI*kd) + 1.L/(12.L*kd) - 1.L/(360.L*pow(kd,3)) + 1.L/(1260.L*pow(kd,5));
			result -= nkd*log(nkd)-nkd + 0.5L*log(2.L*M_PI*nkd) + 1.L/(12.L*nkd) - 1.L/(360.L*pow(nkd,3)) + 1.L/(1260.L*pow(nkd,5));
			result = exp(result);
			result = result*postfix;
			}
			
		else result = (1.L/sqrt(2.L*M_PI*(double)n*p*(1.L-p)))*exp(-pow(k-n*p,2)/(2.L*(double)n*p*(1.L-p)));
		}
	if(result!=result) result=0.L;
	if(isinf(result)) result=0.L;
	return result;
	}
