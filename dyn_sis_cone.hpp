#ifndef DYN_SIS_CONE_HPP_INCLUDED
#define DYN_SIS_CONE_HPP_INCLUDED

#include <boost/multi_array.hpp>

using namespace std;

/**
 * @file   dyn_sis_cone.hpp
 * @brief  ODE system for SIS on CONE networks
 *
 * @author  LHD
 * @since   2021-05-30
 */

struct Sparam {
    const double lambda;
    const vector<int> c;
    const map< tuple<int,int,int,int>, int> L;
    const vector< tuple<double,double,double> > degrees;
    const vector< vector<double> > NormsRGB;
    const vector< pair < tuple<int,int,int>,double > > NodeTypes;   
    const int maxn;
    const int dim;
}; // parameter structure

struct classcomp_types {
  bool operator() (const pair< pair< int,int >, int>& lhs, const pair< pair< int,int >, int>& rhs) const
  {return lhs.first<rhs.first;}
};

//********** function dydt definition **************************
int dydt(double t, const double y[], double f[], void * param) {

    // Cast parameters
    Sparam& p = *static_cast<Sparam* >(param);

    // Create multi_array reference to y and f
    typedef boost::multi_array_ref<const double,2> CSTmatref_type;
    typedef boost::multi_array_ref<double,2> matref_type;
    typedef CSTmatref_type::index indexref;
    CSTmatref_type yref(y,boost::extents[p.maxn][p.dim]);
    matref_type fref(f,boost::extents[p.maxn][p.dim]);

    // Loop over node types
    for(int it = 0; it < p.NodeTypes.size(); ++it) {
        tuple<int,int,int> nodetype = p.NodeTypes[it].first; double Ptemp = p.NodeTypes[it].second;
        int n = get<0>(nodetype); int k = get<1>(nodetype); int l = get<2>(nodetype);
        double kr = get<0>(p.degrees[it]); double kg = get<1>(p.degrees[it]); double kb = get<2>(p.degrees[it]);
        //calculate mean-fields
        double MF0 = 0.0; double MF1 = 0.0; double MF2 = 0.0;
        for(int it2 = 0; it2 < p.NodeTypes.size(); ++it2) {
          tuple<int,int,int> nodetype2 = p.NodeTypes[it2].first; double Ptemp2 = p.NodeTypes[it2].second;
          int n2 = get<0>(nodetype2); int k2 = get<1>(nodetype2); int l2 = get<2>(nodetype2);
          double kr2 = get<0>(p.degrees[it2]); double kg2 = get<1>(p.degrees[it2]); double kb2 = get<2>(p.degrees[it2]);
          try{
             p.L.at(make_tuple(l,k,l2,k2)) == 0;
          } catch(const out_of_range &e) {
             continue;
          }
          for(int i=0; i<=n2; ++i) {
              if(l2<l-1 && p.NormsRGB[it][0]!=0) MF0 += 1.0*kg*(1.0*i/n2)*p.lambda*p.L.at(make_tuple(l,k,l2,k2))*(yref[i][it2]/Ptemp2)/p.NormsRGB[it][0];
              if(l2==l-1 && p.NormsRGB[it][1]!=0) MF1 += 1.0*kb*(1.0*i/n2)*p.lambda*p.L.at(make_tuple(l,k,l2,k2))*(yref[i][it2]/Ptemp2)/p.NormsRGB[it][1];
              if(l2>=l && p.NormsRGB[it][2]!=0) MF2 += 1.0*kr*(1.0*i/n2)*p.lambda*p.L.at(make_tuple(l,k,l2,k2))*(yref[i][it2]/Ptemp2)/p.NormsRGB[it][2];
          } 
        }
        double MF = MF0 + MF1 + MF2;
        //calculate derivatives
        for(int i=0; i<=n; ++i) {
            fref[i][it] = -1.0*(n-i)*(p.lambda*i+MF/(1.0*n))*yref[i][it] - 1.0*i*yref[i][it];
            if(i>0) fref[i][it] += 1.0*(n-i+1)*(p.lambda*(i-1)+MF/(1.0*n))*yref[i-1][it];
            if(i<n) fref[i][it] += 1.0*(i+1)*yref[i+1][it];
        }
    }//end for node types

    return GSL_SUCCESS;

} //********** end function dydt definition ********************************************************

#endif // DYN_SI_CONE_HPP_INCLUDED
