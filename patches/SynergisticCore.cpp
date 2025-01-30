// [[Rcpp::plugins(cpp11)]]
// Standard library headers
#include <array>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <unordered_map>
#include <unordered_set>

// External library headers
#include <Rcpp.h>
#include <boost/functional/hash.hpp>
using namespace Rcpp;

template <class T>
inline void hash_combine(std::size_t& seed, T const& v)
{
  seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{  
  template<typename T, size_t N>
  struct hash<array<T, N> >
  {
    typedef array<T, N> argument_type;
    typedef size_t result_type;
    result_type operator()(const argument_type& a) const
    {
      hash<T> hasher;
      result_type h = 0;
      for (result_type i = 0; i < N; ++i)
	{
	  h = h * 31 + hasher(a[i]);
	}
      return h;
    }
  };
  
  // How to use std::vector? //
  template<typename T>
  struct hash<vector<T> >
  {
    typedef vector<T> argument_type;
    typedef size_t result_type;
    result_type operator()(const argument_type& a) const
    {
      size_t size = a.size();
      size_t seed = 0;
      for (size_t i = 0; i < size; i++)
	//Combine the hash of the current vector with the hashes of the previous ones
	hash_combine(seed, a[i]);
      return seed;
    }    
  };
}

// [[Rcpp::export]]
NumericVector naomitInt(NumericVector x){
  std::vector<int> r(x.size());
  int k=0;
  for(int i = 0; i < x.size(); ++i) {
    if(x[i]==x[i]) {
      r[k] = x[i];
      k++;
    }
  }
  r.resize(k);
  return Rcpp::wrap(r);    
}

// [[Rcpp::export]]
NumericVector percentile_rcpp_sorted(NumericVector& x, NumericVector& percentile){
  int size_per = percentile.size();
  NumericVector percentile_vec = no_init(size_per);
  for (int ii = 0; ii < size_per; ii++){
    uint temp = 0.5 + ( double ) (x.size() - 1) * percentile[ii];
    percentile_vec[ii] = x[temp];
  }
  return percentile_vec;
}

// [[Rcpp::export]]
NumericVector computeJHandTC(NumericMatrix& y3_sub, NumericMatrix& H_sub){  
  float nr = y3_sub.nrow();
  int nl = y3_sub.ncol();
  NumericVector out(4);
  out[1] = 1.0; //jh_max
  
  // FREQUENCY TABLE
  std::unordered_map< std::vector<int>, float> freqs;
  freqs.reserve(nr+1);
  for(int i=0; i < nr; ++i){
    std::vector<int> row;
    int vr = nr*nl+1;
    row.reserve(vr); 
    for(int j=0; j < nl; ++j)
      row.emplace_back(y3_sub(i,j));
    freqs[row] += 1/nr;    
  }
  
  // Entropy
  for (auto it : freqs){
    if (it.second != 0.0){
      out[0] += it.second * log2(it.second);
    }
  }
  out[0] = -out[0];
  if(out[0] < 0){
    out[0] = 0.0;
    out[1] = 0.0;
  } else {  
    // Normalize by max jh
    for(int i=0; i < nl; ++i){
      out[1] *= unique(y3_sub(_,i)).size();
    }
    if(out[1] > 1){
      out[1] = out[0] / log2(out[1]);
    }else{
      out[1] = out[0] / log2(out[1]+1);
    }
  }
  
  // Compute TC
  for(int i=0; i < H_sub.nrow(); ++i){
    out[2] += H_sub(i,0); //H
    out[3] += H_sub(i,2); //H.maxNormalized
  }
  out[2] = out[2] - out[0]; // TC
  out[3] = out[3] - out[1]; // TC.maxNormalized
    
  return out;  
}

///////////////////////////////////////////////////////////////////////
///////////////////// HEURISTIC ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector myseq(int first, int last){
  NumericVector y(abs(last - first) + 1);
  if (first < last) 
    std::iota(y.begin(), y.end(), first);
  else {
    std::iota(y.begin(), y.end(), last);
    std::reverse(y.begin(), y.end());
  }
  return y;
}

// [[Rcpp::export]]
List compute_MMI_JSD(NumericMatrix& newCOM, NumericMatrix& y3, NumericMatrix& H, StringVector& GEN)
{
  // for each newCOM, compute JH and TC and store in TC    
  // use std::vector<int, newCOM[i,]> as the key
  std::cout << "starting c++ part" << "\n";  
  std::unordered_map<std::vector<int>, std::array<float, 10> > TC;    
  int nr = newCOM.nrow();
  int nl = newCOM.ncol();
  TC.reserve(nr+10000); //initial combis + 10 variable combinations
  int topN_TC_heuristic = 9; //not in use
  float cutoff_synergy = 0;
  //int max_add = 11;
  int max_add = 2; //+3
  for(int i=0; i < nr; ++i){
    NumericVector id = naomitInt(newCOM(i,_));
    int id_size = id.size();
    NumericMatrix y3_sub(y3.nrow(),id.size());
    NumericMatrix H_sub(id_size,H.nrow());
    for(int j=0; j < id_size; ++j){
      y3_sub(_,j) = y3(_,id(j));
      H_sub(j,_) = H(id(j),_);
    }

    // compute JH and TC
    std::array<float, 10> jh; 
    NumericVector temp = computeJHandTC(y3_sub, H_sub);
    for(int j=0; j < temp.size(); ++j){
      jh[j] = temp[j];
    }    
    std::vector<int> id_vec;
    int vr = id_size+1;
    id_vec.reserve(vr); 
    for(int j=0; j < id_size; ++j)
      id_vec.emplace_back(id(j));
    TC.insert(std::make_pair(id_vec, jh));
  }
  
  // Compute MMI, II
  std::map<float, std::vector<std::vector<int> > > orderedMMI; //for sorting
  for(int i=0; i < nr; ++i){
    NumericVector id = naomitInt(newCOM(i,_)); 
    int id_size = id.size();
    NumericVector MMI(4);
    
    // MMI, MMI.maxNormalized, II, II.maxNormalized
    // Add single entropies
    for(int j=0; j < id_size; ++j){
      MMI[0] += H(id(j),0); //MMI
      MMI[1] += H(id(j),2); //MMI.maxNormalized
    }
    MMI[2] = -(pow(-1, (id_size-1))) * MMI[0]; //II
    MMI[3] = -(pow(-1, (id_size-1))) * MMI[1]; //II.maxNomalized
	       
    // compute MMI and II
    int n = id.size();  
    for(int k=2; k <= n; ++k){
      std::vector<bool> v(n);
      std::fill(v.end() - k, v.end(), true);
      do {
	std::vector<int> id_vec;
	int vr = n+1;
	id_vec.reserve(vr); 
	for (int i=0; i < n; ++i) {
	  if (v[i]) { 
	    id_vec.emplace_back(id(i)); 	    
	  }
	}
	// 1. MMI
	MMI[0] += -(pow(-1, k)) * TC[id_vec][0];
	MMI[1] += -(pow(-1, k)) * TC[id_vec][1];

	// 2. II
	MMI[2] += -(pow(-1, (id_size-k))) * TC[id_vec][0];
	MMI[3] += -(pow(-1, (id_size-k))) * TC[id_vec][1];
	
      } while (std::next_permutation(v.begin(), v.end()));
    }
    // Add elements to TC
    std::vector<int> id_vec;
    id_vec.reserve(id_size+1); 
    for(int j=0; j < id_size; ++j)
      id_vec.emplace_back(id(j));
    for(int j=0; j < MMI.size(); ++j)
      TC[id_vec][j+4] = MMI[j];
    if(id_size == 3){ //max length
      orderedMMI[ MMI[1] ].emplace_back(id_vec);
    }
  }
  
  // get all value of MMI.maxNormalized
  NumericVector sorted_MMI(orderedMMI.size());
  int temp = 0;
  for(std::map<float, std::vector<std::vector<int> > >::iterator it= orderedMMI.begin(); it!=orderedMMI.end(); ++it){
    sorted_MMI[temp] = it->first;
    temp++;
  }
  NumericVector percentiles = NumericVector::create(1); //100%
  NumericVector percent = percentile_rcpp_sorted(sorted_MMI, percentiles);

  // stop iteratioin if it->first is larger than the 1 percentile value
  std::map<float, std::vector<std::vector<int> > > orderedTC; //sorting
  int size_orderedTC = 0;
  for (std::map<float,std::vector<std::vector<int> > >::iterator it= orderedMMI.begin(); it!=orderedMMI.end(); ++it){
    if(it->first <= percent[0]){
      for (std::vector<std::vector<int> >::iterator it2= it->second.begin(); it2!=it->second.end(); ++it2){
	std::vector<int> id_vec;
	id_vec.reserve(20); //just arbitrary reserving 
	for(int ii : (*it2)){
	  // use the value as the key for TC
	  id_vec.emplace_back(ii);
	}
	// auto-sort with map
	orderedTC[ -TC[id_vec][3] ].emplace_back(id_vec);  // -TC.maxNormalized
	size_orderedTC++;
      }
    } else{
      break;
    }
  }
  
  // Select top N TC.maxNormalized combinations
  std::map<float, std::vector<std::vector<int> > > newCOM_next0;  
  float current_TC = 0;
  int n = 0;
  int size_stored_ids = 0;  // used for output matrix row number  
  //int topN_TC = round(orderedTC.size() * 0.01);
  int topN_TC = round(orderedTC.size() * 1);
  if(topN_TC < 3) 
    topN_TC = 2;
//if(topN_TC >20) 
//  topN_TC = 20;  
  for (std::map<float, std::vector<std::vector<int> > >::iterator it= orderedTC.begin(); it!=orderedTC.end(); ++it){
    // only top N combis
    if(n > topN_TC)
      break;
    if(it->first == current_TC){
      for (std::vector<std::vector<int> >::iterator it2= it->second.begin(); it2!=it->second.end(); ++it2){
	// add to final_set	  	  
	std::vector<int> id_vec;
	id_vec.reserve(20); //just arbitrary reserving 
	for(int ii : (*it2)){
	  id_vec.emplace_back(ii);
	}
	newCOM_next0[it->first].emplace_back(id_vec);
	size_stored_ids++;
      }
      continue; // do not update n++	
    }else{
      for (std::vector<std::vector<int> >::iterator it2= it->second.begin(); it2!=it->second.end(); ++it2){
	// add to final_set	  	  
	std::vector<int> id_vec;
	id_vec.reserve(20); //just arbitrary reserving 
	for(int ii : (*it2)){
	  id_vec.emplace_back(ii);
	}
	newCOM_next0[it->first].emplace_back(id_vec);
	size_stored_ids++;
      }
      current_TC = it->first;
      n++;
    }
  }
  
  // Convert into matrix
  NumericMatrix final_set(size_stored_ids, 3);
  std::fill(final_set.begin(), final_set.end(), std::nan("1"));
  int current_row = 0;
  for (std::map<float, std::vector<std::vector<int> > >::iterator it= newCOM_next0.begin(); it!=newCOM_next0.end(); ++it){      
    for (std::vector<std::vector<int> >::iterator it2= it->second.begin(); it2!=it->second.end(); ++it2){
      int tmp=0;
      for(int ii : (*it2)){
	final_set(current_row,tmp) = ii;
	tmp++;
      }
      current_row++;
    }
  }
  
  ///////////////////////////////////////////
  // repeat until no more synergy increase //
  ///////////////////////////////////////////  
  NumericVector gen = myseq(0, GEN.size()-1);  // GEN size is the same as the total number of indices
  for(int ii=0; ii < max_add; ++ii){
    std::cout << ">" << ii << "\n";
    nr = final_set.nrow();
    nl = final_set.ncol();    
    // remove rowVector elements from gen and make the added matrix
    int nr_newCOM2 = (GEN.size() - nl) * nr;  
    NumericMatrix newCOM2(nr_newCOM2, nl+1);
    for(int i=0; i < nr; ++i){
      NumericVector gen_copy(gen.size());
      std::copy(std::begin(gen), std::end(gen), std::begin(gen_copy));      
      for(int j=0; j < nl; ++j){	
	gen_copy(final_set(i,j)) = std::nan("1");
      }
      // remove NA elements
      gen_copy = naomitInt(gen_copy);
      
      // add each element to the row
      for(int j_ind=0; j_ind < gen_copy.size(); ++j_ind){
	int j = (i * gen_copy.size()) + j_ind;
	for(int k=0; k < nl; ++k){
	  newCOM2(j, k) = final_set(i,k);
	}
	newCOM2(j, nl) = gen_copy(j_ind);
	NumericVector temp = newCOM2(j,_);
	newCOM2(j,_) = temp;
      }
    }
    
    // define the matrix size and initialize with nan  
    std::unordered_set<std::vector<int> > IDs;        
    int ncombi = newCOM2.nrow();
    int nco = newCOM2.ncol();
    IDs.reserve(ncombi*1000); //1000 is selecting all combis among 10
    for(int l=0; l < ncombi; ++l){
      NumericVector id = newCOM2(l,_);
      int n = id.size();
      for(int k=3; k <= nco; ++k){
	std::vector<bool> v(n);
	std::fill(v.end() - k, v.end(), true);        
	do { 
	  std::vector<int> id_vec;
	  id_vec.reserve(n+1); 
	  for (int i=0; i < n; ++i) {
	    if (v[i]) { 
	      id_vec.emplace_back(id(i));	 
	    }
	  }
	  sort(id_vec.begin(), id_vec.end() );
	  // check if already exist in TC	
	  if (TC.find(id_vec) == TC.end()){	
	    IDs.insert(id_vec);
	  }
	} while (std::next_permutation(v.begin(), v.end()));
      }
    }  
    NumericMatrix nCOM(IDs.size(), nl+1);
    std::fill(nCOM.begin(), nCOM.end(), std::nan("1") );
    int t1 = 0;
    int t2 = 0;
    for(std::vector<int> vec : IDs){
      int t2 = 0;
      for (std::vector<int >::iterator it= vec.begin(); it!=vec.end(); ++it){
	nCOM(t1,t2) = *it;
	t2++;
      }    
      t1++;
    }

    // Overwrite newCOM2
    newCOM2 = nCOM;
    nr = newCOM2.nrow();
    nl = newCOM2.ncol();
  
    // Check if subsets already exist in TC or not
    // if exists, continue and remove from newCOM 
    for(int i=0; i < nr; ++i){
      NumericVector id = naomitInt(newCOM2(i,_));
      NumericMatrix y3_sub(y3.nrow(),id.size());
      NumericMatrix H_sub(id.size(),H.nrow());
      for(int j=0; j < id.size(); ++j){
	y3_sub(_,j) = y3(_,id(j));
	H_sub(j,_) = H(id(j),_);
      }
    
      // compute JH and TC
      std::array<float, 10> jh; 
      NumericVector temp = computeJHandTC(y3_sub, H_sub);
      for(int j=0; j < temp.size(); ++j){
	jh[j] = temp[j];
      }      
      std::vector<int> id_vec;
      id_vec.reserve(id.size()+1); 
      for(int j=0; j < id.size(); ++j)
	id_vec.emplace_back(id(j));
      TC.insert(std::make_pair(id_vec, jh));
    }
    
    // Compute MMI, II
    std::map<float, std::vector<std::vector<int> > > orderedMMI2; //for sorting
    std::map<float, std::vector<std::vector<int> > > orderedTC2; //for sorting
    int size_stored_ids = 0;
    for(int i=0; i < nr; ++i){
      NumericVector id = naomitInt(newCOM2(i,_)); 
      int id_size = id.size();
      NumericVector MMI(4);
      
      // MMI, MMI.maxNormalized, II, II.maxNormalized
      // Add single entropies
      for(int j=0; j < id_size; ++j){
	MMI[0] += H(id(j),0); //MMI
	MMI[1] += H(id(j),2); //MMI.maxNormalized
      }
      MMI[2] = -(pow(-1, (id_size-1))) * MMI[0]; //II
      MMI[3] = -(pow(-1, (id_size-1))) * MMI[1]; //II.maxNomalized
	       
      // compute MMI and II
      int n = id.size();
      std::vector<std::vector<int> > seeds;
      seeds.reserve(500);
      for(int k=2; k <= n; ++k){
	std::vector<bool> v(n);
	std::fill(v.end() - k, v.end(), true);
	do {
	  std::vector<int> id_vec;
	  id_vec.reserve(n+1);
	  for (int i=0; i < n; ++i) {
	    if (v[i]) { 
	      id_vec.emplace_back(id(i));
	    }	   
	  }	  
	  // 1. MMI
	  MMI[0] += -(pow(-1, k)) * TC[id_vec][0];
	  MMI[1] += -(pow(-1, k)) * TC[id_vec][1];
	  // 2. II
	  MMI[2] += -(pow(-1, (id_size-k))) * TC[id_vec][0];
	  MMI[3] += -(pow(-1, (id_size-k))) * TC[id_vec][1];
	  
	  if(k==(n-1)){ // save when only 1 variable less (used for comparions)
	    seeds.emplace_back(id_vec);
	  }
	} while (std::next_permutation(v.begin(), v.end()));
      }
      std::vector<int> id_vec; id_vec.reserve(id_size+1);
      for(int j=0; j < id_size; ++j)
	id_vec.emplace_back(id(j));
      for(int j=0; j < MMI.size(); ++j)  // adding to TC
	TC[id_vec][j+4] = MMI[j];      
      if(id_size == nl){ //max length
	if(MMI[0] < 0){ // make sure synergy
	  // for each seed_vec, compute the difference and if not, continue
	  bool all_negative = true;

	  // Option 1: all MMI.maxNormalized difference is negative (original version)
	  bool delta = true;
	  delta = false;
	  if(delta){
	  for (std::vector<std::vector<int> >::iterator it= seeds.begin(); it!=seeds.end(); ++it){
	    if( (MMI[1] - TC[*it][5]) > cutoff_synergy ){ // new - previous < 0.05 (if more negative)
	      all_negative = false;
	      break;
	    }
	  }
	  }
	  if(all_negative){
	    orderedMMI2[MMI[1]].emplace_back(id_vec);
	    orderedTC2[ -TC[id_vec][3] ].emplace_back(id_vec);	  
	    size_stored_ids++;
	  }
	}
      }
    }

    // compute z-score among all entries and remove those below the z-score
    /*    float mean = std::accumulate(mmi1.begin(), mmi1.end(), 0.0) / mmi1.size();
    std::vector<float> diff(mmi1.size());
    std::transform(mmi1.begin(), mmi1.end(), diff.begin(), std::bind2nd(std::minus<float>(), mean));
    float sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    float stdev = std::sqrt(sq_sum / mmi1.size());
    float zcutoff1 = 100;
    //float zcutoff1 = -1.96;
    for (std::map<float, std::vector<std::vector<int> > >::iterator it= orderedMMI2.begin(); it!=orderedMMI2.end(); ++it){
      float z = (it->first - mean) / stdev;
      if(z >= zcutoff1){
	for (std::vector<std::vector<int> >::iterator it2= it->second.begin(); it2!=it->second.end(); ++it2){
	  std::vector<int> id_vec;
	  id_vec.reserve(20); //just arbitrary reserving
	  for(int ii : (*it2)){
	    id_vec.emplace_back(ii);
	  }
	  orderedTC2.erase( -TC[ id_vec ][3] ); // remove all id_vec 
	}
      }
      } */    
    /*    // Option 2: top N MMI
    int count1 = 1;
    int top_N_MMI1 = 1;
    for (std::map<float, std::vector<std::vector<int> > >::iterator it= orderedMMI2.begin(); it!=orderedMMI2.end(); ++it){
      if(count1 > top_N_MMI1){
	std::cout << it->first << "\n";
	for (std::vector<std::vector<int> >::iterator it2= it->second.begin(); it2!=it->second.end(); ++it2){
	  std::vector<int> id_vec;
	  id_vec.reserve(20); //just arbitrary reserving
	  for(int ii : (*it2)){
	    id_vec.emplace_back(ii);
	  }
	  orderedTC2.erase( -TC[ id_vec ][3] ); // remove all id_vec 
	}
      }
      count1++;
      } */   
    
    
    std::map<float, std::vector<std::vector<int> > > newCOM_next;  
    float current_TC = 0;
    topN_TC_heuristic = round(orderedTC2.size() * 1); //exhaustive
    int n = 0;
    size_stored_ids = 0;  // used for output matrix row number
    for (std::map<float, std::vector<std::vector<int> > >::iterator it= orderedTC2.begin(); it!=orderedTC2.end(); ++it){
      if(n > topN_TC_heuristic)
	break;
      if(it->first == current_TC){
	for (std::vector<std::vector<int> >::iterator it2= it->second.begin(); it2!=it->second.end(); ++it2){
	  // add to final_set	  	  
	  std::vector<int> id_vec;
	  id_vec.reserve(20); //just arbitrary reserving
	  for(int ii : (*it2)){
	    id_vec.emplace_back(ii);
	  }
	  newCOM_next[it->first].emplace_back(id_vec);
	  size_stored_ids++;
	}
	continue; // do not update n++	
      }else{
	for (std::vector<std::vector<int> >::iterator it2= it->second.begin(); it2!=it->second.end(); ++it2){
	  // add to final_set	  	  
	  std::vector<int> id_vec;
	  id_vec.reserve(20); //just arbitrary reserving
	  for(int ii : (*it2)){
	    id_vec.emplace_back(ii);
	  }
	  newCOM_next[it->first].emplace_back(id_vec);
	  size_stored_ids++;
	}
	current_TC = it->first;
	n++;
      }
    }    
    if(size_stored_ids==0) // stop iterating
      break;    

    // Convert into matrix
    NumericMatrix newCOM_final(size_stored_ids, nl);
    std::fill(newCOM_final.begin(), newCOM_final.end(), std::nan("1"));
    int current_row = 0;
    for (std::map<float, std::vector<std::vector<int> > >::iterator it= newCOM_next.begin(); it!=newCOM_next.end(); ++it){      
      for (std::vector<std::vector<int> >::iterator it2= it->second.begin(); it2!=it->second.end(); ++it2){
	int tmp=0;
	for(int ii : (*it2)){
	  newCOM_final(current_row,tmp) = ii;
	  tmp++;
	}
	current_row++;
      }
    }
    final_set = newCOM_final; // pass to global variable
  }

  // Option 1: extract only top X combis
  bool t1 = true;
  //t1 = false;
  if(t1){
  std::map<float, std::vector<std::vector<int> > > orderedMMI3; //for sorting
  for(int j=0; j < final_set.nrow(); ++j){
    std::vector<int> id_vec;    
    for(int k=0; k < final_set.ncol(); ++k){
      id_vec.emplace_back(final_set(j,k));	
    }
    //std::cout << TC[id_vec][3] << " " << TC[id_vec][5] << "\n";
    orderedMMI3[ TC[id_vec][5] ].emplace_back(id_vec);
  }
  int topX = 10;
  if(final_set.nrow() < topX )
    topX = final_set.nrow();
  NumericMatrix final_set2(topX, final_set.ncol());
  std::fill(final_set2.begin(), final_set2.end(), std::nan("1"));
  current_row = 0;
  int count = 1;
  for (std::map<float, std::vector<std::vector<int> > >::iterator it= orderedMMI3.begin(); it!=orderedMMI3.end(); ++it){
    if(count <= topX){
      for (std::vector<std::vector<int> >::iterator it2= it->second.begin(); it2!=it->second.end(); ++it2){
	if(count <= topX){
	  std::vector<int> id_vec;
	  id_vec.reserve(20); //just arbitrary reserving
	  int tmp=0;
	  for(int ii : (*it2)){
	    final_set2(current_row,tmp) = ii;
	    tmp++;
	  }
	  current_row++;
	  count++;
	}
      }
    }
  }
  final_set = final_set2;
  }

  
  // make output lists
  List tempList(final_set.nrow());
  List tempTC(final_set.nrow());
  for(int j=0; j < final_set.nrow(); ++j){
    StringVector st(final_set.ncol());
    std::vector<int> id_vec;    
    for(int k=0; k < final_set.ncol(); ++k){
      st[k] = GEN(final_set(j,k));
      id_vec.emplace_back(final_set(j,k));	
    }
    tempList[j] = st;
    //NumericVector tc(TC[id_vec].size());    
    //std::copy(std::begin(TC[id_vec]), std::end(TC[id_vec]), std::begin(tc));
    NumericVector tc(2);
    tc(0) = TC[id_vec][3]; //TC.maxNormalized
    tc(1) = TC[id_vec][5]; //MMI.maxNormalized
    tempTC[j] = tc;
  }
  return Rcpp::List::create(Rcpp::Named("final_set") = final_set, Rcpp::Named("bestCombis") = tempList, Rcpp::Named("TCs") = tempTC, Rcpp::Named("topN_TC") = topN_TC);    
}


// [[Rcpp::export]]
List compute_MMI_finalCombis2(NumericMatrix& newCOM, NumericMatrix& y3, NumericMatrix& H, StringVector& GEN, int& nJSD)
{
  std::cout << "starting c++ compute_MMI_finalCombis..." << "\n";  
  // use std::vector<int, newCOM[i,]> as the key
  std::unordered_map<std::vector<int>, std::array<float, 10> > TC;  
  NumericVector gen = myseq(0, GEN.size()-1);
  int ncombi = newCOM.nrow();
  int nco = newCOM.ncol();
  TC.reserve(ncombi+1000000); //initial combis + 10 variable combinations  

  // add all subsets (from 2)
  // define the matrix size and initialize with nan
  std::vector<std::vector<int> > input_vecs;
  input_vecs.reserve(ncombi+1);
  std::unordered_set<std::vector<int> > IDs;        
  IDs.reserve(ncombi*1000); //1000 is selecting all combis among 10
  for(int l=0; l < ncombi; ++l){
    NumericVector id = naomitInt(newCOM(l,_));
    std::vector<int> id_vec;
    id_vec.reserve(id.size()+1); 
    for (int i=0; i < id.size(); ++i) {
      id_vec.emplace_back(id(i));	 
    }
    input_vecs.emplace_back(id_vec);
    
    int n = id.size();
    for(int k=2; k <= n; ++k){
      std::vector<bool> v(n);
      std::fill(v.end() - k, v.end(), true);        
      do { 
	std::vector<int> id_vec;
	id_vec.reserve(n+1); 
	for (int i=0; i < n; ++i) {
	  if (v[i]) { 
	    id_vec.emplace_back(id(i));	 
	  }
	}
	sort(id_vec.begin(), id_vec.end() );
	// check if already exist in TC	
	if (TC.find(id_vec) == TC.end()){	
	  IDs.insert(id_vec);
	}
      } while (std::next_permutation(v.begin(), v.end()));
    }
  }  
  NumericMatrix nCOM(IDs.size(), nco);
  std::fill(nCOM.begin(), nCOM.end(), std::nan("1") );
  int t1 = 0;
  int t2 = 0;
  for(std::vector<int> vec : IDs){
    int t2 = 0;
    for (std::vector<int >::iterator it= vec.begin(); it!=vec.end(); ++it){
      nCOM(t1,t2) = *it;
      t2++;
    }    
    t1++;
  }
  // Overwrite newCOM2
  newCOM = nCOM;
    
  int nr = newCOM.nrow();
  int nl = newCOM.ncol(); 
  for(int i=0; i < nr; ++i){    
    NumericVector id = naomitInt(newCOM(i,_));    
    int id_size = id.size();
    NumericMatrix y3_sub(y3.nrow(),id.size());
    NumericMatrix H_sub(id_size,H.nrow());
    for(int j=0; j < id_size; ++j){
      y3_sub(_,j) = y3(_,id(j));
      H_sub(j,_) = H(id(j),_);
    }
    
    // compute JH and TC
    std::array<float, 10> jh; 
    NumericVector temp = computeJHandTC(y3_sub, H_sub);
    for(int j=0; j < temp.size(); ++j){
      jh[j] = temp[j];
    }
    std::vector<int> id_vec;
    int vr = id_size+1;
    id_vec.reserve(vr); 
    for(int j=0; j < id_size; ++j)
      id_vec.emplace_back(id(j));
    TC.insert(std::make_pair(id_vec, jh));
  }
  
  // Compute MMI, II
  for(int i=0; i < nr; ++i){
    NumericVector id = naomitInt(newCOM(i,_)); 
    int id_size = id.size();
    NumericVector MMI(4);
    
    // MMI, MMI.maxNormalized, II, II.maxNormalized
    // Add single entropies
    for(int j=0; j < id_size; ++j){
      MMI[0] += H(id(j),0); //MMI
      MMI[1] += H(id(j),2); //MMI.maxNormalized
    }
    MMI[2] = -(pow(-1, (id_size-1))) * MMI[0]; //II
    MMI[3] = -(pow(-1, (id_size-1))) * MMI[1]; //II.maxNomalized

    // compute MMI and II
    int n = id.size();  
    for(int k=2; k <= n; ++k){
      std::vector<bool> v(n);
      std::fill(v.end() - k, v.end(), true);
      do {
	std::vector<int> id_vec;
	id_vec.reserve(n+1);
	for (int i=0; i < n; ++i) {
	  if (v[i]) { 
	    id_vec.emplace_back(id(i));
	  }	   
	}	  
	// 1. MMI
	MMI[0] += -(pow(-1, k)) * TC[id_vec][0];
	MMI[1] += -(pow(-1, k)) * TC[id_vec][1];
	// 2. II
	MMI[2] += -(pow(-1, (id_size-k))) * TC[id_vec][0];
	MMI[3] += -(pow(-1, (id_size-k))) * TC[id_vec][1];
      } while (std::next_permutation(v.begin(), v.end()));
    }  
    std::vector<int> id_vec;
    id_vec.reserve(id_size+1);
    for(int j=0; j < id_size; ++j)
      id_vec.emplace_back(id(j));
    for(int j=0; j < MMI.size(); ++j)  // adding to TC
      TC[id_vec][j+4] = MMI[j];
  }

  // Check if each combi is more synergistic than before
  // for each input vector, compute the difference and if not, continue
  std::map<float, std::vector<std::vector<int> > > orderedMMI; 
  std::map<float, std::vector<std::vector<int> > > orderedTC;
  int size_stored_ids = 0;
  for (std::vector<std::vector<int> >::iterator it= input_vecs.begin(); it!=input_vecs.end(); ++it){
    bool all_negative = true;
    int vec_size = (*it).size();    
    int k = vec_size -1; // one variable less
    std::vector<bool> v(vec_size);
    std::fill(v.end() - k, v.end(), true);
    do {
      std::vector<int> seed_vec;
      seed_vec.reserve(vec_size+1);
      for (int i=0; i < vec_size; ++i) {
	if (v[i]) { 
	  seed_vec.emplace_back( (*it)[i] );	    
	}
      }
      if( TC[*it][4] > 0 ){ //MMI check
	all_negative = false;	
	break;
      }
    } while (std::next_permutation(v.begin(), v.end()));
    
    if(all_negative){
      orderedTC[ -TC[(*it)][3] ].emplace_back( *it );
      size_stored_ids++;
    }
  } 
  
  // Convert into matrix
  NumericMatrix final_set(size_stored_ids, nco);
  std::fill(final_set.begin(), final_set.end(), std::nan("1"));
  int current_row = 0;
  for (std::map<float, std::vector<std::vector<int> > >::iterator it=orderedTC.begin(); it!=orderedTC.end(); ++it){
    for (std::vector<std::vector<int> >::iterator it2= it->second.begin(); it2!=it->second.end(); ++it2){
      int tmp=0;
      for(int ii : (*it2)){
	final_set(current_row,tmp) = ii;
	tmp++;
      }
      current_row++;
    }
  }
  
  // make output lists
  List tempList(final_set.nrow());
  //  List tempTC(final_set.nrow());
  NumericMatrix tempTC(final_set.nrow(),2); //normalized MMI and TC
  for(int j=0; j < final_set.nrow(); ++j){
    NumericVector id = naomitInt(final_set(j,_));
    StringVector st(id.size());
    std::vector<int> id_vec;
    for(int k=0; k < id.size(); ++k){
      st[k] = GEN(id(k));
      id_vec.emplace_back(id(k));	
    }
    tempList[j] = st;
    //NumericVector tc(TC[id_vec].size());
    //std::copy(std::begin(TC[id_vec]), std::end(TC[id_vec]), std::begin(tc));
    NumericVector tc(2);
    //std::cout << TC[id_vec][2] << " " << TC[id_vec][3] << " " << TC[id_vec][4] << " " << TC[id_vec][5] << "\n";
    tc(0) = TC[id_vec][3];
    tc(1) = TC[id_vec][5];    
    tempTC(j,_) = tc;    
  }

  // only the nonJSD final_set
  //int nJSD = 5;
  NumericMatrix fset;
  if(final_set.nrow() > 0 && (final_set.ncol()==nco)){
    fset = final_set(_, Range(nJSD,final_set.ncol()-1));
  }else{ // do not return any value
    tempTC = fset;
    //fset = final_set;
  }
  
  return Rcpp::List::create(Rcpp::Named("final_set") = fset, Rcpp::Named("TCs") = tempTC);
}


