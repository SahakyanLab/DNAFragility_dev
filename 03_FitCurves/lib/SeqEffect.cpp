// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <map>
#include <zlib.h>
#include <chrono>
#include <algorithm>
#include <future>

// for fast loading of fasta/q files
#include "kseq.h"

// fast hash map
#include <boost/unordered_map.hpp>

// // parallel for loop
// #include "/usr/local/Cellar/libomp/15.0.7/include/omp.h"

using namespace Rcpp;

// create new kseq_t type with file type gzFile 
// and read function gzread.
KSEQ_INIT(gzFile, gzread)

// create alias for gzFile as file_t
typedef gzFile file_t;
// alias for time stamps
typedef std::chrono::duration<float> float_seconds;
// mappings
typedef boost::unordered::unordered_map<unsigned int, std::pair<int, std::string>> boost_umap;

/**
 * Reads compressed fasta file.
 * @param f file type gzFile.
 * @param buf buffer.
 * @param len length of buffer.
 * @return number of bytes from file.
*/
int readfn(file_t f, char *buf, int len) {
    return gzread(f, buf, len);
}

/**
 * Reads compressed fasta file.
 * @param filename name and location of fasta.gz file.
 * @return List of fasta sequences.
*/
// [[Rcpp::export]]
std::vector<std::string> read_compressed_fasta(const std::string &filename) {
    std::vector<std::string> ref_seq;

    // open file in read mode and returns pointer to file
    file_t fp = gzopen(filename.c_str(), "r");

    // init new kseq_t type using fp as pointer
    kseq_t *seq = kseq_init(fp);

    // loop over reads in file
    while(kseq_read(seq) >= 0){
        ref_seq.push_back(seq->seq.s);
    }

    // clear files in memory
    kseq_destroy(seq);
    gzclose(fp);
    return ref_seq;
}

/**
 * Define the encoding for each base ATGC, which represents 
 * the encoding of each character in the ASCII character set.
 * 
 * The array is indexed by the ASCII value of each character, 
 * with the first 128 values covering the standard ASCII character set.
 * 
 * The encoding of each character is represented by a single byte, 
 * with the possible values being 0, 1, 2, or 3.
*/
constexpr char base_to_encoding[(unsigned char)128] = {
  /* 0-127 */
  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
  0, 0, 2, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0,
  0, 0, 0, 0,  3, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
  0, 2, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
};

// Initialize the hash value
unsigned int hash = 5381;

/**
 * Encode genome sequence.
 * @param ref_seq reference sequence.
 * @param k size of kmer.
 * @return encoded sequence.
*/
std::vector<int> encode_ref(const std::string &ref_seq, int k) {
  // Reserve space for the number of kmers in the reference sequence
  int num_kmers = ref_seq.length() - k + 1;
  std::vector<int> ref_encodings(num_kmers);

  // Loop through the reference sequence
  for(int i = 0; i <= ref_seq.length() - k; i++){
    
    // Update the hash value for the current k-mer
    for(int j = 0; j < k; j++){
      hash = ((hash << 5) + hash) + base_to_encoding[(unsigned char)ref_seq[i+j]];
    }
    
    // Add the k-mer encoding and its starting position to the map
    ref_encodings[i] = hash;
    
    // Update the hash value for the current k-mer
    hash = 5381;
  }
  return ref_encodings;
}

/**
 * Encode each kmer sequence.
 * @param results vector of kmers.
 * @param k size of kmer.
 * @return encoded kmers.
*/
boost_umap encode_kmers(const std::vector<std::string> &results, int k){
  boost_umap kmer_encodings;

  // Loop through each kmer
  for(int i = 0; i < results.size(); i++){
    
    // Update the hash value for the current k-mer
    for(int j = 0; j < k; j++){
      hash = ((hash << 5) + hash) + base_to_encoding[(unsigned char)results[i][j]];
    }

    // Add the k-mer encoding and its starting position to the map
    kmer_encodings[hash] = std::make_pair(0, results[i]);

    // Update the hash value for the current k-mer
    hash = 5381;
  }
  return kmer_encodings;
}

std::string kmers(int i, const int &kmer){
  static const char* DNA = "ATGC";
  if(kmer == 0){
    return std::string();
  }
  return kmers(i/4, kmer-1) + DNA[i%4];
}

/**
 * Generate all possible kmer combinations.
 * @param kmer size of kmer
 * @return all kmers.
*/
// [[Rcpp::export]]
std::vector<std::string> generate_kmers(const int &kmer){
  std::vector<std::string> results(1 << (2*kmer));
  
  for(int i = 0; i < (1 << (2*kmer)); i++){
    results[i] = kmers(i, kmer);
  }

  // sort kmer vector alphabetically
  std::sort(results.begin(), results.end());

  return results;
}

/**
 * Get reverse complement of a DNA string.
 * @param sequence kmer broken after aligning read with de novo sequence/
 * @return reverse_seq reverse complement of sequence.
*/
std::string reverse_complement(const std::string &sequence){
    // reverse complement map
    boost::unordered::unordered_map<char, char> complement = {
        {'A', 'T'}, 
        {'C', 'G'}, 
        {'G', 'C'}, 
        {'T', 'A'}
    };
    std::string reverse_seq(sequence.rbegin(), sequence.rend());
    for(char &c : reverse_seq){
        c = complement[c];
    }
    return reverse_seq;
}

/**
 * Calulate rmsd between two adjacent positions.
*/
double rmsd(const std::vector<double>& a, 
            const std::vector<double>& b){
    std::vector<double> diff(a.size());
    std::transform(
      a.begin(), a.end(), 
      b.begin(), diff.begin(),
      [](double x, double y){ 
        return x - y; 
    });

    double sum_sq = std::inner_product(
      diff.begin(), 
      diff.end(), 
      diff.begin(), 
      0.0
    );
    return std::sqrt(sum_sq/a.size());
}

// [[Rcpp::export]]
std::vector<double> calc_kmer_freq(std::vector<int> &bp_pos,
                                   const std::vector<std::string> &results,
                                   std::string &ref_seq,
                                   const int &kmer,
                                   const std::vector<std::string> &fwd_kmer_map,
                                   const std::vector<std::string> &rc_kmer_map,
                                   const int &num_threads,
                                   const std::vector<int> &rmsd_range){
    // generate hash maps of encodings
    boost_umap hash_kmers = encode_kmers(results, kmer);
    std::vector<int> hash_ref = encode_ref(ref_seq, kmer);
    ref_seq.clear();

    // hash map of fwd and rc kmers
    boost::unordered::unordered_map<std::string, std::pair<std::string, int>> kmer_matrix;
    for(int i = 0; i < fwd_kmer_map.size(); i++){
        kmer_matrix[fwd_kmer_map[i]] = std::make_pair(rc_kmer_map[i], 0);
    }

    // init average relative kmer frequency
    int rmsd_len = rmsd_range.size();
    std::map<int, std::vector<double>> kmer_mean_all;
    for(int i = 0; i < rmsd_len; i++){
      kmer_mean_all[i] = std::vector<double>(fwd_kmer_map.size(), 0);
    }

    Rcout << "Hashing done!" << "\n";

    // #pragma omp parallel for num_threads(4)    
    for(int outer_ind = 0; outer_ind < rmsd_len; outer_ind++){

      for(int j = 0; j < bp_pos.size(); j++){
        // expand breakpoint positions into kmers
        int start_pos = (int)bp_pos[j]+(int)rmsd_range[outer_ind]-(kmer/2);

        // find encoded string in hash_ref
        int encoded_val = hash_ref[start_pos];

        // update value count in hash_kmers
        hash_kmers[encoded_val].first += 1;
      }

      // get relative kmer frequency count
      for(const auto &kv : hash_kmers){
        // loop over all kmers in encoded kmer map
        const std::string& kmer_str = kv.second.second;
        const int count = kv.second.first;

        // get reverse complement of kmer
        const std::string rc_kmer = reverse_complement(kmer_str); 

        // check if fwd kmer
        auto fwd_it = kmer_matrix.find(kmer_str);
        if(fwd_it != kmer_matrix.end()){
          std::get<1>(fwd_it->second) += count;

          // check if palindrome
          if(kmer_str == rc_kmer){
            std::get<1>(fwd_it->second) += count;
          }
        } else {
          auto rc_it = kmer_matrix.find(rc_kmer);
          std::get<1>(rc_it->second) += count;
        }
      }

      // get total encode matching counts
      int total_encode_matches = 0.0;
      for(const auto &kv : kmer_matrix){
        total_encode_matches += kv.second.second;
      }

      std::vector<double> kmer_counts;
      kmer_counts.reserve(kmer_matrix.size());      
      for(auto &kv : kmer_matrix){
        double rel_count = (int)kv.second.second/(double)total_encode_matches;
        kmer_counts.push_back(rel_count);

        // reset counter
        kv.second.second = 0;
      }

      // reset encoded kmer counts
      for(auto &kv : hash_kmers){
        kv.second.first = 0;
      }

      // push vector value map
      kmer_mean_all[outer_ind] = kmer_counts;
    }

    // calculate rmsd between adjacent positions
    std::vector<double> rmsd_vec(kmer_mean_all.size()-1);
    for(int i = 0; i < kmer_mean_all.size()-1; i++){
      std::vector<double> a = kmer_mean_all[i];
      std::vector<double> b = kmer_mean_all[i+1];

      double rmsd_result = rmsd(a, b);
      rmsd_vec[i] = rmsd_result;
    }
    return rmsd_vec;
}