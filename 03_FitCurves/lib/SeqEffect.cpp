// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <zlib.h>
#include <chrono>
#include <algorithm>

// for fast loading of fasta/q files
#include "kseq.h"

// fast hash map
#include <gtl/include/gtl/phmap.hpp>

using namespace Rcpp;

// create new kseq_t type with file type gzFile 
// and read function gzread.
KSEQ_INIT(gzFile, gzread)

// create alias for gzFile as file_t
typedef gzFile file_t;
// alias for time stamps
typedef std::chrono::duration<float> float_seconds;
// mappings
typedef gtl::flat_hash_map<unsigned int, std::pair<int, std::string>> gtl_hash_kmers;

/**
 * Struct to store first lexicographically occurring kmers.
*/
struct KmerTable {
  std::vector<std::string> fwd_kmers;
  std::vector<std::string> rc_kmers;
};

/**
 * Struct for storing hash kmers, fwd and rc kmers.
*/
struct KmerEncoding {
    gtl_hash_kmers kmer_encodings;
    std::vector<std::string> k_string_vec;
    std::vector<std::string> k_rc_string_vec;
};

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

/**
 * Hashes a kmer string.
 * @param vec string to hash.
 * @param k encoding the string into size k.
 * @param i current index of the iterations.
 * @return hashed string
*/
inline unsigned int hash_func(const std::string &vec, int &k, int i = 0){
    // reset hash
    unsigned int hash = 5381;

    // Update the hash value for the current k-mer
    for(int j = 0; j < k; j++){
      hash = ((hash << 5) + hash) + base_to_encoding[(unsigned char)vec[i+j]];
    }

    return hash;
}

std::string kmers(int i, const int &k){
    static const char* DNA = "ATGC";
    if(k == 0){
        return std::string();
    }
    return kmers(i/4, k-1) + DNA[i%4];
}

/**
 * Generate all possible kmer combinations.
 * @param k size of kmer
 * @return all kmers.
*/
// [[Rcpp::export]]
std::vector<std::string> generate_kmers(const int &k){
  std::vector<std::string> all_kmers(1 << (2*k));
  
  for(int i = 0; i < (1 << (2*k)); i++){
    all_kmers[i] = kmers(i, k);
  }

  // sort kmer vector alphabetically
  std::sort(all_kmers.begin(), all_kmers.end());

  return all_kmers;
}

/**
 * Get reverse complement of a DNA string.
 * @param sequence kmer broken after aligning read with de novo sequence/
 * @return reverse_seq reverse complement of sequence.
*/
std::string reverse_complement(const std::string &sequence){
    // reverse complement map
    gtl::flat_hash_map<char, char> complement = {
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
 * Generates a table of forward and reverse complement kmers
 * of first lexicographically occuring kmers.
 * @param k size of kmer.
 * @return kmer table struct.
*/
KmerTable get_lexicographic_kmertable(const int &k){
  // init structure
  KmerTable kmer_table;

  // generate kmers
  kmer_table.fwd_kmers = generate_kmers(k);
  kmer_table.rc_kmers.resize(kmer_table.fwd_kmers.size());
  for(int i = 0; i < kmer_table.fwd_kmers.size(); i++){
    kmer_table.rc_kmers[i] = reverse_complement(kmer_table.fwd_kmers[i]);
  }

  // only retain first lexicographically occurring kmer
  for(int i = 0; i < kmer_table.fwd_kmers.size(); i++){
    bool is_fwd_first = (kmer_table.fwd_kmers[i] <= kmer_table.rc_kmers[i]);
    if(!is_fwd_first){
        // Remove the corresponding element from both fwd_kmers and rc_kmers
        kmer_table.fwd_kmers.erase(kmer_table.fwd_kmers.begin()+i);
        kmer_table.rc_kmers.erase(kmer_table.rc_kmers.begin()+i);
        
        // Decrement i by 1 to avoid skipping the next element
        i--;
    }
  }

  return kmer_table;
}

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
 * Encode genome sequence.
 * @param ref_seq reference sequence.
 * @param hash_kmers encoded kmers.
 * @param k size of kmer.
 * @return encoded sequence.
*/
std::vector<int> encode_refseq(
  const std::string &ref_seq, 
  gtl_hash_kmers &hash_kmers, 
  int k){

  // Reserve space for the number of kmers in the reference sequence
  int num_kmers = ref_seq.length()-k+1;
  std::vector<int> ref_encodings(num_kmers);

  // Loop through the reference sequence
  for(int i = 0; i <= ref_seq.length()-k; i++){
    unsigned int hash = hash_func(ref_seq, k, i);

    // find index of encoded kmer in hash_kmers
    int ind = hash_kmers[hash].first;
    
    // Add the k-mer encoding and its starting position to the map
    ref_encodings[i] = ind;
  }
  return ref_encodings;
}

/**
 * Encode each kmer sequence.
 * @param all_kmers vector of kmers.
 * @param k size of kmer.
 * @return encoded kmers.
*/
KmerEncoding encode_kmers(const std::vector<std::string> &all_kmers, int k){
    KmerEncoding kmer_encoding;

    // Loop through each kmer
    for(int i = 0; i < all_kmers.size(); i++){
      unsigned int hash = hash_func(all_kmers[i], k);

      // key: k-mer encoding. value: [init index at 0, kmer]
      kmer_encoding.kmer_encodings[hash] = std::make_pair(0, all_kmers[i]);
    }

    // get index of encoded kmers
    kmer_encoding.k_string_vec.resize(kmer_encoding.kmer_encodings.size());
    kmer_encoding.k_rc_string_vec.resize(kmer_encoding.kmer_encodings.size());
    int ind = 0;
    for(auto &kv : kmer_encoding.kmer_encodings){
        // get kmer
        std::string hash_kmer = kv.second.second;

        // and push into vector
        kmer_encoding.k_string_vec[ind] = hash_kmer;
        kmer_encoding.k_rc_string_vec[ind] = reverse_complement(hash_kmer);

        // update value by its index in map
        kv.second.first = ind;

        ind++;
    }

    return kmer_encoding;
}

/**
 * Calulate rmsd between two adjacent positions.
 * @param a normalised kmer frequencies.
 * @param b normalised kmer frequencies.
 * @return root-mean-square deviation between a and b.
*/
double rmsd(const std::vector<double> &a, 
            const std::vector<double> &b){
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

/**
 * Calculate kmer frequencies at each shifted breakpoint position.
 * @param bp_pos all breakpoints at single nucleotide resolution.
 * @param filename file name of the reference genome.
 * @param k size of breakage kmer.
 * @param rmsd_range range to calculate sequence context influence.
 * @return rmsd between adjacent breakpoint positions.
*/
// [[Rcpp::export]]
std::vector<double> calc_kmer_freq(std::vector<int> &bp_pos,
                                   const std::string &filename,
                                   const int &k,
                                   const std::vector<int> &rmsd_range){

    // read compressed fasta file
    std::vector<std::string> ref_seq_vec = read_compressed_fasta(filename);
    std::string ref_seq = ref_seq_vec[0];
    ref_seq_vec.clear();

    // generate kmers
    const std::vector<std::string> all_kmers = generate_kmers(k);
    
    // generate hash maps of kmer encodings
    KmerEncoding kmer_encoding = encode_kmers(all_kmers, k);
    gtl_hash_kmers hash_kmers = kmer_encoding.kmer_encodings;
    std::vector<std::string> k_string_vec = kmer_encoding.k_string_vec;
    std::vector<std::string> k_rc_string_vec = kmer_encoding.k_rc_string_vec;

    // generate hash maps of kmer encodings for reference sequence
    std::vector<int> hash_refseq = encode_refseq(ref_seq, hash_kmers, k);
    ref_seq.clear();
    
    // hash map of fwd and rc kmers
    KmerTable kmer_table = get_lexicographic_kmertable(k);
    gtl::flat_hash_map<std::string, std::pair<std::string, int>> kmer_matrix;
    for(int i = 0; i < kmer_table.fwd_kmers.size(); i++){
        kmer_matrix[kmer_table.fwd_kmers[i]] = std::make_pair(kmer_table.rc_kmers[i], 0);
    }

    // init average relative kmer frequency
    int rmsd_len = rmsd_range.size();
    gtl::flat_hash_map<int, std::vector<double>> kmer_mean_all;
    for(int i = 0; i < rmsd_len; i++){
      kmer_mean_all[i] = std::vector<double>(kmer_table.fwd_kmers.size(), 0);
    }

    // init kmer count vector
    std::vector<int> k_count_vec(hash_kmers.size(), 0);

    for(int outer_ind = 0; outer_ind < rmsd_len; outer_ind++){

      for(int j = 0; j < bp_pos.size(); j++){
        // expand breakpoint positions into kmers
        int start_pos = (int)bp_pos[j]+(int)rmsd_range[outer_ind]-(k/2);

        // find encoded string in hash_refseq
        int encoded_val_ind = hash_refseq[start_pos];

        // update value count in k_count_vec
        k_count_vec[encoded_val_ind]++;
      }

      // get relative kmer frequency count
      for(int i = 0; i < k_count_vec.size(); i++){
        // loop over all kmers in encoded kmer map
        const std::string &kmer_str = k_string_vec[i];
        const std::string &rc_kmer = k_rc_string_vec[i];
        const int count = k_count_vec[i];

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
      for(int i = 0; i < k_count_vec.size(); i++){
        k_count_vec[i] = 0;
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