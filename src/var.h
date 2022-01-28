#ifndef VAR_H
#define VAR_H

#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <htslib/sam.h>


namespace lorax
{

  struct BiallelicVariant {
    int32_t pos;
    char ref;
    char alt;
    bool hap;
    uint32_t rsup;
    uint32_t asup;

    explicit BiallelicVariant(int32_t const p) : pos(p), ref('#'), alt('#'), hap(0), rsup(0), asup(0) {}
    BiallelicVariant(int32_t const p, char const r, char const a, bool const h) : pos(p), ref(r), alt(a), hap(h), rsup(0), asup(0) {}
  };


  template<typename TRecord>
  struct SortVariants : public std::binary_function<TRecord, TRecord, bool> {
    inline bool operator()(TRecord const& s1, TRecord const& s2) const {
      return s1.pos < s2.pos;
    }
  };

  template<typename TVariants>
  inline bool
  _loadVariants(htsFile* ifile, hts_idx_t* bcfidx, bcf_hdr_t* hdr, std::string const& sample, std::string const& chrom, TVariants& pV) {
    typedef typename TVariants::value_type TVariant;
    
    int32_t sampleIndex = -1;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
      if (hdr->samples[i] == sample) sampleIndex = i;
    if (sampleIndex < 0) return false;
        
    // Genotypes
    int ngt = 0;
    int32_t* gt = NULL;

    // Collect het. bi-allelic variants for this chromosome
    int32_t chrid = bcf_hdr_name2id(hdr, chrom.c_str());
    int32_t lastpos = -1;
    if (chrid < 0) return false;
    hts_itr_t* itervcf = bcf_itr_querys(bcfidx, hdr, chrom.c_str());
    if (itervcf != NULL) {
      bcf1_t* rec = bcf_init1();
      while (bcf_itr_next(ifile, itervcf, rec) >= 0) {
	// Only bi-allelic variants
	if (rec->n_allele == 2) {
	  bcf_unpack(rec, BCF_UN_ALL);
	  bcf_get_genotypes(hdr, rec, &gt, &ngt);
	  if ((bcf_gt_allele(gt[sampleIndex*2]) != -1) && (bcf_gt_allele(gt[sampleIndex*2 + 1]) != -1) && (!bcf_gt_is_missing(gt[sampleIndex*2])) && (!bcf_gt_is_missing(gt[sampleIndex*2 + 1]))) {
	    int gt_type = bcf_gt_allele(gt[sampleIndex*2]) + bcf_gt_allele(gt[sampleIndex*2 + 1]);
	    if (gt_type == 1) {
	      if (rec->pos != lastpos) {
		// Only one variant per position
		std::string ref = std::string(rec->d.allele[0]);
		std::string alt = std::string(rec->d.allele[1]);
		if ((ref.size() == 1) && (alt.size() == 1)) {
		  // Only SNPs
		  pV.push_back(TVariant(rec->pos, ref[0], alt[0], bcf_gt_allele(gt[sampleIndex*2])));
		  lastpos = rec->pos;
		}
	      }
	    }
	  }
	}
      }
      bcf_destroy(rec);
      hts_itr_destroy(itervcf);
    }
    if (gt != NULL) free(gt);
    return true;
  }

  
  template<typename TVariants>
  inline bool
  _loadVariants(std::string const& sample, std::string const& chrom, std::string const& bcffile, TVariants& pV) {
    // Load BCF file
    htsFile* ifile = bcf_open(bcffile.c_str(), "r");
    hts_idx_t* bcfidx = bcf_index_load(bcffile.c_str());
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);

    bool retVal = _loadVariants(ifile, bcfidx, hdr, sample, chrom, pV);
    
    // Close BCF
    bcf_hdr_destroy(hdr);
    hts_idx_destroy(bcfidx);
    bcf_close(ifile);
    
    return retVal;
  }

}

#endif
