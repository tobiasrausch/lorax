#ifndef UTIL_H
#define UTIL_H

#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/tokenizer.hpp>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>

namespace lorax
{

  #ifndef INVALID
  #define INVALID 4294967295
  #endif


  inline void
  reverseComplement(std::string& sequence) {
    std::string rev = boost::to_upper_copy(std::string(sequence.rbegin(), sequence.rend()));
    std::size_t i = 0;
    for(std::string::iterator revIt = rev.begin(); revIt != rev.end(); ++revIt, ++i) {
      switch (*revIt) {
      case 'A': sequence[i]='T'; break;
      case 'C': sequence[i]='G'; break;
      case 'G': sequence[i]='C'; break;
      case 'T': sequence[i]='A'; break;
      case 'N': sequence[i]='N'; break;
      default: break;
      }
    }
  }


  template<typename TConfig>
  inline void
  createRepeatMotifs(TConfig& c, bool const mix) {
    // Fetch repeats
    typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
    boost::char_separator<char> sep(",");
    Tokenizer tokens(c.repeatStr, sep);
    typedef std::vector<std::string> TRepeatSet;
    TRepeatSet inputRepeats;
    for(Tokenizer::iterator tokIter = tokens.begin(); tokIter != tokens.end(); ++tokIter) {
      if (inputRepeats.empty()) c.replen = tokIter->size();
      if (c.replen == tokIter->size()) inputRepeats.push_back(*tokIter);
      else {
	std::cerr << "Input repeats must have the same length!" << std::endl;
      }
    }

    // Repeat period
    TRepeatSet repset = inputRepeats;
    // Compute one additional period
    for(uint32_t i = 0; i < c.period; ++i) {
      uint32_t stopsize = repset.size();
      for(uint32_t k = 0; k < stopsize; ++k) {
	std::string sourceStr = repset[k];
	for(uint32_t j = 0; j < inputRepeats.size(); ++j) {
	  if (k % inputRepeats.size() == j) repset[k] += inputRepeats[j];
	  else {
	    if (mix) repset.push_back(sourceStr + inputRepeats[j]);
	  }
	}
      }
    }

    // Subset to replen * period substrings
    for(uint32_t k = 0; k < repset.size(); ++k) {
      for(uint32_t i = 0; i < c.replen; ++i) {
	c.fwdmotif.insert(boost::to_upper_copy(repset[k].substr(i, c.replen * c.period)));
      }
    }
    
    // Create reverse motifs
    for(std::set<std::string>::iterator itr = c.fwdmotif.begin(); itr != c.fwdmotif.end(); ++itr) {
      std::string s(*itr);
      reverseComplement(s);
      c.revmotif.insert(s);
      //std::cerr << *itr << ',' << s << std::endl;
    }
  }
    
  inline double
  entropy(std::string const& st) {
    typedef double TPrecision;
    std::vector<char> stvec(st.begin(), st.end());
    std::set<char> alphabet(stvec.begin(), stvec.end());
    TPrecision ent = 0;
    for(std::set<char>::const_iterator c = alphabet.begin(); c != alphabet.end(); ++c) {
      int ctr = 0;
      for (std::vector<char>::const_iterator s = stvec.begin(); s != stvec.end(); ++s)
	if (*s == *c) ++ctr;
      TPrecision freq = (TPrecision) ctr / (TPrecision) stvec.size();
      ent += (freq) * log(freq)/log(2);
    }
    return -ent;
  }

  
  inline void
  estimateCoverage(std::string const& filename, std::string const& genomefile, uint32_t const win, uint64_t& avgcov, uint64_t& sdcov) {
    samFile* samfile = sam_open(filename.c_str(), "r");
    hts_set_fai_filename(samfile, genomefile.c_str());
    hts_idx_t* idx = sam_index_load(samfile, filename.c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Sort chromosomes by length
    typedef std::pair<uint32_t, uint32_t> TLengthChr;
    std::vector<TLengthChr> chrByLength;
    for(int32_t i = 0; i < hdr->n_targets; ++i) chrByLength.push_back(std::make_pair(hdr->target_len[i], i));
    std::sort(chrByLength.begin(), chrByLength.end(), std::greater<TLengthChr>());

    // Collect coverage estimates
    std::vector<uint64_t> vcov;
    
    // Sample windows from the 10 longest chromosomes
    for(uint32_t k = 0; ((k < 10) && (k < chrByLength.size())); ++k) {
      for(uint32_t i = 10; i <= 90; ++i) {
	if ((i >= 20) && (i <=80)) continue;
	uint32_t regstart = ((float) i / 100.0) * chrByLength[k].first;
	uint32_t regend = regstart + win;
	if (regend > hdr->target_len[chrByLength[k].second]) continue;
	std::vector<uint32_t> cov(win, 0);
	// Read alignments
	hts_itr_t* iter = sam_itr_queryi(idx, chrByLength[k].second, regstart, regend);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile, iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FUNMAP)) continue;
	  uint32_t rp = rec->core.pos; // reference pointer
	  uint32_t sp = 0; // sequence pointer
	  uint32_t* cigar = bam_get_cigar(rec);
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
		if ((rp >= regstart) && (rp < regend)) ++cov[rp-regstart];
		++rp;
		++sp;
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      rp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	      // Nop
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      rp += bam_cigar_oplen(cigar[i]);
	    } else {
	      std::cerr << "Warning: Unknown Cigar operation!" << std::endl;
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
	uint64_t covsum = 0;
	for(uint32_t j = 0; j < cov.size(); ++j) covsum += cov[j];
	vcov.push_back(covsum);
      }
    }

    // Drop lowest and highest 25%
    uint32_t ist = 0;
    uint32_t ien = vcov.size();
    if (ien > 100) {
      ist = (uint32_t) (0.25 * vcov.size());
      ien = (uint32_t) (0.75 * vcov.size());
      std::sort(vcov.begin(), vcov.end());
    }
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance> > acc;
    for(uint32_t i = ist; i < ien; ++i) acc(vcov[i]);
    sdcov = sqrt(boost::accumulators::variance(acc));
    avgcov = boost::accumulators::mean(acc);
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }

  inline uint32_t sequenceLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t slen = 0;
    for (uint32_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CINS) || (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) slen += bam_cigar_oplen(cigar[i]);
    return slen;
  }  

  inline std::size_t
  hash_string(const char *s) {
    std::size_t h = 37;
    while (*s) {
      h = (h * 54059) ^ (s[0] * 76963);
      s++;
    }
    return h;
  }

  inline std::size_t hash_lr(bam1_t* rec) {
    boost::hash<std::string> string_hash;
    std::string qname = bam_get_qname(rec);
    std::size_t seed = hash_string(qname.c_str());
    boost::hash_combine(seed, string_hash(qname));
    return seed;
  }

  inline uint32_t
  alignmentLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (uint32_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) alen += bam_cigar_oplen(cigar[i]);
    return alen;
  }
  
  inline uint32_t
  lastAlignedPosition(bam1_t const* rec) {
    return rec->core.pos + alignmentLength(rec);
  }

}

#endif
