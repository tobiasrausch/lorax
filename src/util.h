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
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>

namespace lorax
{

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
