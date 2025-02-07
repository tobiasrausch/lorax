#ifndef JUNCTION_H
#define JUNCTION_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace lorax
{

  struct Junction {
    bool forward;
    bool scleft;
    bool telLeft;
    bool telRight;
    int32_t refidx;
    int32_t refpos;
    int32_t seqpos;
    uint16_t qual;
    std::size_t seed;

    Junction(std::size_t const s) : refidx(0), refpos(0), seed(s) {}
    
    Junction(bool const fw, bool const cl, int32_t const idx, int32_t const r, int32_t const s, uint16_t const qval, std::size_t seedval) : forward(fw), scleft(cl), telLeft(false), telRight(false), refidx(idx), refpos(r), seqpos(s), qual(qval), seed(seedval) {}

    bool operator<(const Junction& j2) {
      return ((seed<j2.seed) || ((seed == j2.seed) && (refidx<j2.refidx)) || ((seed == j2.seed) && (refidx==j2.refidx) && (refpos<j2.refpos)));
    }
  };

  template<typename TReadBp>
  inline void
  _insertJunction(TReadBp& readBp, std::size_t const seed, bam1_t* rec, int32_t const rp, int32_t const sp, bool const scleft) {
    bool fw = true;
    if (rec->core.flag & BAM_FREVERSE) fw = false;
    int32_t seqlen = sequenceLength(rec);
    if (sp <= seqlen) {
      if (rec->core.flag & BAM_FREVERSE) readBp.push_back(Junction(fw, scleft, rec->core.tid, rp, seqlen - sp, rec->core.qual, seed));
      else readBp.push_back(Junction(fw, scleft, rec->core.tid, rp, sp, rec->core.qual, seed));
    }
  }


  template<typename TConfig>
  inline void
  findJunctions(TConfig const& c, std::vector<Junction>& readBp) {
    // Open file handles
    samFile* samfile = sam_open(c.sample.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.sample.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Parse genome chr-by-chr
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Split-read scanning" << std::endl;

    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      if (hdr->target_len[refIndex] < c.minChrLen) continue;
      
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Processing chromosome " << hdr->target_name[refIndex] << std::endl;
      
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {

	// Keep secondary alignments
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	std::size_t seed = hash_lr(rec);
	uint32_t rp = rec->core.pos; // reference pointer
	uint32_t sp = 0; // sequence pointer

	// Parse the CIGAR
	uint32_t* cigar = bam_get_cigar(rec);
	for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	  if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	    sp += bam_cigar_oplen(cigar[i]);
	    rp += bam_cigar_oplen(cigar[i]);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	    rp += bam_cigar_oplen(cigar[i]);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	    if (bam_cigar_oplen(cigar[i]) > c.minClip) _insertJunction(readBp, seed, rec, rp, sp, false);
	    sp += bam_cigar_oplen(cigar[i]);
	    if (bam_cigar_oplen(cigar[i]) > c.minClip) _insertJunction(readBp, seed, rec, rp, sp, true);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	    rp += bam_cigar_oplen(cigar[i]);
	  } else if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
	    int32_t finalsp = sp;
	    bool scleft = false;
	    if (sp == 0) {
	      finalsp += bam_cigar_oplen(cigar[i]); // Leading soft-clip / hard-clip
	      scleft = true;
	    }
	    sp += bam_cigar_oplen(cigar[i]);
	    if (bam_cigar_oplen(cigar[i]) > c.minClip) _insertJunction(readBp, seed, rec, rp, finalsp, scleft);
	  } else {
	    std::cerr << "Unknown Cigar options" << std::endl;
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }

    // Sort junctions
    std::sort(readBp.begin(), readBp.end());

    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }


}

#endif
