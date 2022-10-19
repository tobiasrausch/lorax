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
    int32_t refidx;
    int32_t refpos;
    int32_t seqpos;
    uint16_t qual;
    
    Junction(bool const fw, bool const cl, int32_t const idx, int32_t const r, int32_t const s, uint16_t const qval) : forward(fw), scleft(cl), refidx(idx), refpos(r), seqpos(s), qual(qval) {}
  };


  template<typename TReadBp>
  inline void
  _insertJunction(TReadBp& readBp, std::size_t const seed, bam1_t* rec, int32_t const rp, int32_t const sp, bool const scleft) {
    bool fw = true;
    if (rec->core.flag & BAM_FREVERSE) fw = false;
    typedef typename TReadBp::mapped_type TJunctionVector;
    typename TReadBp::iterator it = readBp.find(seed);
    int32_t seqlen = sequenceLength(rec);
    if (sp <= seqlen) {
      if (rec->core.flag & BAM_FREVERSE) {
	if (it != readBp.end()) it->second.push_back(Junction(fw, scleft, rec->core.tid, rp, seqlen - sp, rec->core.qual));
	else readBp.insert(std::make_pair(seed, TJunctionVector(1, Junction(fw, scleft, rec->core.tid, rp, seqlen - sp, rec->core.qual))));
      } else {
	if (it != readBp.end()) it->second.push_back(Junction(fw, scleft, rec->core.tid, rp, sp, rec->core.qual));
	else readBp.insert(std::make_pair(seed, TJunctionVector(1, Junction(fw, scleft, rec->core.tid, rp, sp, rec->core.qual))));
      }
    }
  }

  template<typename TJunction>
  struct SortJunction : public std::binary_function<TJunction, TJunction, bool>
  {
    inline bool operator()(TJunction const& j1, TJunction const& j2) {
      return ((j1.seqpos<j2.seqpos) || ((j1.seqpos==j2.seqpos) && (j1.refidx<j2.refidx)) || ((j1.seqpos==j2.seqpos) && (j1.refidx==j2.refidx) && (j1.refpos<j2.refpos)));
    }
  };


  template<typename TConfig, typename TReadBp>
  inline void
  findJunctions(TConfig const& c, TReadBp& readBp) {

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
	if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;

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
	    sp += bam_cigar_oplen(cigar[i]);
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
	    if (bam_cigar_oplen(cigar[i]) > c.minClip) {
	      // Is this a non-telomere junction?
	      if ((rp > c.maxTelLen) && (rp + c.maxTelLen < hdr->target_len[refIndex])) _insertJunction(readBp, seed, rec, rp, finalsp, scleft);
	    }
	  } else {
	    std::cerr << "Unknown Cigar options" << std::endl;
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }

    // Sort junctions
    for(typename TReadBp::iterator it = readBp.begin(); it != readBp.end(); ++it) std::sort(it->second.begin(), it->second.end(), SortJunction<Junction>());


    // Debug
    //for(typename TReadBp::iterator it = readBp.begin(); it != readBp.end(); ++it) {
    //for(uint32_t i = 0; i < it->second.size(); ++i) {
    //std::cerr << it->first << '\t' << hdr->target_name[it->second[i].refidx] << ':' << it->second[i].refpos << '\t' << it->second[i].seqpos << '\t' << it->second[i].forward << '\t' << '\t' << it->second[i].qual << std::endl;
    //}
    //}
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }

}

#endif
