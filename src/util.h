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
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>

namespace lorax
{

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
