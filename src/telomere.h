#ifndef TELOMERE_H
#define TELOMERE_H

#include <fstream>
#include <iomanip>

#include <boost/dynamic_bitset.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/icl/split_interval_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

namespace lorax
{

  struct TelomereConfig {
    uint16_t minSeqQual;
    uint32_t maxOffset;
    uint32_t minChrEndDist;
    uint32_t minChrLen;
    uint32_t minTelmoSize;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
    boost::filesystem::path tumor;
    boost::filesystem::path control;
  };


  struct Mapping {
    int32_t cid;
    int32_t tid;
    int32_t gstart;
    int32_t gend;
    int32_t rstart;
    int32_t rend;
    int32_t telmo;
    bool fwd;
    bool chrend;
    uint16_t qual;
    std::size_t seed;
    std::string qname;

    Mapping(int32_t const cd, int32_t const t, int32_t const gs, int32_t const ge, int32_t const rs, int32_t const re, int32_t const telm, bool const val, bool const cre, uint16_t const qval, std::size_t s, std::string const& qn) : cid(cd), tid(t), gstart(gs), gend(ge), rstart(rs), rend(re), telmo(telm), fwd(val), chrend(cre), qual(qval), seed(s), qname(qn) {}
  };


  template<typename TConfig, typename TEdges>
  inline void
  concomp(TConfig const&, std::vector<Mapping>& mp, TEdges& es) {
    // Segments are nodes, split-reads are edges
    for(uint32_t id1 = 0; id1 < mp.size(); ++id1) {
      for(uint32_t id2 = id1 + 1; id2 < mp.size(); ++id2) {
	if (es.find(std::make_pair(id1, id2)) != es.end()) {
	  // Vertices in different components?
	  if (mp[id1].cid != mp[id2].cid) {
	    int32_t oldid = mp[id2].cid;
	    for(uint32_t i = 0; i < mp.size(); ++i) {
	      if (mp[i].cid == oldid) mp[i].cid = mp[id1].cid;
	    }
	  }
	}
      }
    }
  }
  
  
  template<typename TConfig, typename TEdges>
  inline void
  links(TConfig const& c, std::vector<Mapping> const& mp, TEdges& es) {
    // Valid reads
    std::set<std::size_t> valid_read;
    // Ignore reads that lack a non-telomere mapping
    for(uint32_t i = 0; i < mp.size(); ++i) {
      if (!mp[i].chrend) valid_read.insert(mp[i].seed);
    }    
    
    // Segments are nodes, read connections define edges
    for(uint32_t id1 = 0; id1 < mp.size(); ++id1) {
      if (valid_read.find(mp[id1].seed) != valid_read.end()) {
	for(uint32_t id2 = id1 + 1; id2 < mp.size(); ++id2) {
	  if (valid_read.find(mp[id2].seed) != valid_read.end()) {
	    if ((mp[id1].seed == mp[id2].seed)) {
	      // Same read, check intra-read offset
	      if ((std::abs(mp[id1].rend - mp[id2].rstart) < c.maxOffset) || (std::abs(mp[id1].rstart - mp[id2].rend) < c.maxOffset)) {
		es.insert(std::make_pair(id1, id2));
	      }
	    } else {
	      // Different read, only link via non-telomere mappings
	      if (mp[id1].tid == mp[id2].tid) {
		if ((!mp[id1].chrend) && (!mp[id2].chrend)) {
		  if (!((mp[id1].gend < mp[id2].gstart) or (mp[id1].gstart > mp[id2].gend))) {
		    es.insert(std::make_pair(id1, id2));
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  

  template<typename TConfig>
  inline int32_t
  collectCandidates(TConfig const& c, std::string const& filename, std::vector<std::string> const& motifs, std::set<std::size_t>& candidates) {
    bool tumor_run = false;
    if (filename == c.tumor.string()) tumor_run = true;
    int32_t seqMotifSize = motifs[0].size();

    // Open file handles
    samFile* samfile = sam_open(filename.c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, filename.c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Candidate reads
    std::set<std::size_t> telomere_reads;
    std::set<std::size_t> nontel_reads;
    
    // Parse BAM
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[refIndex] << std::endl;

      // Iterate alignments 
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FUNMAP)) continue;
	if (rec->core.qual <= c.minSeqQual) continue;
	std::size_t seed = hash_string(bam_get_qname(rec));
      
	// Load sequence and quality
	typedef std::vector<uint8_t> TQuality;
	TQuality quality(rec->core.l_qseq);
	std::string sequence(rec->core.l_qseq, 'N');
	uint8_t* seqptr = bam_get_seq(rec);
	uint8_t* qualptr = bam_get_qual(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) {
	  quality[i] = qualptr[i];
	  sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	}
	
	// Telomere motif length
	uint32_t telmo = 0;
	for(uint32_t i = 0; i < motifs.size(); ++i) {
	  uint32_t telmo_i = 0;
	  std::size_t pos = sequence.find(motifs[i], 0);
	  while (pos != std::string::npos) {
	    // Sufficient quality?
	    int32_t avgqual = 0;
	    for(uint32_t k = pos; ((k < pos + seqMotifSize) && (k < quality.size())); ++k) avgqual += (int32_t) quality[k];
	    avgqual /= seqMotifSize;
	    if (avgqual > c.minSeqQual) telmo_i += seqMotifSize;
	    pos = sequence.find(motifs[i], pos+seqMotifSize);
	  }
	  if (telmo_i > telmo) {
	    telmo = telmo_i;
	    if ((telmo >= c.minTelmoSize) || ((!tumor_run) && ((int) telmo >= seqMotifSize))) {
	      telomere_reads.insert(seed);
	      break;
	    }
	  }
	}
	
	// Check for telomere distal mapping
	if ((hdr->target_len[rec->core.tid] > c.minChrLen)  && ((hdr->target_len[rec->core.tid] - rec->core.pos) > c.minChrEndDist) && (rec->core.pos > c.minChrEndDist)) nontel_reads.insert(seed);
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    // Take intersection
    std::set_intersection(telomere_reads.begin(), telomere_reads.end(), nontel_reads.begin(), nontel_reads.end(), std::inserter(candidates, candidates.begin()));

    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    return candidates.size();
  }


  template<typename TConfig>
  inline void
  segmentOut(TConfig const& c, std::vector<Mapping> const& mp) {
    // Open file handles
    samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.tumor.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Pre-compute support
    typedef std::pair<std::size_t, uint32_t> TReadComp;
    typedef std::set<TReadComp> TRCSet;
    TRCSet rcs;
    for(uint32_t i = 0; i < mp.size(); ++i) rcs.insert(std::make_pair(mp[i].seed, mp[i].cid));
    std::vector<uint32_t> sup(mp.size(), 0);
    for(TRCSet::iterator it = rcs.begin(); it != rcs.end(); ++it) ++sup[it->second];
    
    // Output file
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    dataOut << "chr\trefstart\trefend\treadname\treadstart\treadend\tcomponentid\tsupport\tforward\ttelmotiflen\tchrend" << std::endl;

    for(uint32_t i = 0; i < mp.size(); ++i) {
      // Ignore singletons
      if (sup[mp[i].cid] > 1) {
	dataOut << hdr->target_name[mp[i].tid] << '\t' << mp[i].gstart << '\t' << mp[i].gend << '\t' << mp[i].qname << '\t' << mp[i].rstart << '\t' << mp[i].rend << '\t' << mp[i].cid << '\t' << sup[mp[i].cid] << '\t' << (int) (mp[i].fwd) << '\t' << mp[i].telmo << '\t' << (int) (mp[i].chrend) << std::endl;
      }
    }

    // Close output file
    dataOut.pop();
    dataOut.pop();
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }

  template<typename TConfig>
  inline void
  mappings(TConfig const& c, std::string const& filename, std::vector<std::string> const& motifs, std::set<std::size_t> const& candidates, std::vector<Mapping>& mp, std::vector<Mapping> const& ctrl) {
    bool tumor_run = false;
    if (filename == c.tumor.string()) tumor_run = true;
    int32_t seqMotifSize = motifs[0].size();

    // Open file handles
    samFile* samfile = sam_open(filename.c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, filename.c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Parse BAM
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[refIndex] << std::endl;

      // Mask control regions
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet mask;
      if (tumor_run) {
	mask.resize(hdr->target_len[refIndex], 0);
	for(uint32_t i = 0; i < ctrl.size(); ++i) {
	  if (ctrl[i].tid == refIndex) {
	    for(int32_t k = ctrl[i].gstart; ((k < ctrl[i].gend) && (k < (int) hdr->target_len[refIndex])); ++k) mask[k] = 1;
	  }
	}
      }
      
      // Iterate alignments 
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY)) continue;
	if (rec->core.qual <= c.minSeqQual) continue;
	std::size_t seed = hash_string(bam_get_qname(rec));
	if (candidates.find(seed) != candidates.end()) {
      
	  // Get read sequence
	  std::string sequence;
	  sequence.resize(rec->core.l_qseq);
	  uint8_t* seqptr = bam_get_seq(rec);
	  for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	
	  // Parse CIGAR
	  uint32_t* cigar = bam_get_cigar(rec);
	  int32_t gp = rec->core.pos; // Genomic position
	  int32_t gpStart = -1; //Match start
	  int32_t gpEnd = -1; //Match end
	  int32_t sp = 0; // Sequence position
	  int32_t seqStart = -1;  // Match start
	  int32_t seqEnd = -1; // Match end
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      if (seqStart == -1) {
		seqStart = sp;
		gpStart = gp;
	      }
	      gp += bam_cigar_oplen(cigar[i]);
	      sp += bam_cigar_oplen(cigar[i]);
	      seqEnd = sp;
	      gpEnd = gp;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      if (seqStart == -1) {
		seqStart = sp;
		gpStart = gp;
	      }
	      sp += bam_cigar_oplen(cigar[i]);
	      seqEnd = sp;
	      gpEnd = gp;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      if (seqStart == -1) {
		seqStart = sp;
		gpStart = gp;
	      }
	      gp += bam_cigar_oplen(cigar[i]);
	      seqEnd = sp;
	      gpEnd = gp;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      gp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else {
	      std::cerr << "Warning: Unknown Cigar options!" << std::endl;
	    }
	  }
	  
	  // Telomere motif length
	  uint32_t telmo = 0;
	  for(uint32_t i = 0; i < motifs.size(); ++i) {
	    uint32_t telmo_i = 0;
	    std::size_t pos = sequence.find(motifs[i], 0);
	    while (pos != std::string::npos) {
	      telmo_i += seqMotifSize;
	      pos = sequence.find(motifs[i], pos+seqMotifSize);
	    }
	    if (telmo_i > telmo) telmo = telmo_i;
	  }
	  bool dir = true;
	  if (rec->core.flag & BAM_FREVERSE) {
	    dir = false;
	    int32_t seqTmp = seqStart;
	    seqStart = sp - seqEnd;
	    seqEnd = sp - seqTmp;
	  }
	  bool novelSegment = true;
	  int32_t wiggle = 50;
	  for(uint32_t i = 0; i < mp.size(); ++i) {
	    if (mp[i].seed == seed) {
	      if ((seqStart + wiggle >= mp[i].rstart) && (seqEnd <= mp[i].rend + wiggle)) {
		novelSegment = false;
		break;
	      } else if ((mp[i].rstart + wiggle >= seqStart) && (mp[i].rend <= seqEnd + wiggle)) {
		// Replace entry
		mp[i].tid = rec->core.tid;
		mp[i].gstart = gpStart;
		mp[i].gend = gpEnd;
		mp[i].rstart = seqStart;
		mp[i].rend = seqEnd;
		mp[i].telmo = telmo;
		mp[i].fwd = dir;
		novelSegment = false;
		break;
	      }
	    }
	  }
	  if ((novelSegment) && (gpStart < gpEnd)) {
	    bool chrend = true;
	    if (((hdr->target_len[rec->core.tid] - gpEnd) > c.minChrEndDist) && (gpStart > (int32_t) c.minChrEndDist)) chrend = false;
	    bool incl_mapping = false;
	    if (tumor_run) {
	      // For the tumor, check that non-telomere mappings are absent in the control
	      incl_mapping = true;
	      if (!chrend) {
		for(int32_t k = gpStart; ((k < gpEnd) && (k < (int) hdr->target_len[refIndex])); ++k) {
		  if (mask[k]) {
		    incl_mapping = false;
		    break;
		  }
		}
	      }
	    } else {
	      // For the control, only include non-telomere mappings
	      incl_mapping = !chrend;
	    }
	    if (incl_mapping) {
	      //std::cerr << bam_get_qname(rec) << ',' << hdr->target_name[rec->core.tid] << ',' << gpStart << std::endl;
	      mp.push_back(Mapping(mp.size(), rec->core.tid, gpStart, gpEnd, seqStart, seqEnd, telmo, dir, chrend, rec->core.qual, seed, bam_get_qname(rec)));
	    }
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }
    

  template<typename TConfig>
  inline int32_t
  runTelomere(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    std::vector<std::string> motifs = { "CCCTCACCCTAACCCTCA", "CCCTGACCCTGACCCCAA", "CCCCAACCCTAACCCTCA", "CCCCAACCCTAACCCTGA", "TCAGGGTTAGGGTTGGGG", "TGAGGGTGAGGGTCAGGG", "CCCTAACCCTGACCCTAA", "TTGGGGTTGGGGTTGGGG", "CCCTAACCCTCACCCTGA", "TTGGGGTCAGGGTTGGGG", "CCCTCACCCCAACCCCAA", "TTAGGGTCAGGGTTAGGG", "CCCTGACCCTAACCCTAA", "TGAGGGTCAGGGTTAGGG", "CCCCAACCCTCACCCTAA", "CCCTCACCCCAACCCTCA", "CCCCAACCCTCACCCTGA", "TCAGGGTCAGGGTTGGGG", "TGAGGGTCAGGGTTGGGG", "TTGGGGTTAGGGTCAGGG", "TTGGGGTTAGGGTTAGGG", "CCCTCACCCTCACCCTCA", "CCCTGACCCTAACCCCAA", "CCCCAACCCTGACCCTCA", "TTAGGGTTGGGGTCAGGG", "CCCTGACCCTGACCCTAA", "CCCCAACCCTAACCCTAA", "TGAGGGTTGGGGTTGGGG", "CCCTAACCCTAACCCTGA", "CCCTGACCCTCACCCTAA", "TTAGGGTTGGGGTTGGGG", "TCAGGGTTGGGGTGAGGG", "CCCTCACCCTAACCCCAA", "CCCCAACCCCAACCCTCA", "TCAGGGTTGGGGTTAGGG", "TCAGGGTTAGGGTCAGGG", "TTAGGGTGAGGGTTGGGG", "TTGGGGTTAGGGTTGGGG", "CCCCAACCCTCACCCCAA", "TTAGGGTGAGGGTGAGGG", "CCCTGACCCCAACCCTCA", "CCCTCACCCTAACCCTAA", "CCCTGACCCCAACCCCAA", "TCAGGGTCAGGGTTAGGG", "TGAGGGTTGGGGTGAGGG", "CCCTCACCCTCACCCCAA", "CCCTAACCCTGACCCTGA", "CCCCAACCCTGACCCTGA", "CCCTAACCCTAACCCTAA", "CCCTGACCCTAACCCTGA", "CCCTCACCCCAACCCTAA", "CCCCAACCCTGACCCCAA", "TGAGGGTTGGGGTCAGGG", "CCCTAACCCCAACCCTGA", "TCAGGGTGAGGGTCAGGG", "TTGGGGTTAGGGTGAGGG", "TTAGGGTGAGGGTCAGGG", "TGAGGGTTAGGGTTGGGG", "TTGGGGTGAGGGTTGGGG", "TGAGGGTTGGGGTTAGGG", "TTAGGGTTGGGGTGAGGG", "CCCCAACCCTGACCCTAA", "TTAGGGTCAGGGTGAGGG", "CCCTGACCCTCACCCTCA", "TGAGGGTTAGGGTTAGGG", "TTAGGGTTGGGGTTAGGG", "TCAGGGTCAGGGTGAGGG", "CCCCAACCCTCACCCTCA", "TTGGGGTGAGGGTTAGGG", "TCAGGGTTGGGGTCAGGG", "TTGGGGTTGGGGTTAGGG", "TTAGGGTTAGGGTTGGGG", "CCCTAACCCTCACCCTCA", "CCCTCACCCCAACCCTGA", "TTAGGGTTAGGGTGAGGG", "TTGGGGTCAGGGTGAGGG", "TTGGGGTTGGGGTCAGGG", "CCCTCACCCTAACCCTGA", "TTGGGGTCAGGGTCAGGG", "TCAGGGTGAGGGTGAGGG", "CCCTGACCCTGACCCTGA", "TTGGGGTCAGGGTTAGGG", "CCCTGACCCCAACCCTAA", "TTAGGGTTAGGGTTAGGG", "TGAGGGTCAGGGTCAGGG", "CCCTCACCCTGACCCTGA", "TTAGGGTCAGGGTTGGGG", "CCCTAACCCTAACCCTCA", "CCCTAACCCCAACCCTAA", "CCCTAACCCTGACCCCAA", "TCAGGGTCAGGGTCAGGG", "CCCTCACCCTCACCCTGA", "TTAGGGTGAGGGTTAGGG", "CCCCAACCCCAACCCTGA", "CCCTAACCCTCACCCCAA", "CCCTGACCCTCACCCTGA", "TTGGGGTGAGGGTGAGGG", "TTGGGGTGAGGGTCAGGG", "TCAGGGTTGGGGTTGGGG", "TGAGGGTTAGGGTGAGGG", "TCAGGGTTAGGGTTAGGG", "TGAGGGTCAGGGTGAGGG", "CCCCAACCCTAACCCCAA", "TGAGGGTGAGGGTTAGGG", "CCCTAACCCTCACCCTAA", "TGAGGGTGAGGGTGAGGG", "TTAGGGTCAGGGTCAGGG", "CCCTAACCCTAACCCCAA", "CCCTAACCCTGACCCTCA", "TCAGGGTTAGGGTGAGGG", "CCCTAACCCCAACCCCAA", "CCCTGACCCTAACCCTCA", "TCAGGGTGAGGGTTGGGG", "TTAGGGTTAGGGTCAGGG", "CCCTGACCCCAACCCTGA", "CCCCAACCCCAACCCCAA", "CCCTCACCCTCACCCTAA", "CCCTCACCCTGACCCCAA", "TGAGGGTGAGGGTTGGGG", "CCCTCACCCTGACCCTAA", "TCAGGGTGAGGGTTAGGG", "CCCTAACCCCAACCCTCA", "TGAGGGTTAGGGTCAGGG", "CCCTGACCCTGACCCTCA", "CCCCAACCCCAACCCTAA", "TTGGGGTTGGGGTGAGGG", "CCCTGACCCTCACCCCAA", "CCCTCACCCTGACCCTCA" };

    std::vector<Mapping> tumor_mp;
    if (tumor_mp.empty()) {
      // Collect control regions
      std::vector<Mapping> contr_mp;
    
      // Parse matched control
      std::set<std::size_t> control_reads;
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Collecting control BAM telomere reads." << std::endl;
      int32_t contr_size = collectCandidates(c, c.control.string(), motifs, control_reads);

      // Get alignments
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Extracting alignments for " << contr_size << " reads." << std::endl;
      mappings(c, c.control.string(), motifs, control_reads, contr_mp, tumor_mp);
    
      // Candidate tumor reads
      std::set<std::size_t> tumor_reads;
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Collecting tumor BAM telomere reads." << std::endl;
      int32_t tum_size = collectCandidates(c, c.tumor.string(), motifs, tumor_reads);
    
      // Get alignments
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Extracting alignments for " << tum_size << " reads." << std::endl;
      mappings(c, c.tumor.string(), motifs, tumor_reads, tumor_mp, contr_mp);
    }

    // Connect mappings
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing " << tumor_mp.size() << " mappings." << std::endl;
    typedef std::pair<uint32_t, uint32_t> TEdge;
    std::set<TEdge> es;
    links(c, tumor_mp, es);

    // Compute connected components
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Connected components using " << es.size() << " edges." << std::endl;
    concomp(c, tumor_mp, es);

    // Data out
    segmentOut(c, tumor_mp);
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    return 0;
  }


  
  int telomere(int argc, char** argv) {
    TelomereConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Options");
    generic.add_options()
      ("help,?", "show help message")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minSeqQual)->default_value(10), "min. mapping quality")
      ("offset,f", boost::program_options::value<uint32_t>(&c.maxOffset)->default_value(5000), "max. basepair offset to connect read segments")
      ("chrend,c", boost::program_options::value<uint32_t>(&c.minChrEndDist)->default_value(1000000), "min. distance to chromosome end for fusion part")
      ("chrlen,l", boost::program_options::value<uint32_t>(&c.minChrLen)->default_value(40000000), "min. chromosome length")
      ("telsize,t", boost::program_options::value<uint32_t>(&c.minTelmoSize)->default_value(36), "min. telomere motif length")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("matched,m", boost::program_options::value<boost::filesystem::path>(&c.control), "matched control BAM")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.bed.gz"), "gzipped output file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.tumor), "input file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome")) || (!vm.count("matched"))) {
      std::cout << "Usage: lorax " << argv[0] << " [OPTIONS] -g <ref.fa> -m <control.bam> <tumor.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runTelomere(c);
  }

}

#endif
