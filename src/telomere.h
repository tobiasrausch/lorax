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
    uint32_t segdev;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
    boost::filesystem::path tumor;
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
    std::size_t seed;
    std::string qname;

    Mapping(int32_t const cd, int32_t const t, int32_t const gs, int32_t const ge, int32_t const rs, int32_t const re, int32_t const telm, bool const val, std::size_t s, std::string const& qn) : cid(cd), tid(t), gstart(gs), gend(ge), rstart(rs), rend(re), telmo(telm), fwd(val), seed(s), qname(qn) {}
  };


  template<typename TConfig, typename TEdges>
  inline void
  concomp(TConfig const& c, std::vector<Mapping>& mp, TEdges& es) {
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
    // Segments are nodes, read connections define edges
    for(uint32_t id1 = 0; id1 < mp.size(); ++id1) {
      for(uint32_t id2 = id1 + 1; id2 < mp.size(); ++id2) {
	// Same read?
	if ((mp[id1].seed == mp[id2].seed)) {
	  // Check intra-read offset
	  if ((std::abs(mp[id1].rend - mp[id2].rstart) < c.maxOffset) || (std::abs(mp[id1].rstart - mp[id2].rend) < c.maxOffset)) {
	    es.insert(std::make_pair(id1, id2));
	  }
	} else {
	  // Any segment overlap
	  if (mp[id1].tid == mp[id2].tid) {
	    if (!((mp[id1].gend < mp[id2].gstart) or (mp[id1].gstart > mp[id2].gend))) {
	      // Check segment boundaries
	      if ((std::abs(mp[id1].gstart - mp[id2].gstart) < c.segdev) || (std::abs(mp[id1].gend - mp[id2].gend) < c.segdev)) {
		es.insert(std::make_pair(id1, id2));
	      }
	    }
	  }
	}
      }
    }
  }
  

  template<typename TConfig>
  inline int32_t
  collectCandidates(TConfig const& c, std::vector<std::string> const& motifs, std::set<std::size_t>& candidates) {
    int32_t seqMotifSize = motifs[0].size();

    // Open file handles
    samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.tumor.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Candidate reads
    std::set<std::size_t> telomere_reads;
    std::set<std::size_t> inter_chr_reads;
    std::map<std::size_t, int32_t> reads_chr;

    // Parse BAM
    int32_t oldId = -1;
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      if (rec->core.tid != oldId) {
	oldId = rec->core.tid;
	if ((oldId >= 0) && (oldId < hdr->n_targets)) {
	  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[oldId] << std::endl;
	}
      }
      if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP)) continue;
      std::size_t seed = hash_string(bam_get_qname(rec));
      
      // Check for inter-chromosomal mapping
      if (!(rec->core.flag & BAM_FUNMAP)) {
	if (reads_chr.find(seed) != reads_chr.end()) {
	  if (reads_chr[seed] != rec->core.tid) inter_chr_reads.insert(seed);
	  } else {
	  reads_chr.insert(std::make_pair(seed, rec->core.tid));
	}
      }
      
      // Load sequence and quality
      if (rec->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
      typedef std::vector<uint8_t> TQuality;
      TQuality quality(rec->core.l_qseq);
      std::string sequence(rec->core.l_qseq, 'N');
      uint8_t* seqptr = bam_get_seq(rec);
      uint8_t* qualptr = bam_get_qual(rec);
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) {
	quality[i] = qualptr[i];
	sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
      }
      
      // Search motifs, count only first hit
      for(uint32_t i = 0; i < motifs.size(); ++i) {
	std::size_t pos = sequence.find(motifs[i]);
	if (pos != std::string::npos) {
	  // Sufficient quality?
	  int32_t avgqual = 0;
	  for(uint32_t k = pos; ((k < pos + seqMotifSize) && (k < quality.size())); ++k) avgqual += (int32_t) quality[k];
	  avgqual /= seqMotifSize;
	  if (avgqual > c.minSeqQual) {
	    telomere_reads.insert(seed);
	    break; // Max. one hit per sequence
	  }
	}
      }
    }
    bam_destroy1(rec);

    // Take intersection
    std::set_intersection(telomere_reads.begin(), telomere_reads.end(), inter_chr_reads.begin(), inter_chr_reads.end(), std::inserter(candidates, candidates.begin()));

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
    dataOut << "component\tsupport\tchr\trefstart\trefend\treadname\treadstart\treadend\tforward\ttelmotiflen" << std::endl;

    for(uint32_t i = 0; i < mp.size(); ++i) {
      dataOut << mp[i].cid << '\t' << sup[mp[i].cid] << '\t' << hdr->target_name[mp[i].tid] << '\t' << mp[i].gstart << '\t' << mp[i].gend << '\t' << mp[i].qname << '\t' << mp[i].rstart << '\t' << mp[i].rend << '\t' << (int) (mp[i].fwd) << '\t' << mp[i].telmo << std::endl;
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
  mappings(TConfig const& c, std::vector<std::string> const& motifs, std::set<std::size_t> const& candidates, std::vector<Mapping>& mp) {
    int32_t seqMotifSize = motifs[0].size();

    // Open file handles
    samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.tumor.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Parse BAM
    int32_t oldId = -1;
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      if (rec->core.tid != oldId) {
	oldId = rec->core.tid;
	if ((oldId >= 0) && (oldId < hdr->n_targets)) {
	  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[oldId] << std::endl;
	}
      }
      if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
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
	uint32_t subSeqStart = 0; // Aligned sequence start
	uint32_t subSeqEnd = 0; // Aligned sequence end
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
	    subSeqEnd += bam_cigar_oplen(cigar[i]);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	    if (seqStart == -1) {
	      seqStart = sp;
	      gpStart = gp;
	    }
	    sp += bam_cigar_oplen(cigar[i]);
	    seqEnd = sp;
	    gpEnd = gp;
	    subSeqEnd += bam_cigar_oplen(cigar[i]);
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
	    if (subSeqEnd == 0) {
	      // Leading soft-clip
	      subSeqStart += bam_cigar_oplen(cigar[i]);
	      subSeqEnd = subSeqStart;
	    }
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
	  std::size_t pos = sequence.find(motifs[i], subSeqStart);
	  while ((pos != std::string::npos) && (pos < subSeqEnd)) {
	    telmo_i += seqMotifSize;
	    pos = sequence.find(motifs[i], pos+seqMotifSize);
	  }
	  if (telmo_i > telmo) telmo = telmo_i;
	}
	
	bool dir = true;
	if (rec->core.flag & BAM_FREVERSE) dir = false;
	mp.push_back(Mapping(mp.size(), rec->core.tid, gpStart, gpEnd, seqStart, seqEnd, telmo, dir, seed, bam_get_qname(rec)));
      }
    }
    // Clean-up
    bam_destroy1(rec);
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
    
    // Candidate reads
    std::set<std::size_t> candidates;
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Collecting candidate reads." << std::endl;
    int32_t candsize = collectCandidates(c, motifs, candidates);
    
    // Get alignments
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Extracting alignments for " << candsize << " reads." << std::endl;
    std::vector<Mapping> mp;
    mappings(c, motifs, candidates, mp);

    // Connect mappings
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Process " << mp.size() << " mappings." << std::endl;
    typedef std::pair<uint32_t, uint32_t> TEdge;
    std::set<TEdge> es;
    links(c, mp, es);

    // Compute connected components
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Connected components using " << es.size() << " edges." << std::endl;
    concomp(c, mp, es);

    // Data out
    segmentOut(c, mp);
    
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
      ("quality,q", boost::program_options::value<uint16_t>(&c.minSeqQual)->default_value(10), "min. sequence quality")
      ("offset,m", boost::program_options::value<uint32_t>(&c.maxOffset)->default_value(5000), "max. basepair offset to connect read segments")
      ("segment,s", boost::program_options::value<uint32_t>(&c.segdev)->default_value(100), "segment breakpoint deviation")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
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
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
      std::cout << "Usage: lorax " << argv[0] << " [OPTIONS] -g <ref.fa> <tumor.bam>" << std::endl;
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
