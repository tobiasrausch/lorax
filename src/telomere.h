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


#include "junction.h"

namespace lorax
{

  struct TelomereConfig {
    typedef std::set<std::string> TStringSet;

    uint16_t minMapQual;
    uint32_t minClip;
    uint32_t minSupport;
    uint32_t period;
    uint32_t replen;
    uint32_t delta;
    uint32_t maxTelLen;
    uint32_t minChrLen;

    std::string outprefix;    
    std::string repeatStr;
    TStringSet fwdmotif;
    TStringSet revmotif;
    boost::filesystem::path genome;
    boost::filesystem::path sample;
  };


  struct TelomereRecord {
    uint32_t telfwd;
    uint32_t telrev;
    std::string qname;
    std::string sequence;

    TelomereRecord() {}
    
    TelomereRecord(uint32_t const tf, uint32_t const tr, std::string const& qn, std::string const& s) : telfwd(tf), telrev(tr), qname(qn), sequence(s) {}

    TelomereRecord(TelomereRecord const& tr) {
      telfwd = tr.telfwd;
      telrev = tr.telrev;
      qname = tr.qname;
      sequence = tr.sequence;
    }
  };
  
  template<typename TConfig>
  inline void
  collectCandidates(TConfig const& c, std::set<std::size_t> const& reads, std::map<std::size_t, TelomereRecord>& candidates) {
    uint32_t motiflen = c.period * c.replen;
    
    // Open file handles
    samFile* samfile = sam_open(c.sample.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.sample.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse BAM
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[refIndex] << std::endl;

      // Iterate alignments 
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	std::size_t seed = hash_lr(rec);
	if (reads.find(seed) != reads.end()) {

	  // Load sequence
	  std::string sequence;
	  sequence.resize(rec->core.l_qseq, 'N');
	  uint8_t* seqptr = bam_get_seq(rec);
	  for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	  if (rec->core.flag & BAM_FREVERSE) reverseComplement(sequence);  // Junctions are always reported on the forward strand

	  // Bitsets
	  typedef boost::dynamic_bitset<> TBitSet;
	  TBitSet telfwd;
	  TBitSet telrev;
	  telfwd.resize(rec->core.l_qseq, 0);
	  telrev.resize(rec->core.l_qseq, 0);
	
	  // Mark telomere hits
	  for(std::set<std::string>::const_iterator it = c.fwdmotif.begin(); it != c.fwdmotif.end(); ++it) {
	    std::size_t index = 0;
	    while ((index = sequence.find(*it, index)) != std::string::npos) telfwd[index++] = 1;
	  }
	  for(std::set<std::string>::const_iterator it = c.revmotif.begin(); it != c.revmotif.end(); ++it) {
	    std::size_t index = 0;
	    while ((index = sequence.find(*it, index)) != std::string::npos) telrev[index++] = 1;
	  }

	  // Candidate telomere read
	  uint32_t tf = telfwd.count();
	  uint32_t tr = telrev.count();
	  if (tf + tr >= motiflen) candidates.insert(std::make_pair(seed, TelomereRecord(tf, tr, std::string(bam_get_qname(rec)), sequence)));
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


  template<typename TConfig, typename TTelomereReads>
  inline void
  outputTelomereSVs(TConfig const& c, std::vector<Junction>& junctions, TTelomereReads& telreads) {
    // Open file handles
    samFile* samfile = sam_open(c.sample.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Find representative junction
    std::vector<uint32_t> jctidx;
    uint32_t idx = 0;
    while(idx < junctions.size()) {
      uint32_t idxEnd = idx;
      if (junctions[idx].support >= c.minSupport) {
	// Extend
	uint32_t bestIdxStart = idx;
	uint32_t bestIdxEnd = idx;
	for(uint32_t j = idx + 1; ((j < junctions.size()) && (junctions[j].refidx == junctions[idx].refidx) && (std::abs(junctions[j].refpos - junctions[idx].refpos) < 2 * c.delta)); ++j) {
	  if (junctions[j].support > junctions[bestIdxStart].support) bestIdxStart = j;
	  if (junctions[j].support == junctions[bestIdxStart].support) bestIdxEnd = j;
	  idxEnd = j;
	}
	//std::cerr << "(" << idx << ',' << idxEnd << ") (" << bestIdxStart << ',' << bestIdxEnd << ")" << std::endl;

	// Midpoint
	jctidx.push_back(bestIdxStart + ((bestIdxEnd - bestIdxStart) / 2));
      }
      idx = idxEnd + 1;
    }

    // Output telomere-associated SVs
    std::string filename = c.outprefix + ".svs.tsv";
    std::ofstream rfile(filename.c_str());
    rfile << "chr\tpos\tid\tsupport\tqual\ttelfwd\ttelrev" << std::endl;
    for(uint32_t compId = 0; compId < jctidx.size(); ++compId) {
      std::string id("TEL");
      std::string padNumber = boost::lexical_cast<std::string>(compId + 1);
      padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
      id += padNumber;

      // FASTA file of all junction supporting reads
      uint32_t supp = 0;
      uint32_t qsup = 0;
      uint32_t telfwd = 0;
      uint32_t telrev = 0;
      std::string fname = c.outprefix + "." + id  + ".fa";
      std::ofstream ffile(fname.c_str());
      for(uint32_t i = 0; i < junctions.size(); ++i) {
	if ((junctions[i].refidx == junctions[jctidx[compId]].refidx) && (std::abs(junctions[i].refpos - junctions[jctidx[compId]].refpos) < c.delta)) {
	  ffile << ">" << telreads[junctions[i].seed].qname << " " << hdr->target_name[junctions[i].refidx] << ":" << junctions[i].refpos << " forward:" << (int) junctions[i].forward << " scleft:" << (int) junctions[i].scleft << " telfwd:" << telreads[junctions[i].seed].telfwd << " telrev:" << telreads[junctions[i].seed].telrev << " seqlen:" << telreads[junctions[i].seed].sequence.size() << " seqpos:" << junctions[i].seqpos << " qual:" << junctions[i].qual << std::endl;
	  ffile << telreads[junctions[i].seed].sequence << std::endl;
	  ++supp;
	  qsup += junctions[i].qual;
	  if (junctions[i].forward) {
	    telfwd += telreads[junctions[i].seed].telfwd;
	    telrev += telreads[junctions[i].seed].telrev;
	  } else {
	    telfwd += telreads[junctions[i].seed].telrev;
	    telrev += telreads[junctions[i].seed].telfwd;
	  }
	}
      }
      ffile.close();
      // Average mapping quality
      qsup /= supp;

      // Summarize
      rfile << hdr->target_name[junctions[jctidx[compId]].refidx] << '\t' << junctions[jctidx[compId]].refpos << '\t' << id << '\t' << supp << '\t' << qsup << '\t' << telfwd << '\t' << telrev << std::endl;
    }
    rfile.close();

    // Clean-up
    bam_hdr_destroy(hdr);
    sam_close(samfile);
  }

  template<typename TConfig>
  inline int32_t
  runTelomere(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Telomere reads
    typedef std::map<std::size_t, TelomereRecord> TReadTelomere;
    TReadTelomere telreads;

    // Junctions
    typedef std::vector<Junction> TJunctionVector;
    TJunctionVector junctions;

    if (junctions.empty()) {
      // Parse SV breakpoints
      TJunctionVector readBp;
      findJunctions(c, readBp);

      // Cluster breakpoints
      clusterJunctions(c, readBp);

      // Junction reads
      std::set<std::size_t> reads;
      for(uint32_t i = 0; i < readBp.size(); ++i) {
	if (readBp[i].support >= c.minSupport) reads.insert(readBp[i].seed);
      }

      // Telomere reads
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Telomere content screening for " << reads.size() << " reads." << std::endl;
      collectCandidates(c, reads, telreads);

      // Final junctions
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Intersection of telomere and SV reads" << std::endl;
      for(uint32_t i = 0; i < readBp.size(); ++i) {
	if (telreads.find(readBp[i].seed) != telreads.end()) {
	  readBp[i].support = 1;
	  junctions.push_back(readBp[i]);
	}
      }
    }

    // Cluster telomere associated SV junctions
    clusterJunctions(c, junctions);

    // Data out
    outputTelomereSVs(c, junctions, telreads);

#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  
  int telomere(int argc, char** argv) {
    TelomereConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Options");
    generic.add_options()
      ("help,?", "show help message")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(10), "min. mapping quality")
      ("minclip,c", boost::program_options::value<uint32_t>(&c.minClip)->default_value(25), "min. clipping length")
      ("support,s", boost::program_options::value<uint32_t>(&c.minSupport)->default_value(3), "min. telomere read support")
      ("delta,d", boost::program_options::value<uint32_t>(&c.delta)->default_value(250), "max. breakpoint junction offset")
      ("tellen,t", boost::program_options::value<uint32_t>(&c.maxTelLen)->default_value(20000), "max. telomere length")
      ("chrlen,l", boost::program_options::value<uint32_t>(&c.minChrLen)->default_value(40000000), "min. chromosome length")
      ("repeats,r", boost::program_options::value<std::string>(&c.repeatStr)->default_value("TTAGGG,TCAGGG,TGAGGG,TTGGGG"), "repeat units")
      ("period,p", boost::program_options::value<uint32_t>(&c.period)->default_value(3), "repeat period")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("out"), "output prefix")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.sample), "input file")
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

    // Create telomere motifs
    createRepeatMotifs(c, true);

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
