#ifndef PCT_H
#define PCT_H

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

  struct PctConfig {
    float pct;
    uint32_t len;
    boost::filesystem::path genome;
    boost::filesystem::path sample;
    boost::filesystem::path outfile;
  };

  struct PctRead {
  };

  template<typename TConfig>
  inline void
  pctscreen(TConfig const& c, std::map<std::size_t, PctRead>& reads) {
    typedef typename std::map<std::size_t, PctRead> TReadMap;
    
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
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY)) continue;
	std::size_t seed = hash_lr(rec);


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
	bool dir = true;
	if (rec->core.flag & BAM_FREVERSE) {
	  dir = false;
	  int32_t seqTmp = seqStart;
	  seqStart = sp - seqEnd;
	  seqEnd = sp - seqTmp;
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
  runPct(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Screen mappings" << std::endl;
    std::map<std::size_t, PctRead> reads;
    pctscreen(c, reads);
    

#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  
  int pct(int argc, char** argv) {
    PctConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Options");
    generic.add_options()
      ("help,?", "show help message")
      ("pct,p", boost::program_options::value<float>(&c.pct)->default_value(0.95), "percent identity threshold")
      ("len,l", boost::program_options::value<uint32_t>(&c.len)->default_value(5000), "read length threshold")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.fa"), "output fasta file")
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
      std::cout << "Usage: lorax " << argv[0] << " [OPTIONS] -g <ref.fa> <sample.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runPct(c);
  }

}

#endif
