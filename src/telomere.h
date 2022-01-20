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
    uint16_t minMapQual;
    uint16_t minClip;
    uint16_t minSplit;
    uint32_t minSegmentSize;
    uint32_t maxSegmentSize;
    uint32_t minChrLen;
    float contam;
    boost::filesystem::path genome;
    boost::filesystem::path outfile;
    boost::filesystem::path tumor;
    boost::filesystem::path control;
  };


  template<typename TConfig>
  inline int32_t
  runTelomere(TConfig& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Open file handles
    samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.tumor.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    samFile* cfile = sam_open(c.control.string().c_str(), "r");
    hts_set_fai_filename(cfile, c.genome.string().c_str());
    hts_idx_t* cidx = sam_index_load(cfile, c.control.string().c_str());
    faidx_t* fai = fai_load(c.genome.string().c_str());
    
    // Parse genome, process chromosome by chromosome
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      // Any data?
      if ((!mappedReads(idx, refIndex, c.tumor.string())) || (!mappedReads(idx, refIndex, c.control.string()))) continue;

      // Large enough chromosome?
      if (hdr->target_len[refIndex] > c.minChrLen) {
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();	  
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Parsing " << hdr->target_name[refIndex] << std::endl;
      } else continue;
    }
    
    // Clean-up
    bam_hdr_destroy(hdr);
    fai_destroy(fai);
    hts_idx_destroy(idx);
    sam_close(samfile);
    hts_idx_destroy(cidx);
    sam_close(cfile);
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    return 0;
  }


  
  int telomere(int argc, char** argv) {
    TelomereConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Options");
    generic.add_options()
      ("help,?", "show help message")
      ("qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("clip,c", boost::program_options::value<uint16_t>(&c.minClip)->default_value(25), "min. clipping length")
      ("split,s", boost::program_options::value<uint16_t>(&c.minSplit)->default_value(3), "min. split-read support")
      ("chrlen,l", boost::program_options::value<uint32_t>(&c.minChrLen)->default_value(40000000), "min. chromosome length")
      ("minsize,i", boost::program_options::value<uint32_t>(&c.minSegmentSize)->default_value(100), "min. segment size")
      ("maxsize,j", boost::program_options::value<uint32_t>(&c.maxSegmentSize)->default_value(10000), "max. segment size")
      ("contam,n", boost::program_options::value<float>(&c.contam)->default_value(0), "max. fractional tumor-in-normal contamination")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("matched,m", boost::program_options::value<boost::filesystem::path>(&c.control), "matched control BAM")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.bed"), "BED output file")
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
