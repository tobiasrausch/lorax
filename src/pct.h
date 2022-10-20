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
    bool hasOutfile;
    float pct;
    uint32_t len;
    uint32_t delcut;
    boost::filesystem::path genome;
    boost::filesystem::path sample;
    boost::filesystem::path outfile;
    boost::filesystem::path outfastq;
  };

  template<typename TConfig>
  inline void
  pctscreen(TConfig const& c) {
    // Open file handles
    samFile* samfile = sam_open(c.sample.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.sample.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Output file
    std::ofstream ofile;
    if (c.hasOutfile) {
      ofile.open(c.outfile.string().c_str(), std::ofstream::out | std::ofstream::trunc);
      ofile << "qname\treadlen\tpctidentity\tlargestdel\tmapped\tmatches\tmismatches\tdeletions\tdelsize\tinsertions\tinssize\tsoftclips\tsoftclipsize\thardclips\thardclipsize" << std::endl;
    }

    // Output FASTQ file
    std::ofstream ffile(c.outfastq.string().c_str());
    
    // Parse BAM
    int32_t refIndex = -1;
    char* seq = NULL;
    faidx_t* fai = fai_load(c.genome.string().c_str());
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      // New chromosome?
      if ((!(rec->core.flag & BAM_FUNMAP)) && (rec->core.tid != refIndex)) {
	// Load new chromosome
	if (seq != NULL) free(seq);
	refIndex = rec->core.tid;
	int32_t seqlen = -1;
	std::string tname(hdr->target_name[refIndex]);
	seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);

	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[refIndex] << std::endl;
      }
	
      // Parse only primary alignments
      if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;

      // Get the read sequence
      std::string sequence;
      sequence.resize(rec->core.l_qseq);
      typedef std::vector<uint8_t> TQuality;
      TQuality quality;
      quality.resize(rec->core.l_qseq);
      uint8_t* seqptr = bam_get_seq(rec);
      uint8_t* qualptr = bam_get_qual(rec);
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) {
	quality[i] = qualptr[i];
	sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
      }

      // Unmapped read
      if (rec->core.flag & (BAM_FUNMAP)) {
	if (c.hasOutfile) ofile << bam_get_qname(rec) << '\t' << sequence.size() << "\t0\t0\tunmapped" << std::endl;

	// Output FASTQ record
	ffile << "@" << bam_get_qname(rec) << std::endl;
	ffile << sequence << std::endl;
	ffile << "+" << std::endl;
	for(uint32_t i = 0; i < quality.size(); ++i) {
	  char c = 33 + quality[i];
	  ffile << c;
	}
	ffile << std::endl;
	continue;
      }
	
      // Parse cigar
      uint32_t largestdel = 0;
      uint32_t rp = rec->core.pos; // reference pointer
      uint32_t sp = 0; // sequence pointer
      uint32_t mismatch = 0;
      uint32_t match = 0;
      uint32_t del = 0;
      uint32_t delsize = 0;
      uint32_t ins = 0;
      uint32_t inssize = 0;
      uint32_t sc = 0;
      uint32_t scsize = 0;
      uint32_t hc = 0;
      uint32_t hcsize = 0;
      uint32_t* cigar = bam_get_cigar(rec);
      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	  for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
	    if (sequence[sp] == seq[rp]) ++match;
	    else ++mismatch;
	    ++sp;
	    ++rp;
	  }
	} else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	  if (bam_cigar_oplen(cigar[i]) > largestdel) largestdel = bam_cigar_oplen(cigar[i]);
	  ++del;
	  delsize += bam_cigar_oplen(cigar[i]);
	  rp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	  ++ins;
	  inssize += bam_cigar_oplen(cigar[i]);
	  sp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	  ++sc;
	  scsize += bam_cigar_oplen(cigar[i]);
	  sp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	  rp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	  ++hc;
	  hcsize += bam_cigar_oplen(cigar[i]);
	  sp += bam_cigar_oplen(cigar[i]);
	} else {
	  std::cerr << "Warning: Unknown Cigar options!" << std::endl;
	}
      }

      // Percent identity
      double pctval = (double) (match) / (double) sequence.size();
      if (c.hasOutfile) ofile << bam_get_qname(rec) << '\t' << sequence.size() << '\t' << pctval << '\t' << largestdel << "\taligned\t" << match << '\t' << mismatch << '\t' << del << '\t' << delsize << '\t' << ins << '\t' << inssize << '\t' << sc << '\t' << scsize << '\t' << hc << '\t' << hcsize << std::endl;


      // Selected read?
      if ((sequence.size() >= c.len) && ((pctval <= c.pct) || (largestdel >= c.delcut))) {
	// Output FASTQ record
	ffile << "@" << bam_get_qname(rec) << std::endl;
	ffile << sequence << std::endl;
	ffile << "+" << std::endl;
	for(uint32_t i = 0; i < quality.size(); ++i) {
	  char c = 33 + quality[i];
	  ffile << c;
	}
	ffile << std::endl;
      }
    }
    if (seq != NULL) free(seq);
    if (c.hasOutfile) ofile.close();
    ffile.close();
    
    // Clean-up
    bam_destroy1(rec);
    fai_destroy(fai);
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
    pctscreen(c);
    

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
      ("del,d", boost::program_options::value<uint32_t>(&c.delcut)->default_value(50), "min. deletion size")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output tsv file")
      ("outfastq,f", boost::program_options::value<boost::filesystem::path>(&c.outfastq)->default_value("out.fq"), "output fastq file")
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

    if (vm.count("outfile")) c.hasOutfile = true;
    else c.hasOutfile = false;

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
