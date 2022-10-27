#ifndef PCT_H
#define PCT_H

#define BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX

#include <fstream>
#include <iomanip>

#include <boost/dynamic_bitset.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
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

#include "gfa.h"
#include "gaf.h"

namespace lorax
{

  struct PctConfig {
    bool gfaMode;
    bool hasOutfile;
    bool hasFastq;
    float pct;
    uint32_t len;
    uint32_t delcut;
    boost::filesystem::path genome;
    boost::filesystem::path gfafile;
    boost::filesystem::path seqfile;
    boost::filesystem::path sample;
    boost::filesystem::path outfile;
    boost::filesystem::path readsfile;
    boost::filesystem::path outfastq;
  };


  template<typename TConfig>
  inline void
  parseReads(TConfig const& c, Graph const& g, std::vector<AlignRecord> const& aln) {
    typedef std::vector<AlignRecord> TAlignRecords;
    
    faidx_t* fai = fai_load(c.readsfile.string().c_str());
    for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
      std::string tname(faidx_iseq(fai, refIndex));
      int32_t seqlen = 0;
      char* seq = faidx_fetch_seq(fai, tname.c_str(), 0, faidx_seq_len(fai, tname.c_str()), &seqlen);
      std::size_t seed = hash_lr(tname);

      // Find alignment records
      typename TAlignRecords::const_iterator iter = std::lower_bound(aln.begin(), aln.end(), AlignRecord(0, seed), SortAlignRecord<AlignRecord>());
      for(; ((iter != aln.end()) && (iter->seed == seed)); ++iter) {
	std::cerr << seed << ',' << tname << std::endl;
      }
      
      free(seq);
    }
    fai_destroy(fai);
    
    /*
    	// Evaluate alignment record

	// Parse only primary alignments
	//if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;

	// Unmapped read
	//if (rec->core.flag & (BAM_FUNMAP)) {
	//if (c.hasOutfile) ofile << bam_get_qname(rec) << '\t' << sequence.size() << "\t0\t0\tunmapped" << std::endl;
	
	// Parse cigar
	uint32_t largestdel = 0;
	//uint32_t rp = rec->core.pos; // reference pointer
	uint32_t rp = 0;
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
	for (uint32_t i = 0; i < ar.cigarop.size(); ++i) {
	  if (ar.cigarop[i] == BAM_CEQUAL) {
	    match += ar.cigarlen[i];
	    sp += ar.cigarlen[i];
	    rp += ar.cigarlen[i];
	  }
	  else if (ar.cigarop[i] == BAM_CDIFF) {
	    mismatch += ar.cigarlen[i];
	    sp += ar.cigarlen[i];
	    rp += ar.cigarlen[i];
	  }
	  else if (ar.cigarop[i] == BAM_CDEL) {
	    if (ar.cigarlen[i] > largestdel) largestdel = ar.cigarlen[i];
	    ++del;
	    delsize += ar.cigarlen[i];
	    rp += ar.cigarlen[i];
	  }
	  else if (ar.cigarop[i] == BAM_CINS) {
	    ++ins;
	    inssize += ar.cigarlen[i];
	    sp += ar.cigarlen[i];
	  }
	  else if (ar.cigarop[i] == BAM_CSOFT_CLIP) {
	    ++sc;
	    scsize += ar.cigarlen[i];
	    sp += ar.cigarlen[i];
	  }
	  else if (ar.cigarop[i] == BAM_CREF_SKIP) {
	    rp += ar.cigarlen[i];
	  }
	  else if (ar.cigarop[i] == BAM_CHARD_CLIP) {
	    ++hc;
	    hcsize += ar.cigarlen[i];
	    sp += ar.cigarlen[i];
	  }
	  else {
	    std::cerr << "Warning: Unknown Cigar option " << ar.cigarop[i] << std::endl;
	    exit(-1);
	  }
	}

	// Percent identity
	double pctval = (double) (match) / (double) ar.qlen;
	if (c.hasOutfile) ofile << ar.qname << '\t' << ar.qlen << '\t' << pctval << '\t' << largestdel << "\taligned\t" << match << '\t' << mismatch << '\t' << del << '\t' << delsize << '\t' << ins << '\t' << inssize << '\t' << sc << '\t' << scsize << '\t' << hc << '\t' << hcsize << std::endl;
      }
  */
  }
    
  template<typename TConfig>
  inline void
  pctBam(TConfig const& c) {
    // Open file handles
    samFile* samfile = sam_open(c.sample.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Output file
    std::ofstream ofile;
    if (c.hasOutfile) {
      ofile.open(c.outfile.string().c_str(), std::ofstream::out | std::ofstream::trunc);
      ofile << "qname\treadlen\tpctidentity\tlargestdel\tmapped\tmatches\tmismatches\tdeletions\tdelsize\tinsertions\tinssize\tsoftclips\tsoftclipsize\thardclips\thardclipsize" << std::endl;
    }

    // Output FASTQ file
    boost::iostreams::filtering_ostream dataOut;
    if (c.hasFastq) {
      dataOut.push(boost::iostreams::gzip_compressor());
      dataOut.push(boost::iostreams::file_sink(c.outfastq.string().c_str(), std::ios_base::out | std::ios_base::binary));
    }
    
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
      uint8_t* seqptr = bam_get_seq(rec);
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

      // Get read qualities
      typedef std::vector<uint8_t> TQuality;
      TQuality quality;
      quality.resize(rec->core.l_qseq, (uint8_t) 33);
      uint8_t* qualptr = bam_get_qual(rec);
      if (!(qualptr[0] == 0xff)) {
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
      }

      // Unmapped read
      if (rec->core.flag & (BAM_FUNMAP)) {
	if (c.hasOutfile) ofile << bam_get_qname(rec) << '\t' << sequence.size() << "\t0\t0\tunmapped" << std::endl;

	// Output FASTQ record
	if (c.hasFastq) {
	  dataOut << "@" << bam_get_qname(rec) << std::endl;
	  dataOut << sequence << std::endl;
	  dataOut << "+" << std::endl;
	  for(uint32_t i = 0; i < quality.size(); ++i) {
	    char c = 33 + quality[i];
	    dataOut << c;
	  }
	  dataOut << std::endl;
	}
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
	    if (sequence[sp] == upper(seq[rp])) ++match;
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
	if (c.hasFastq) {
	  dataOut << "@" << bam_get_qname(rec) << std::endl;
	  dataOut << sequence << std::endl;
	  dataOut << "+" << std::endl;
	  for(uint32_t i = 0; i < quality.size(); ++i) {
	    char c = 33 + quality[i];
	    dataOut << c;
	  }
	  dataOut << std::endl;
	}
      }
    }
    if (seq != NULL) free(seq);
    if (c.hasOutfile) ofile.close();
    if (c.hasFastq) {
      dataOut.pop();
      dataOut.pop();
    }
    
    // Clean-up
    bam_destroy1(rec);
    fai_destroy(fai);
    bam_hdr_destroy(hdr);
    sam_close(samfile);
  }
    

  template<typename TConfig>
  inline int32_t
  runPct(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    if (c.gfaMode) {
      // Load pan-genome graph
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Load pan-genome graph" << std::endl;
      Graph g;
      parseGfa(c, g);

      // Parse alignments
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse alignments" << std::endl;
      std::vector<AlignRecord> aln;
      parseGaf(c, g, aln);

      // Parse reads
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse reads" << std::endl;
      parseReads(c, g, aln);
      exit(-1);
      


      
      // Write pan-genome graph
      //writeGfa(c, g);
    } else {
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse alignments" << std::endl;
      pctBam(c);
    }
    

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
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("graph,g", boost::program_options::value<boost::filesystem::path>(&c.gfafile), "GFA pan-genome graph")
      ("reads,e", boost::program_options::value<boost::filesystem::path>(&c.readsfile), "reads in FASTA format")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output tsv file")
      ;

    boost::program_options::options_description rdsopt("Read selection");
    rdsopt.add_options()
      ("pct,p", boost::program_options::value<float>(&c.pct)->default_value(0.95), "max. percent identity")
      ("len,l", boost::program_options::value<uint32_t>(&c.len)->default_value(5000), "min. read length")
      ("del,d", boost::program_options::value<uint32_t>(&c.delcut)->default_value(50), "min. deletion size")
      ("outfastq,f", boost::program_options::value<boost::filesystem::path>(&c.outfastq), "gzipped output fastq file")
      ;
      
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.sample), "input file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(rdsopt).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(rdsopt);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || ((!vm.count("reference")) && (!vm.count("graph")))) {
      std::cout << "Usage:" << std::endl;
      std::cout << "Linear reference genome: lorax " << argv[0] << " [OPTIONS] -r <ref.fa> <sample.bam>" << std::endl;
      std::cout << "Pan-genome graph: lorax " << argv[0] << " [OPTIONS] -g <pangenome.hg38.gfa.gz> -e <reads.fasta> <sample.gaf>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    if (vm.count("outfile")) c.hasOutfile = true;
    else c.hasOutfile = false;
    if (vm.count("outfastq")) c.hasFastq = true;
    else c.hasFastq = false;
    if (vm.count("graph")) {
      c.gfaMode = true;
      if (!vm.count("reads")) {
	std::cerr << "Reads are required for pan-genome graphs!" << std::endl;
	return -1;
      }

      // Random name for temporary file
      boost::uuids::uuid uuid = boost::uuids::random_generator()();
      c.seqfile = "tmp." + boost::lexical_cast<std::string>(uuid) + ".fa";
    } else c.gfaMode = false;

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
