#ifndef AMPLICON_H
#define AMPLICON_H

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

#include "bed.h"
#include "var.h"
#include "util.h"

namespace lorax
{

  struct AmpliconConfig {
    uint16_t minQual;
    std::string sample;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
    boost::filesystem::path tumor;
    boost::filesystem::path bedfile;
    boost::filesystem::path vcffile;
  };


  template<typename TConfig>
  inline int32_t
  selectReads(TConfig const& c, std::set<std::size_t> const& candidates) {
    // Open file handles
    samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.tumor.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Load bcf file
    htsFile* ibcffile = bcf_open(c.vcffile.string().c_str(), "r");
    hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
    bcf_hdr_t* bcfhdr = bcf_hdr_read(ibcffile);

    // Parse amplicons
    typedef boost::icl::interval_set<uint32_t> TChrIntervals;
    typedef std::vector<TChrIntervals> TRegionsGenome;
    TRegionsGenome scanRegions;
    if (!_parseBedIntervals(c.bedfile.string(), hdr, scanRegions)) {
      std::cerr << "Warning: Couldn't parse BED intervals. Do the chromosome names match?" << std::endl;
      return 1;
    }
    
    // Assign reads to SNPs
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      std::string chrName(hdr->target_name[refIndex]);

      // Load het. markers
      typedef std::vector<BiallelicVariant> TPhasedVariants;
      TPhasedVariants pv;
      if (!_loadVariants(ibcffile, bcfidx, bcfhdr, c.sample, chrName, pv)) continue;
      if (pv.empty()) continue;
      
      // Sort variants
      std::sort(pv.begin(), pv.end(), SortVariants<BiallelicVariant>());

      // Load reference
      int32_t seqlen = -1;
      char* seq = NULL;
      seq = faidx_fetch_seq(fai, chrName.c_str(), 0, hdr->target_len[refIndex], &seqlen);    
    
      // Assign reads to haplotypes
      std::set<std::size_t> h1;
      std::set<std::size_t> h2;
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	if ((rec->core.qual < c.minQual) || (rec->core.tid<0)) continue;
	std::size_t seed = hash_string(bam_get_qname(rec));
	
	uint32_t hp1votes = 0;
	uint32_t hp2votes = 0;
	TPhasedVariants::const_iterator vIt = std::lower_bound(pv.begin(), pv.end(), BiallelicVariant(rec->core.pos), SortVariants<BiallelicVariant>());
	TPhasedVariants::const_iterator vItEnd = std::upper_bound(pv.begin(), pv.end(), BiallelicVariant(lastAlignedPosition(rec)), SortVariants<BiallelicVariant>());

	// Get read sequence
	std::string sequence;
	sequence.resize(rec->core.l_qseq);
	uint8_t* seqptr = bam_get_seq(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	  
	// Parse CIGAR
	uint32_t* cigar = bam_get_cigar(rec);
	for(;vIt != vItEnd; ++vIt) {
	  int32_t gp = rec->core.pos; // Genomic position
	  int32_t sp = 0; // Sequence position
	  bool varFound = false;
	  for (std::size_t i = 0; ((i < rec->core.n_cigar) && (!varFound)); ++i) {
	    if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) gp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
		//Nop
	    }
	    else if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      if (gp + (int32_t) bam_cigar_oplen(cigar[i]) < vIt->pos) {
		gp += bam_cigar_oplen(cigar[i]);
		sp += bam_cigar_oplen(cigar[i]);
	      } else {
		for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp) {
		  if (gp == vIt->pos) {
		    varFound = true;
		    // Check REF allele
		    if (vIt->ref == std::string(seq + gp, seq + gp + vIt->ref.size())) {
		      // Check ALT allele
		      if ((sp + vIt->alt.size() < sequence.size()) && (sp + vIt->ref.size() < sequence.size())) {
			if (vIt->ref.size() == vIt->alt.size()) {
			  // SNP
			  if ((sequence.substr(sp, vIt->alt.size()) == vIt->alt) && (sequence.substr(sp, vIt->ref.size()) != vIt->ref)) {
			    // ALT supporting read
			    if (vIt->hap) ++hp1votes;
			    else ++hp2votes;
			  } else if ((sequence.substr(sp, vIt->alt.size()) != vIt->alt) && (sequence.substr(sp, vIt->ref.size()) == vIt->ref)) {
			    // REF supporting read
			    if (vIt->hap) ++hp2votes;
			    else ++hp1votes;
			  }
			} else if (vIt->ref.size() < vIt->alt.size()) {
			  // Insertion
			  int32_t diff = vIt->alt.size() - vIt->ref.size();
			  std::string refProbe = vIt->ref + std::string(seq + gp + vIt->ref.size(), seq + gp + vIt->ref.size() + diff);
			  if ((sequence.substr(sp, vIt->alt.size()) == vIt->alt) && (sequence.substr(sp, vIt->alt.size()) != refProbe)) {
			    // ALT supporting read
			    if (vIt->hap) ++hp1votes;
			    else ++hp2votes;
			  } else if ((sequence.substr(sp, vIt->alt.size()) != vIt->alt) && (sequence.substr(sp, vIt->alt.size()) == refProbe)) {
			    // REF supporting read
			    if (vIt->hap) ++hp2votes;
			    else ++hp1votes;
			  }
			} else {
			  // Deletion
			  int32_t diff = vIt->ref.size() - vIt->alt.size();
			  std::string altProbe = vIt->alt + std::string(seq + gp + vIt->ref.size(), seq + gp + vIt->ref.size() + diff);
			  if ((sequence.substr(sp, vIt->ref.size()) == altProbe) && (sequence.substr(sp, vIt->ref.size()) != vIt->ref)) {
			    // ALT supporting read
			    if (vIt->hap) ++hp1votes;
			    else ++hp2votes;
			  } else if ((sequence.substr(sp, vIt->ref.size()) != altProbe) && (sequence.substr(sp, vIt->ref.size()) == vIt->ref)) {
			    // REF supporting read
			    if (vIt->hap) ++hp2votes;
			    else ++hp1votes;
			  }
			}
		      }
		    }
		  }
		}
	      }
	    } else {
	      std::cerr << "Unknown Cigar options" << std::endl;
	      return 1;
	    }
	  }
	  int32_t hp = 0;
	  if (hp1votes > 2*hp2votes) hp = 1;
	  else if (hp2votes > 2*hp1votes) hp = 2;
	  if (hp) {
	    if (hp == 1) h1.insert(seed);
	    else if (hp == 2) h2.insert(seed);
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
      if (seq != NULL) free(seq);
    }
    fai_destroy(fai);

    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    return candidates.size();
  }
    

  template<typename TConfig>
  inline int32_t
  runAmplicon(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Select amplicon reads
    std::set<std::size_t> candidates;
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Selecting reads." << std::endl;
    int32_t candsize = selectReads(c, candidates);
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    return 0;
  }


  
  int amplicon(int argc, char** argv) {
    AmpliconConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Options");
    generic.add_options()
      ("help,?", "show help message")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minQual)->default_value(10), "min. sequence quality")
      ("sample,s", boost::program_options::value<std::string>(&c.sample)->default_value("NA12878"), "sample name (as in VCF/BCF file)")
      ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file")
      ("bedfile,b", boost::program_options::value<boost::filesystem::path>(&c.bedfile), "amplicon regions in BED format")
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
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome")) || (!vm.count("vcffile"))) {
      std::cout << "Usage: lorax " << argv[0] << " [OPTIONS] -g <ref.fa> -s NA12878 -v <snps.bcf> <tumor.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Check input BAM file
    if (vm.count("input-file")) {
      if (!(boost::filesystem::exists(c.tumor) && boost::filesystem::is_regular_file(c.tumor) && boost::filesystem::file_size(c.tumor))) {
	std::cerr << "Input BAM file is missing: " << c.tumor.string() << std::endl;
	return 1;
      }
      samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.tumor.string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.tumor.string().c_str());
      if (idx == NULL) {
	std::cerr << "Fail to open index for " << c.tumor.string() << std::endl;
	return 1;
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.tumor.string() << std::endl;
	return 1;
      }
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }
  
    // Check VCF/BCF file
    if (vm.count("vcffile")) {
      if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
	std::cerr << "Input VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
	return 1;
      }
      htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
      if (ifile == NULL) {
	std::cerr << "Fail to open file " << c.vcffile.string() << std::endl;
	return 1;
      }
      hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
      if (bcfidx == NULL) {
	std::cerr << "Fail to open index file for " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_t* hdr = bcf_hdr_read(ifile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_destroy(hdr);
      hts_idx_destroy(bcfidx);
      bcf_close(ifile);
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runAmplicon(c);
  }

}

#endif
